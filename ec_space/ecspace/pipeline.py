"""Layered, multi-process, resumable enumeration with on-disk (hash-sharded) de-duplication.

out/
  closed/L{n}/s{k}.smi.gz      unique closed-shell neutral molecules with n heavy atoms (sorted)
  closed/L{n}/DONE.json        layer counts (written last -> resume marker)
  species/L{n}/part-*.tsv.gz   all species owned by parents with n heavy atoms (n <= derive_max_heavy)
  reactions/L{n}/part-*.tsv.gz elementary reactions owned by those parents
  derive/L{n}.json             per-layer species/reaction counts (exact, or estimated from a sample)
"""
from __future__ import annotations

import gzip
import json
import math
import multiprocessing as mp
import os
import shutil
import time
import zlib
from collections import Counter
from pathlib import Path

from . import chem

T0 = time.time()


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')} +{time.time() - T0:8.0f}s]", *a, flush=True)


def heavy(smi: str) -> int:
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi)
    return sum(1 for a in m.GetAtoms() if a.GetSymbol() in ("C", "O"))


def shard_of(s: str, S: int) -> int:
    return zlib.crc32(s.encode()) % S


def read_lines(path: Path):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                yield line


def layer_files(out: Path, n: int):
    return sorted((out / "closed" / f"L{n}").glob("s*.smi.gz"))


# ------------------------------------------------------------------ generation workers
def _gen_worker(wid, cfgd, n_next, S, tmpdir, q_in, q_out):
    F = chem.Filters(cfgd["allow_peroxide"], cfgd["strain_filter"], cfgd["max_C"], cfgd["max_O"])
    buf = [[] for _ in range(S)]
    nbuf = 0
    tmpdir = Path(tmpdir)

    def flush():
        nonlocal nbuf
        for k in range(S):
            if buf[k]:
                with open(tmpdir / f"w{wid:04d}.s{k:05d}.gz", "ab") as fh:
                    fh.write(gzip.compress(("\n".join(buf[k]) + "\n").encode(), 1))
                buf[k] = []
        nbuf = 0

    while True:
        task = q_in.get()
        if task is None:
            break
        local = set()
        for p in task:
            local.update(chem.children(p, F))
        for c in local:
            buf[shard_of(c, S)].append(c)
        nbuf += len(local)
        if nbuf > 400_000:
            flush()
        q_out.put((len(task), len(local)))
    flush()
    q_out.put(None)


def _dedupe_shard(args):
    tmpdir, k, dst, level = args
    seen = set()
    for f in Path(tmpdir).glob(f"w*.s{k:05d}.gz"):
        with gzip.open(f, "rt") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line:
                    seen.add(line)
    data = sorted(seen)
    with gzip.open(dst, "wt", compresslevel=level) as fh:
        if data:
            fh.write("\n".join(data) + "\n")
    return k, len(data)


def auto_shards(cfg, est):
    if cfg.shards != "auto":
        return int(cfg.shards)
    return int(min(8192, max(8, 2 ** math.ceil(math.log2(max(1, est / 2_000_000))))))


def enumerate_closed(cfg):
    out = cfg.out
    cfgd = {"allow_peroxide": cfg.allow_peroxide, "strain_filter": cfg.strain_filter, "max_C": cfg.max_C, "max_O": cfg.max_O}
    # layer 1 seeds
    L1 = out / "closed" / "L1"
    if not (L1 / "DONE.json").exists():
        L1.mkdir(parents=True, exist_ok=True)
        seeds = sorted(s for s, el in (("C", "C"), ("O", "O")) if (el == "C" and cfg.max_C >= 1) or (el == "O" and cfg.max_O >= 1))
        with gzip.open(L1 / "s00000.smi.gz", "wt") as fh:
            fh.write("\n".join(seeds) + "\n")
        (L1 / "DONE.json").write_text(json.dumps({"n": 1, "count": len(seeds), "shards": 1}))
    prev = 1
    for n in range(2, cfg.max_heavy + 1):
        Ln = out / "closed" / f"L{n}"
        if (Ln / "DONE.json").exists():
            log(f"L{n}: already done ({json.loads((Ln / 'DONE.json').read_text())['count']:,})")
            prev = n
            continue
        info_prev = json.loads((out / "closed" / f"L{prev}" / "DONE.json").read_text())
        ratio = 10.0 if prev < 6 else 11.0
        S = auto_shards(cfg, info_prev["count"] * ratio * 3)
        tmp = out / "tmp" / f"L{n}"
        if tmp.exists():
            shutil.rmtree(tmp)
        tmp.mkdir(parents=True)
        Ln.mkdir(parents=True, exist_ok=True)
        for f in Ln.glob("*"):
            f.unlink()
        W = max(1, cfg.workers)
        ctx = mp.get_context("fork")
        q_in, q_out = ctx.Queue(maxsize=4 * W), ctx.Queue()
        procs = [ctx.Process(target=_gen_worker, args=(w, cfgd, n, S, str(tmp), q_in, q_out)) for w in range(W)]
        for p in procs:
            p.start()
        total_par = info_prev["count"]
        B = max(20, min(2000, total_par // (W * 20) or 20))
        log(f"L{n}: generating children of {total_par:,} parents  workers={W} shards={S} batch={B}")
        sent = done_par = gen = 0
        t1 = time.time()
        import threading

        stop = threading.Event()

        def drain():
            nonlocal done_par, gen
            finished = 0
            last = time.time()
            while finished < W:
                r = q_out.get()
                if r is None:
                    finished += 1
                    continue
                done_par += r[0]
                gen += r[1]
                if time.time() - last > 60:
                    last = time.time()
                    el = time.time() - t1
                    eta = el / max(done_par, 1) * (total_par - done_par)
                    log(f"  L{n}: {done_par:,}/{total_par:,} parents  {gen:,} children (with duplicates)  ETA {eta / 60:.1f} min")
            stop.set()

        th = threading.Thread(target=drain, daemon=True)
        th.start()
        batch = []
        for f in layer_files(out, prev):
            for s in read_lines(f):
                batch.append(s)
                if len(batch) >= B:
                    q_in.put(batch)
                    sent += len(batch)
                    batch = []
        if batch:
            q_in.put(batch)
        for _ in procs:
            q_in.put(None)
        th.join()
        for p in procs:
            p.join()
            if p.exitcode != 0:
                raise RuntimeError(f"worker exited with {p.exitcode}")
        log(f"L{n}: de-duplicating {S} shards")
        with ctx.Pool(W) as pool:
            res = pool.map(_dedupe_shard, [(str(tmp), k, str(Ln / f"s{k:05d}.smi.gz"), cfg.gzip_level) for k in range(S)])
        cnt = sum(c for _, c in res)
        shutil.rmtree(tmp)
        (Ln / "DONE.json").write_text(json.dumps({"n": n, "count": cnt, "shards": S, "parents": total_par,
                                                  "children_with_duplicates": gen, "seconds": round(time.time() - t1, 1)}))
        log(f"L{n}: {cnt:,} unique closed-shell molecules ({time.time() - t1:.0f}s)")
        prev = n
    tmp_root = out / "tmp"
    if tmp_root.exists() and not any(tmp_root.iterdir()):
        tmp_root.rmdir()


# ------------------------------------------------------------------ derivation workers
_D = None
_CFG = None


def _init_derive(cfg):
    global _D, _CFG
    _CFG = cfg
    _D = chem.Deriver(cfg)


def _derive_task(args):
    n, tid, smis, write_dir_sp, write_dir_rx, write, scale = args
    sp_c, rx_c, hist = Counter(), Counter(), Counter()
    sp_lines, rx_lines = [], []
    for p in smis:
        if _CFG.max_H is not None:
            if chem.formula(chem.Chem.MolFromSmiles(p))[1] > _CFG.max_H:
                continue
        try:
            species, rx = _D.derive(p)
        except Exception as e:          # never lose a whole task on one odd molecule
            sp_c["_errors"] += 1
            continue
        for k, t in species.items():
            sp_c[t] += 1
            if write:
                q, mult, c, h, o, li = chem.info_of(k)
                sp_lines.append(f"{k}\t{t}\t{q}\t{mult}\t{c}\t{h}\t{o}\t{li}\t{c + o}\t{p}")
        for fam, r, pr in rx:
            rx_c[fam] += 1
            if write:
                rx_lines.append(f"{fam}\t{' + '.join(r)}\t{' + '.join(pr)}\t{n}")
    if write:
        level = _CFG.gzip_level
        with gzip.open(Path(write_dir_sp) / f"part-{tid:06d}.tsv.gz", "wt", compresslevel=level) as fh:
            if sp_lines:
                fh.write("\n".join(sp_lines) + "\n")
        with gzip.open(Path(write_dir_rx) / f"part-{tid:06d}.tsv.gz", "wt", compresslevel=level) as fh:
            if rx_lines:
                fh.write("\n".join(rx_lines) + "\n")
    return len(smis), sp_c, rx_c


SPECIES_HEADER = "key\ttype\tcharge\tmultiplicity\tC\tH\tO\tLi\theavy\tparent"
REACTION_HEADER = "family\treactants\tproducts\tparent_heavy"


def derive_all(cfg):
    out = cfg.out
    W = max(1, cfg.workers)
    for n in range(1, cfg.max_heavy + 1):
        dj = out / "derive" / f"L{n}.json"
        if dj.exists():
            log(f"derive L{n}: already done")
            continue
        write = n <= cfg.derive_max_heavy
        frac = 1.0 if write else max(1e-6, cfg.count_sample)
        sp_dir, rx_dir = out / "species" / f"L{n}", out / "reactions" / f"L{n}"
        if write:
            for d in (sp_dir, rx_dir):
                if d.exists():
                    shutil.rmtree(d)
                d.mkdir(parents=True)
        info = json.loads((out / "closed" / f"L{n}" / "DONE.json").read_text())
        total = info["count"]
        B = max(10, min(2000, total // (W * 8) or 10))
        thr = int(frac * 1_000_000)

        def tasks():
            tid = 0
            batch = []
            for f in layer_files(out, n):
                for s in read_lines(f):
                    if frac < 1.0 and zlib.crc32(s.encode()) % 1_000_000 >= thr:
                        continue
                    batch.append(s)
                    if len(batch) >= B:
                        yield (n, tid, batch, str(sp_dir), str(rx_dir), write, frac)
                        tid += 1
                        batch = []
            if n == 2:
                batch.extend(x for x in cfg.special if heavy(x) == 2)
            if batch:
                yield (n, tid, batch, str(sp_dir), str(rx_dir), write, frac)

        t1 = time.time()
        log(f"derive L{n}: {total:,} parents, {'full + write' if write else f'sample {frac:g} (counts only)'}")
        sp_c, rx_c = Counter(), Counter()
        done = 0
        last = time.time()
        ctx = mp.get_context("fork")
        with ctx.Pool(W, initializer=_init_derive, initargs=(cfg,)) as pool:
            for np_, a, b in pool.imap_unordered(_derive_task, tasks(), chunksize=1):
                done += np_
                sp_c.update(a)
                rx_c.update(b)
                if time.time() - last > 60:
                    last = time.time()
                    exp = total * frac
                    log(f"  derive L{n}: {done:,}/{exp:,.0f} parents  species {sum(sp_c.values()):,}  reactions {sum(rx_c.values()):,}")
        scale = 1.0 / frac
        res = {"n": n, "parents": total, "sampled_parents": done, "exact": write,
               "species": {k: round(v * scale) for k, v in sp_c.items()},
               "reactions": {k: round(v * scale) for k, v in rx_c.items()},
               "seconds": round(time.time() - t1, 1)}
        dj.parent.mkdir(parents=True, exist_ok=True)
        dj.write_text(json.dumps(res, indent=1))
        log(f"derive L{n}: species {sum(res['species'].values()):,}  reactions {sum(res['reactions'].values()):,}"
            f"{'' if write else ' (estimated)'}  ({time.time() - t1:.0f}s)")


def write_headers(cfg):
    (cfg.out / "species").mkdir(parents=True, exist_ok=True)
    (cfg.out / "reactions").mkdir(parents=True, exist_ok=True)
    (cfg.out / "species" / "HEADER.tsv").write_text(SPECIES_HEADER + "\n")
    (cfg.out / "reactions" / "HEADER.tsv").write_text(REACTION_HEADER + "\n")
    rows = [("[H]", 0, 2, 1, 0), ("[H+]", 1, 1, 1, 0), ("[H-]", -1, 1, 1, 0), ("[e-]", -1, 2, 0, 0), ("[Li+]", 1, 1, 0, 1)]
    pool = "\n".join(f"{k}\tpool\t{q}\t{m}\t0\t{h}\t0\t{li}\t0\t" for k, q, m, h, li in rows)
    (cfg.out / "species" / "POOL.tsv").write_text(SPECIES_HEADER + "\n" + pool + "\n")
