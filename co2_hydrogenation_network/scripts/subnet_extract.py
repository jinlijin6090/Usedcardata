"""Extract the sub-network between given reactants and products, within (minimum steps + k), and export it.

    python subnet_extract.py <db_dir> --from CO2 --to 对二甲苯 --extra 1 --out 结果_对二甲苯
    python subnet_extract.py ../data/net_C16_species --from "O=C=O" --to "Cc1ccc(C)cc1" "CCCCCCCCCCCCCCCC" --extra 2

<db_dir> is a species-level database written by species_merge.py (one node per species, reversible reactions).
Species can be given as SMILES (any valid writing; resonance writings work), as a label in nodes.tsv, or as one
of the common names in NAMES below.

Definitions (a reaction can run in either direction; each direction is a "step"):
  pool        always-available small species (default H, H2, H2O; --pool to change, --pool none for nothing)
  f(x)        minimum number of steps to make x from the reactants + pool: a step can run once all of its inputs
              are made; it ends at max(f(inputs)) + 1 (same definition as the path search of the web pages)
  D(t)        f(target t): the minimum number of steps to the target
  L(x)        the latest step by which x must be available so that the target is made by step D(t) + k
              (critical-path deadline, propagated back from the target through every input of every step)
  a step e    must fire by T(e) = max L(its useful outputs); it is kept for target t when f(e) <= T(e),
              f(e) = max f(inputs) + 1.  Its slack is T(e) - f(e) (0 = on the critical path of a fastest route).
With k = 0 this is every step of every fastest route, including the side branches that make a co-reactant in
time (e.g. CH3 for a methylation).  Several targets: union of the per-target sub-networks.

Outputs in --out (default ./subnet_<time>):
  物种.csv     node_id smiles label kind formula C H O 最早第几步生成f 最晚第几步需要L 角色
  反应.csv     reaction_id 方程(按使用方向) 反应族 第几步 余量 原始基元反应数 输入 输出
  subnet.json  the same as JSON;  subnet.graphml  bipartite species / reaction graph (Cytoscape, Gephi, networkx)
  网页.html    self-contained interactive view (double-click to open; columns = steps)
  说明.txt     the query, minimum steps per target and the counts
"""
from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

T0 = time.time()
INF = np.int32(1 << 29)
NAMES = {
    "CO2": "O=C=O", "二氧化碳": "O=C=O", "CO": "[C-]#[O+]", "一氧化碳": "[C-]#[O+]", "H2": "[H][H]", "氢气": "[H][H]",
    "H": "[H]", "H2O": "O", "水": "O", "甲醇": "CO", "methanol": "CO", "甲烷": "C", "甲醛": "C=O", "甲酸": "O=CO",
    "二甲醚": "COC", "DME": "COC", "甲基": "[CH3]", "乙烯": "C=C", "丙烯": "C=CC", "丁烯": "C=CCC", "乙烷": "CC", "丙烷": "CCC",
    "苯": "c1ccccc1", "甲苯": "Cc1ccccc1", "对二甲苯": "Cc1ccc(C)cc1", "间二甲苯": "Cc1cccc(C)c1", "邻二甲苯": "Cc1ccccc1C",
    "乙苯": "CCc1ccccc1", "萘": "c1ccc2ccccc2c1", "正十六烷": "CCCCCCCCCCCCCCCC", "正十二烷": "CCCCCCCCCCCC",
    "正辛烷": "CCCCCCCC", "异辛烷": "CC(C)CC(C)(C)C", "环己烷": "C1CCCCC1", "甲基环己烷": "CC1CCCCC1",
}


def log(*a):
    print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)


# ------------------------------------------------------------------ database
def load(db: Path):
    cache = db / "_cache.pkl"
    files = sorted(db.glob("reactions.part*.tsv.gz")) or [db / "reactions.tsv.gz"]
    stamp = [(f.name, f.stat().st_size) for f in files + [db / "nodes.tsv.gz"]]
    if cache.exists():
        try:
            with open(cache, "rb") as fh:
                d = pickle.load(fh)
            if d.get("stamp") == stamp:
                log(f"database from cache: {len(d['smiles']):,} species, {len(d['R']):,} reactions")
                return d
        except Exception:
            pass
    log("reading the database (first time only; a cache file is written next to it)")
    cols = {k: [] for k in ("smiles", "label", "kind", "formula", "C", "H", "O")}
    with gzip.open(db / "nodes.tsv.gz", "rt", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            cols["smiles"].append(p[1]); cols["label"].append(p[2]); cols["kind"].append(p[3]); cols["formula"].append(p[4])
            cols["C"].append(int(p[5])); cols["H"].append(int(p[6])); cols["O"].append(int(p[7]))
    rows, blocks = [], []
    fam_codes = {}
    n = 0
    for f in files:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            next(fh)
            for line in fh:
                p = line.split("\t", 5)
                r = [int(x) for x in p[1].split(",")]
                q = [int(x) for x in p[2].split(",")]
                rows.append(r + [-1] * (3 - len(r)) + q + [-1] * (3 - len(q)) + [fam_codes.setdefault(p[3], len(fam_codes)), int(p[4])])
                n += 1
                if len(rows) >= 1_000_000:
                    blocks.append(np.array(rows, dtype=np.int32))
                    rows = []
        log(f"  {f.name}: {n:,} reactions")
    if rows:
        blocks.append(np.array(rows, dtype=np.int32))
    A = np.concatenate(blocks) if blocks else np.zeros((0, 8), dtype=np.int32)
    del blocks, rows
    d = {"stamp": stamp, **cols, "C": np.array(cols["C"], dtype=np.int16), "H": np.array(cols["H"], dtype=np.int16),
         "O": np.array(cols["O"], dtype=np.int16), "R": np.ascontiguousarray(A[:, 0:3]), "P": np.ascontiguousarray(A[:, 3:6]),
         "fam": np.ascontiguousarray(A[:, 6]), "fam_names": [k for k, _ in sorted(fam_codes.items(), key=lambda x: x[1])],
         "nor": np.ascontiguousarray(A[:, 7])}
    try:
        with open(cache, "wb") as fh:
            pickle.dump(d, fh, protocol=5)
    except OSError:
        pass
    return d


def resolve(d, names, db):
    idx = {}
    for i, s in enumerate(d["smiles"]):
        idx[s] = i
    for i, s in enumerate(d["label"]):
        if s:
            idx.setdefault(s, i)
    out = []
    for n in names:
        s = NAMES.get(n, NAMES.get(n.upper(), n))
        i = idx.get(s)
        if i is None:
            try:
                from rdkit import Chem, RDLogger
                RDLogger.DisableLog("rdApp.*")
                m = Chem.MolFromSmiles(s)
                if m is not None:
                    c = Chem.MolToSmiles(m)
                    i = idx.get(c, idx.get(Chem.MolToSmiles(Chem.MolFromSmiles(c))))
                    if i is None:                       # a non-representative resonance writing
                        with gzip.open(db / "aliases.tsv.gz", "rt", encoding="utf-8") as fh:
                            for line in fh:
                                p = line.rstrip("\n").split("\t")
                                if p[4] == c and p[3]:
                                    i = int(p[3])
                                    break
            except ImportError:
                pass
        if i is None:
            sys.exit(f"找不到物种：{n}（可以用 SMILES、nodes.tsv 里的 label，或常用名：{'、'.join(list(NAMES)[:12])}…）")
        out.append(i)
    return out


# ------------------------------------------------------------------ search
def forward(IN, OUT, N, sources):
    f = np.full(N + 1, INF, dtype=np.int32)
    f[N] = 0                                    # padding column
    f[sources] = 0
    L = 0
    while True:
        fin = f[IN].max(axis=1)
        o = OUT[fin == L].ravel()
        o = o[o != N]
        o = o[f[o] > L + 1]
        f[o] = L + 1
        L += 1
        if not (f[:N] == L).any() or L > 1000:     # nothing new at this level: no step can fire later
            break
    fin = f[IN].max(axis=1)
    return f, np.where(fin < INF, fin + 1, INF).astype(np.int32)


def prepare(IN, OUT, N, fe, budget_max, blocked):
    """index the steps that can fire by budget_max once, for all deadline searches"""
    usable = np.flatnonzero(fe <= budget_max)
    INu, OUTu, feu = IN[usable], OUT[usable], fe[usable]
    ent_node = INu.ravel()
    ent_edge = np.repeat(np.arange(len(usable), dtype=np.int64), IN.shape[1])
    keep = (ent_node != N) & ~blocked[ent_node]
    ent_node, ent_edge = ent_node[keep], ent_edge[keep]
    o = np.argsort(ent_node, kind="stable")
    ent_node, ent_edge = ent_node[o], ent_edge[o]
    starts = np.flatnonzero(np.r_[True, ent_node[1:] != ent_node[:-1]]) if len(ent_node) else np.zeros(0, dtype=np.int64)
    omask = blocked.copy()
    omask[N] = True
    OUTm = np.where(omask[OUTu], N, OUTu)          # outputs that carry no deadline point at the padding slot
    return {"usable": usable, "OUT": OUTm, "fe": feu, "ent_edge": ent_edge, "starts": starts,
            "nodes_g": ent_node[starts], "E": len(IN), "N": N}


def deadlines(S, target, budget):
    """latest step L(x) by which x must be available so that the target is made by step `budget`.
    A step e must fire by T(e) = max L(useful outputs); it is feasible when fe(e) <= T(e); then each input must be
    available by T(e) - 1.  Pool and reactant species are always available and carry no deadline."""
    NEG = np.int32(-1)
    N = S["N"]
    L = np.full(N + 1, NEG, dtype=np.int32)
    L[target] = budget
    OUTm, feu, ent_edge, starts, nodes_g = S["OUT"], S["fe"], S["ent_edge"], S["starts"], S["nodes_g"]
    active = np.zeros(len(feu), dtype=bool)
    while len(starts):
        L[N] = NEG
        T = L[OUTm].max(axis=1)
        feas = (T >= 0) & (feu <= T)
        cand = np.where(feas, T - 1, NEG).astype(np.int32)
        best = np.maximum.reduceat(cand[ent_edge], starts)
        upd = best > L[nodes_g]
        if not upd.any():
            break
        L[nodes_g[upd]] = best[upd]
    L[N] = NEG
    T = L[OUTm].max(axis=1)
    Te = np.full(S["E"], NEG, dtype=np.int32)
    Te[S["usable"]] = T
    return L, Te


# ------------------------------------------------------------------ export
def write_all(out, d, sub, meta):
    out.mkdir(parents=True, exist_ok=True)
    nodes, edges = sub["nodes"], sub["edges"]
    with open(out / "物种.csv", "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["node_id", "smiles", "label", "kind", "formula", "C", "H", "O", "最早第几步生成f", "最晚第几步需要L", "角色"])
        for n in nodes:
            w.writerow([n["id"], n["smiles"], n["label"], n["kind"], n["formula"], n["C"], n["H"], n["O"], n["f"],
                        "" if n["L"] is None else n["L"], n["role"]])
    with open(out / "反应.csv", "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["reaction_id", "方程（按使用方向）", "反应族", "最早第几步", "额外步数", "时间余量", "原始基元反应数", "输入node_id", "输出node_id"])
        for e in edges:
            w.writerow([e["rid"], e["eq"], e["fam"], e["step"], e["extra"], e["slack"], e["n_original"],
                        ",".join(map(str, e["in"])), ",".join(map(str, e["out"]))])
    (out / "subnet.json").write_text(json.dumps({"meta": meta, "nodes": nodes, "edges": edges}, ensure_ascii=False), encoding="utf-8")
    g = ['<?xml version="1.0" encoding="UTF-8"?>',
         '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
         '<key id="type" for="node" attr.name="type" attr.type="string"/>',
         '<key id="smiles" for="node" attr.name="smiles" attr.type="string"/>',
         '<key id="label" for="node" attr.name="label" attr.type="string"/>',
         '<key id="kind" for="node" attr.name="kind" attr.type="string"/>',
         '<key id="step" for="node" attr.name="step" attr.type="int"/>',
         '<key id="role" for="node" attr.name="role" attr.type="string"/>',
         '<key id="family" for="node" attr.name="family" attr.type="string"/>',
         '<graph id="subnet" edgedefault="directed">']
    q = lambda s: html.escape(str(s), quote=True)  # noqa: E731
    for n in nodes:
        g.append(f'<node id="s{n["id"]}"><data key="type">species</data><data key="smiles">{q(n["smiles"])}</data>'
                 f'<data key="label">{q(n["label"] or n["smiles"])}</data><data key="kind">{q(n["kind"])}</data>'
                 f'<data key="step">{n["f"]}</data><data key="role">{q(n["role"])}</data></node>')
    for k, e in enumerate(edges):
        g.append(f'<node id="r{k}"><data key="type">reaction</data><data key="label">{q(e["eq"])}</data>'
                 f'<data key="family">{q(e["fam"])}</data><data key="step">{e["step"]}</data></node>')
        for x in e["in"]:
            g.append(f'<edge source="s{x}" target="r{k}"/>')
        for x in e["out"]:
            g.append(f'<edge source="r{k}" target="s{x}"/>')
    g += ["</graph>", "</graphml>"]
    (out / "subnet.graphml").write_text("\n".join(g), encoding="utf-8")
    tpl = Path(__file__).resolve().parent / "subnet_view.html"
    sd = Path(__file__).resolve().parent.parent / "web" / "vendor" / "smiles-drawer.min.js"
    if tpl.exists():
        page = tpl.read_text(encoding="utf-8")
        page = page.replace("/*__SMILES_DRAWER__*/", sd.read_text(encoding="utf-8") if sd.exists() else "")
        data = json.dumps({"meta": meta, "nodes": nodes, "edges": edges}, ensure_ascii=False).replace("</", "<\\/")
        page = page.replace("/*__DATA__*/null", data)
        (out / "网页.html").write_text(page, encoding="utf-8")
    lines = [f"查询：{meta['query']}", f"数据库：{meta['db']}", f"小分子池：{'、'.join(meta['pool']) or '无'}",
             f"允许比最少步数多：{meta['extra']} 步", ""]
    for t in meta["targets"]:
        lines.append(f"目标 {t['name']}（{t['smiles']}）：最少 {t['D']} 步；放宽 0…{meta['extra']} 步时的反应步数："
                     f"{' / '.join(f'{n:,}' for n in t['n_steps_by_extra'])}" if t["D"] is not None else f"目标 {t['name']}：不可达")
    lines += ["", f"子网络：{len(nodes):,} 个物种，{len(edges):,} 个反应步（按使用方向计）", "",
              "定义见 scripts/subnet_extract.py 开头的说明。每个反应步：最早第几步 = 它最早能发生的步数；额外步数 = 它要在比最少步数多几步时才进入子网络（0 = 在某条最短路线上）；时间余量 = 它最晚可以推迟到第几步发生（在 k = 设定值时）减去最早第几步。"]
    (out / "说明.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description="提取指定反应物到产物、步数不超过最少步数 + k 的子网络，并导出")
    ap.add_argument("db")
    ap.add_argument("--from", dest="src", nargs="+", required=True, help="起始反应物（可多个）")
    ap.add_argument("--to", dest="dst", nargs="+", required=True, help="目标产物（可多个）")
    ap.add_argument("--extra", type=int, default=0, help="允许比最少步数多几步（默认 0：只要最短路径）")
    ap.add_argument("--pool", nargs="*", default=["H", "H2", "H2O"], help="始终可用的小分子（默认 H H2 H2O；写 none 表示不用）")
    ap.add_argument("--max-c", type=int, default=None, help="只用碳数不超过这个值的物种")
    ap.add_argument("--exclude", nargs="*", default=[], help="不使用这些反应族（例如 R:abs R:beta）；一个方程的所有族都被排除才去掉")
    ap.add_argument("--max-reactions", type=int, default=200000, help="子网络反应步数超过这个值时只导出表格、不生成网页")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    db = Path(a.db)
    d = load(db)
    N = len(d["smiles"])
    src = resolve(d, a.src, db)
    dst = resolve(d, a.dst, db)
    pool = [] if a.pool in ([], ["none"]) else resolve(d, a.pool, db)
    R, P = d["R"], d["P"]
    ok = np.ones(len(R), dtype=bool)
    if a.max_c is not None:
        Cn = np.r_[d["C"], 0]
        ok &= (Cn[np.where(R < 0, N, R)] <= a.max_c).all(axis=1) & (Cn[np.where(P < 0, N, P)] <= a.max_c).all(axis=1)
    if a.exclude:
        ex = set(a.exclude)
        bad = np.array([all(x in ex for x in s.split(";")) for s in d["fam_names"]])
        ok &= ~bad[d["fam"]]
    rid = np.flatnonzero(ok)
    Rk, Pk = np.where(R[rid] < 0, N, R[rid]), np.where(P[rid] < 0, N, P[rid])
    IN, OUT = np.vstack([Rk, Pk]), np.vstack([Pk, Rk])          # both directions
    E = len(IN)
    log(f"steps considered: {E:,} ({len(rid):,} reactions x 2 directions)")
    f, fe = forward(IN, OUT, N, src + pool)
    blocked = np.zeros(N + 1, dtype=bool)
    blocked[pool + src] = True
    budgets = [int(f[t]) + a.extra for t in dst if f[t] < INF]
    S = prepare(IN, OUT, N, fe, max(budgets) if budgets else 0, blocked)
    log(f"steps that can fire within the largest budget: {len(S['usable']):,}")
    keep = np.zeros(E, dtype=bool)
    excess = np.full(E, INF, dtype=np.int32)      # smallest k at which the step enters the sub-network
    slack = np.full(E, INF, dtype=np.int32)
    lmax = np.full(N + 1, -1, dtype=np.int32)
    targets = []
    for t, name in zip(dst, a.dst):
        D = int(f[t]) if f[t] < INF else None
        targets.append({"name": name, "smiles": d["smiles"][t], "node_id": int(t), "D": D, "n_steps_by_extra": []})
        if D is None:
            log(f"target {name}: not reachable")
            continue
        for k in range(a.extra + 1):
            L, Te = deadlines(S, t, D + k)
            sel = (Te >= 0) & (fe <= Te) & ~(IN == t).any(axis=1)     # steps that consume the target are not part of it
            excess = np.where(sel, np.minimum(excess, k), excess)
            if k == a.extra:
                keep |= sel
                slack = np.where(sel, np.minimum(slack, Te - fe), slack)
                lmax = np.maximum(lmax, L)
            targets[-1]["n_steps_by_extra"].append(int(sel.sum()))
            log(f"target {name}: minimum {D} steps; steps within +{k}: {int(sel.sum()):,}")
    ks = np.flatnonzero(keep)
    roles = {}
    for x in pool:
        roles[x] = "小分子池"
    for x in src:
        roles[x] = "起始反应物"
    for x in dst:
        roles[x] = "目标产物"
    used = set()
    edges = []
    lab = lambda x: d["label"][x] if d["smiles"][x] in ("[H]", "[H][H]", "O", "O=C=O", "[C-]#[O+]") else d["smiles"][x]  # noqa: E731
    for e in ks.tolist():
        r = int(rid[e % len(rid)])
        ins = [int(x) for x in IN[e] if x != N]
        outs = [int(x) for x in OUT[e] if x != N]
        used.update(ins)
        used.update(outs)
        edges.append({"rid": r, "in": ins, "out": outs, "eq": " + ".join(lab(x) for x in ins) + " → " + " + ".join(lab(x) for x in outs),
                      "fam": d["fam_names"][d["fam"][r]], "step": int(fe[e]), "extra": int(excess[e]), "slack": int(slack[e]),
                      "n_original": int(d["nor"][r])})
    edges.sort(key=lambda e: (e["extra"], e["step"], e["eq"]))
    nodes = []
    for x in sorted(used, key=lambda x: (int(f[x]), d["smiles"][x])):
        nodes.append({"id": x, "smiles": d["smiles"][x], "label": d["label"][x], "kind": d["kind"][x], "formula": d["formula"][x],
                      "C": int(d["C"][x]), "H": int(d["H"][x]), "O": int(d["O"][x]), "f": int(f[x]),
                      "L": None if lmax[x] < 0 else int(lmax[x]), "role": roles.get(x, "中间体")})
    meta = {"query": " ".join(sys.argv[1:] if argv is None else argv), "db": str(db), "extra": a.extra,
            "sources": [{"name": n, "smiles": d["smiles"][x], "node_id": x} for n, x in zip(a.src, src)],
            "pool": [d["label"][x] or d["smiles"][x] for x in pool], "targets": targets,
            "n_species": len(nodes), "n_steps": len(edges), "max_c": a.max_c, "exclude": a.exclude}
    out = Path(a.out or f"subnet_{time.strftime('%Y%m%d_%H%M%S')}")
    if len(edges) > a.max_reactions:
        log(f"{len(edges):,} steps > --max-reactions {a.max_reactions:,}: tables only, no web page")
    write_all(out, d, {"nodes": nodes, "edges": edges}, meta)
    if len(edges) > a.max_reactions and (out / "网页.html").exists():
        (out / "网页.html").unlink()
    log(f"written {out}: {len(nodes):,} species, {len(edges):,} steps")
    print((out / "说明.txt").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
