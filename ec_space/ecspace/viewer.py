"""Browser preview: zoom-to-load tiles (same format as the CO2 network page) built from the database.

Includes every species of the layers <= viewer.max_heavy (or a per-cell sample when there are more than
viewer.max_species), every reaction among them, shortest-path distances / counts / parents from the start
species (EC) with e-, Li+, H., H+, H- always available. Output: <out>/viewer/index.html (+ meta.js, l/*.js, x/*.js).
Open index.html directly from disk (no server needed) or serve the folder.
"""
from __future__ import annotations

import base64
import gzip
import json
import math
import multiprocessing as mp
import shutil
import struct
import zlib
from collections import Counter
from pathlib import Path

import numpy as np

from .pipeline import log

ROWS = [
    ("closed", "闭壳层中性分子", ["closed_shell"], "OX"),
    ("rad", "自由基 R•", ["radical"], "GAS"),
    ("an", "负离子 R⁻ / 二价阴离子", ["anion", "dianion"], "BA"),
    ("cat", "正离子 R⁺ / 质子化 MH⁺", ["cation", "protonated", "zwitterion"], "MET"),
    ("radion", "自由基阴 / 阳离子 M•⁻ / M•⁺", ["radical_anion", "radical_cation"], "BA2"),
    ("dist", "远端自由基离子 / 双自由基", ["distonic_radical_anion", "distonic_radical_cation", "biradical"], "BA2"),
    ("li1", "Li⁺ 配位物（1 个 Li）", [], "AR"),
    ("li2", "Li⁺ 配位物（2 个 Li）", [], "AR"),
]
CLASS_COL = {"closed": "#2F3338", "rad": "#6E7380", "an": "#1F5FBF", "cat": "#D64545", "radion": "#8C5AD8",
             "dist": "#13A39A", "li1": "#E07B00", "li2": "#B85C00", "pool": "#8C95A3"}
CLASS_CN = {k: lab for k, lab, *_ in ROWS} | {"pool": "电子 / 离子池"}
FAMS = ["h_loss", "deprotonation", "hydride_loss", "protonation", "electron_capture", "electron_loss", "site_electron_transfer",
        "li_binding", "homolysis", "heterolysis", "beta_scission", "radical_anion_fragmentation", "radical_cation_fragmentation", "pair_steps"]
FAM_CN = {"h_loss": "失 / 得 H•", "deprotonation": "去质子 / 质子化 C、O", "hydride_loss": "失 / 得 H⁻", "protonation": "O 质子化",
          "electron_capture": "得电子 M•⁻", "electron_loss": "失电子 M•⁺", "site_electron_transfer": "位点电子转移",
          "li_binding": "Li⁺ 配位", "homolysis": "均裂 / 自由基复合", "heterolysis": "异裂 / 离子复合", "beta_scission": "自由基 β-断裂 / 加成",
          "radical_anion_fragmentation": "自由基阴离子断键 / 开环", "radical_cation_fragmentation": "自由基阳离子断键", "pair_steps": "双位点物种逐步形成"}
FAM_COL = {"h_loss": "#8A9099", "deprotonation": "#3A78D4", "hydride_loss": "#C2385E", "protonation": "#D99A00", "electron_capture": "#6A63DB",
           "electron_loss": "#B04848", "site_electron_transfer": "#8C5AD8", "li_binding": "#E07B00", "homolysis": "#5F6670",
           "heterolysis": "#2C5DB8", "beta_scission": "#13A39A", "radical_anion_fragmentation": "#1E9C4B",
           "radical_cation_fragmentation": "#C77A2E", "pair_steps": "#7FCB97"}
FAM_GROUPS = [["电子转移", ["electron_capture", "electron_loss", "site_electron_transfer"]],
              ["H 与质子转移", ["h_loss", "deprotonation", "hydride_loss", "protonation"]],
              ["Li⁺", ["li_binding"]],
              ["断键 / 成键", ["homolysis", "heterolysis", "beta_scission", "radical_anion_fragmentation", "radical_cation_fragmentation"]],
              ["其他", ["pair_steps"]]]
POOL_KEYS = ["[e-]", "[Li+]", "[H]", "[H+]", "[H-]"]
POOL_LABEL = {"[e-]": "e⁻", "[Li+]": "Li⁺", "[H]": "H•", "[H+]": "H⁺", "[H-]": "H⁻"}
ALIAS = {"O=C1OCCO1": "EC 碳酸乙烯酯", "O=C=O": "CO₂", "[C-]#[O+]": "CO", "C=C": "乙烯", "O": "H₂O", "OCCO": "乙二醇",
         "O=CO": "甲酸", "O=C(O)O": "碳酸", "C": "甲烷", "CO": "甲醇", "CCO": "乙醇", "C=O": "甲醛", "CC=O": "乙醛",
         "O=C1OC=CO1": "VC 碳酸亚乙烯酯", "COC(=O)OC": "DMC 碳酸二甲酯", "CCOC(=O)OCC": "DEC 碳酸二乙酯",
         "O=C(O)OCCOC(=O)O": "乙二醇双碳酸（LEDC 母体）"}


def row_of(t):
    if t.startswith("li1_"):
        return 6
    if t.startswith("li2_"):
        return 7
    for i, (_, _, ts, _) in enumerate(ROWS):
        if t in ts:
            return i
    return 0


def _uvar(v, out):
    while v >= 0x80:
        out.append((v & 0x7F) | 0x80)
        v >>= 7
    out.append(v)


def build(cfg):
    out = cfg.out
    vmax = min(cfg.max_heavy, cfg.derive_max_heavy, getattr(cfg, "viewer_max_heavy", 99))
    vdir = out / "viewer"
    if vdir.exists():
        shutil.rmtree(vdir)
    (vdir / "l").mkdir(parents=True)
    (vdir / "x").mkdir(parents=True)
    # ---- species
    keys, types, heavy = [], [], []
    for n in range(1, vmax + 1):
        for f in sorted((out / "species" / f"L{n}").glob("part-*.tsv.gz")):
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 9:
                        continue
                    keys.append(p[0]); types.append(p[1]); heavy.append(int(p[8]))
    log(f"viewer: {len(keys):,} species in layers <= {vmax}")
    # per-cell sampling if needed
    cell_members = {}
    for i, (t, h) in enumerate(zip(types, heavy)):
        cell_members.setdefault((row_of(t), max(1, h)), []).append(i)
    total = len(keys)
    cap = cfg.max_per_cell
    if total > cfg.max_species:
        cap = min(cap, max(50, int(cfg.max_species / max(1, len(cell_members)))))
    start_key = cfg.start
    keep = []
    for k, mem in cell_members.items():
        if len(mem) > cap:
            mem = sorted(mem, key=lambda i: zlib.crc32(keys[i].encode()))
            st = [i for i in mem if keys[i] == start_key]
            mem = st + [i for i in mem[:cap] if keys[i] != start_key]
        cell_members[k] = sorted(mem, key=lambda i: (keys[i].count("("), len(keys[i]), keys[i]))
        keep.extend(cell_members[k])
    sampled = sum(len(v) for v in cell_members.values()) < total
    # ---- ids: pool first, then cells in (row, heavy) order
    cell_list = sorted(cell_members)
    order = []
    gid_of_key = {}
    for k in POOL_KEYS:
        gid_of_key[k] = len(order)
        order.append(("pool", k, 0))
    cell_base = {}
    for ck in cell_list:
        cell_base[ck] = len(order)
        for i in cell_members[ck]:
            gid_of_key[keys[i]] = len(order)
            order.append((types[i], keys[i], heavy[i]))
    NS = len(order)
    del keys, types, heavy
    classes = [r[0] for r in ROWS] + ["pool"]
    g_cls = np.array([classes.index("pool") if t == "pool" else row_of(t) for t, _, _ in order], dtype=np.uint8)
    # ---- reactions among included species
    rf, A = [], []
    fi = {f: i for i, f in enumerate(FAMS)}
    for n in range(1, vmax + 1):
        for f in sorted((out / "reactions" / f"L{n}").glob("part-*.tsv.gz")):
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 3 or p[0] not in fi:
                        continue
                    rs, ps = p[1].split(" + "), p[2].split(" + ")
                    if len(rs) > 2 or len(ps) > 2:
                        continue
                    ids = [gid_of_key.get(x) for x in rs + ps]
                    if any(v is None for v in ids):
                        continue
                    rf.append(fi[p[0]])
                    A.append((ids[0], ids[1] if len(rs) > 1 else -1, ids[len(rs)], ids[len(rs) + 1] if len(ps) > 1 else -1))
    rf = np.array(rf, dtype=np.uint8)
    A = np.array(A, dtype=np.int64).reshape(-1, 4)
    NR = len(rf)
    log(f"viewer: {NS:,} species, {NR:,} reactions{' (sampled per cell)' if sampled else ''}")
    # ---- levels from the start species
    is_pool = g_cls == classes.index("pool")
    lev = np.full(NS, -1, dtype=np.int32)
    lev[is_pool] = 0
    start = gid_of_key.get(start_key, -1)
    if start >= 0:
        lev[start] = 0
    din = np.concatenate([A[:, 0:2], A[:, 2:4]])
    dout = np.concatenate([A[:, 2:4], A[:, 0:2]])
    dfam = np.concatenate([rf, rf])
    elev = np.full(len(din), -1, dtype=np.int32)
    L = 0
    while L < 200:
        l0 = np.where(din[:, 0] >= 0, lev[np.maximum(din[:, 0], 0)], 0)
        l1 = np.where(din[:, 1] >= 0, lev[np.maximum(din[:, 1], 0)], 0)
        ready = (elev < 0) & (l0 >= 0) & (l1 >= 0) & (np.maximum(l0, l1) == L)
        if not ready.any() and not (lev == L).any():
            break
        elev[ready] = L + 1
        o = dout[ready].ravel()
        o = o[o >= 0]
        lev[o[lev[o] < 0]] = L + 1
        L += 1
    cnt_p = np.zeros(NS)
    if start >= 0:
        cnt_p[start] = 1.0
    po, pf, pi_ = [], [], []
    for Lv in range(1, int(lev.max()) + 1 if NS else 1):
        sel = np.where(elev == Lv)[0]
        for j in (0, 1):
            ins = din[sel, j]
            ic = np.maximum(ins, 0)
            ok = (ins >= 0) & ~is_pool[ic] & (lev[ic] == Lv - 1)
            if j == 1:
                ok &= din[sel, 1] != din[sel, 0]
            for k in (0, 1):
                outs = dout[sel, k]
                oc = np.maximum(outs, 0)
                m = ok & (outs >= 0) & (lev[oc] == Lv) & ~is_pool[oc]
                np.add.at(cnt_p, outs[m], cnt_p[ic[m]])
                po.append(outs[m]); pf.append(dfam[sel][m].astype(np.int64)); pi_.append(ins[m])
    lcount = np.where(cnt_p > 0, np.log10(np.maximum(cnt_p, 1e-300)), -1).astype(np.float32)
    if po:
        po, pf, pi_ = np.concatenate(po), np.concatenate(pf), np.concatenate(pi_)
        key_ = (po * 64 + pf) * NS + pi_
        _, first = np.unique(key_, return_index=True)
        po, pf, pi_ = po[first], pf[first], pi_[first]
        o = np.lexsort((pi_, po))
        po, pf, pi_ = po[o], pf[o], pi_[o]
        gs = np.searchsorted(po, po, side="left")
        kp = (np.arange(len(po)) - gs) < 12
        po, pf, pi_ = po[kp], pf[kp], pi_[kp]
    else:
        po = pf = pi_ = np.zeros(0, dtype=np.int64)
    par_start = np.searchsorted(po, np.arange(NS + 1))
    log(f"viewer: max distance from start {int(lev.max())}, unreachable {int((lev < 0).sum()):,}")
    # ---- layout: column 0 = pool panel, columns 1..vmax = heavy atoms
    cnt = {ck: len(cell_members[ck]) for ck in cell_list}
    W = 100.0
    colmax = {c: max([cnt.get((r, c), 0) for r in range(len(ROWS))] + [1]) for c in range(1, vmax + 1)}
    rowmax = {r: max([cnt.get((r, c), 0) for c in range(1, vmax + 1)] + [1]) for r in range(len(ROWS))}
    EXP = 0.5 if NS < 1_000_000 else 1 / 3
    C0W = 5.0
    wraw = {c: max(2.0, colmax[c] ** EXP) for c in range(1, vmax + 1)}
    sc = (W - C0W) / sum(wraw.values())
    colx, x = {0: (0.0, C0W)}, C0W
    for c in range(1, vmax + 1):
        colx[c] = (x, x + wraw[c] * sc)
        x += wraw[c] * sc
    hraw = {r: max(8.0, rowmax[r] ** EXP) for r in range(len(ROWS))}
    H = round(W * 0.6, 3)
    sr = H / sum(hraw.values())
    rowy, y = {}, 0.0
    for r in range(len(ROWS)):
        rowy[r] = (y, y + hraw[r] * sr)
        y += hraw[r] * sr

    def geom(ck):
        r, c = ck
        x0, x1 = colx[c]
        y0, y1 = rowy[r]
        pad = min(0.012 * (x1 - x0) + 0.05, 0.08 * (x1 - x0), 0.08 * (y1 - y0))
        top = min(0.035 * (y1 - y0) + 0.18, 0.18 * (y1 - y0))
        w, h = x1 - x0 - 2 * pad, y1 - y0 - top - pad
        n = cnt[ck]
        ncol = max(1, round(math.sqrt(n * w / h)))
        nrow = math.ceil(n / ncol)
        ncol = math.ceil(n / nrow)
        return {"px": x0 + pad, "py": y0 + top, "dx": w / ncol, "dy": h / nrow, "ncol": ncol, "nrow": nrow}

    special = [{"g": i, "label": POOL_LABEL[k], "cls": "pool", "xy": [C0W * 0.5, H * (0.12 + 0.17 * i)]} for i, k in enumerate(POOL_KEYS)]
    labels = [POOL_LABEL.get(k, ALIAS.get(k.split("|")[0], k) + ({"|-1": " •⁻", "|+1": " •⁺"}.get("|" + k.split("|")[1], "") if "|" in k else ""))
              if t == "pool" else (ALIAS.get(k) or (ALIAS.get(k.split("|")[0]) + {"-1": " •⁻", "+1": " •⁺"}[k.split("|")[1]]
                                                    if "|" in k and k.split("|")[0] in ALIAS and k.split("|")[1] in ("-1", "+1") else k))
              for t, k, _ in order]
    # ---- chunks (.js twins only: work from disk and from static hosting)
    CH = 40000
    nchunks = (NS + CH - 1) // CH
    inc_sp, inc_r, inc_side = [], [], []
    for j, side in ((0, 0), (1, 0), (2, 1), (3, 1)):
        col = A[:, j] if NR else np.zeros(0, dtype=np.int64)
        m = col >= 0
        if j in (1, 3):
            m &= col != A[:, j - 1]
        sp = col[m]
        kp = ~is_pool[sp]
        inc_sp.append(sp[kp]); inc_r.append(np.where(m)[0][kp]); inc_side.append(np.full(int(kp.sum()), side, dtype=np.uint8))
    inc_sp = np.concatenate(inc_sp) if NR else np.zeros(0, dtype=np.int64)
    inc_r = np.concatenate(inc_r) if NR else np.zeros(0, dtype=np.int64)
    inc_side = np.concatenate(inc_side) if NR else np.zeros(0, dtype=np.uint8)
    o = np.lexsort((inc_r, inc_sp))
    inc_sp, inc_r, inc_side = inc_sp[o], inc_r[o], inc_side[o]
    starts = np.searchsorted(inc_sp, np.arange(NS + 1))
    sizes = Counter()
    for k in range(nchunks):
        g0, g1 = k * CH, min(NS, (k + 1) * CH)
        n = g1 - g0
        labs = "\n".join(labels[g0:g1]).encode()
        blob = struct.pack("<II", n, len(labs)) + g_cls[g0:g1].tobytes() + lev[g0:g1].astype(np.int16).tobytes()
        blob += b"\0" * ((-len(blob)) % 4) + lcount[g0:g1].tobytes() + labs
        pb = bytearray()
        for g in range(g0, g1):
            a0, a1 = int(par_start[g]), int(par_start[g + 1])
            _uvar(a1 - a0, pb)
            for q in range(a0, a1):
                d = int(pi_[q]) - g
                pb.append(int(pf[q]))
                _uvar((d << 1) if d >= 0 else ((-d << 1) - 1), pb)
        blob += bytes(pb)
        lb = gzip.compress(blob, 6)
        s0, s1 = int(starts[g0]), int(starts[g1])
        rr = inc_r[s0:s1]
        rows_, fams_, sds = A[rr].tolist(), rf[rr].tolist(), inc_side[s0:s1].tolist()
        cnts = (starts[g0 + 1:g1 + 1] - starts[g0:g1]).tolist()
        xb_ = bytearray(struct.pack("<I", n))
        q = 0
        for i in range(n):
            g = g0 + i
            _uvar(cnts[i], xb_)
            for _ in range(cnts[i]):
                a0, a1, p0, p1 = rows_[q]
                side = sds[q]
                others = [(x_, 0) for x_ in (a0, a1) if x_ >= 0] + [(x_, 1) for x_ in (p0, p1) if x_ >= 0]
                others.remove((g, side))
                xb_.append(fams_[q] | (side << 5))
                xb_.append(len(others))
                for x_, sd in others:
                    d = x_ - g
                    z = (d << 1) if d >= 0 else ((-d << 1) - 1)
                    _uvar((z << 1) | sd, xb_)
                q += 1
        xb = gzip.compress(bytes(xb_), 6)
        for kind, data in (("l", lb), ("x", xb)):
            (vdir / kind / f"{k}.js").write_text(f'window.__lodChunk("{kind}",{k},"{base64.b64encode(data).decode()}");')
            sizes[kind] += len(data)
    # flows between cells
    cell_of = np.full(NS, -1, dtype=np.int32)
    for ci, ck in enumerate(cell_list):
        cell_of[cell_base[ck]:cell_base[ck] + cnt[ck]] = ci
    flows = Counter()
    for ja in (0, 1):
        for jb in (2, 3):
            if not NR:
                continue
            a, b = A[:, ja], A[:, jb]
            m = (a >= 0) & (b >= 0)
            ca, cb = cell_of[np.maximum(a, 0)], cell_of[np.maximum(b, 0)]
            m &= (ca >= 0) & (cb >= 0) & (ca != cb)
            kk = (np.minimum(ca, cb)[m].astype(np.int64) * 4096 + np.maximum(ca, cb)[m]) * 32 + rf[m]
            u, nn = np.unique(kk, return_counts=True)
            for k2, v in zip(u.tolist(), nn.tolist()):
                flows[k2] += v
    flow_list = sorted(([k // 32 // 4096, k // 32 % 4096, k % 32, v] for k, v in flows.items()), key=lambda x: -x[3])[:6000]
    alias = {}
    for k, nm in ALIAS.items():
        if k in gid_of_key:
            alias[nm] = gid_of_key[k]
    meta = {
        "cmax": vmax, "full_upto": vmax, "W": W, "H": H, "chunk": CH, "nchunks": nchunks, "n_species": NS, "n_reactions": NR,
        "n_pool": len(POOL_KEYS), "classes": classes, "fams": FAMS,
        "rows": [{"key": k, "label": lab, "block": blk, "cls": [k], "y0": rowy[r][0], "y1": rowy[r][1],
                  "n": sum(cnt.get((r, c), 0) for c in range(1, vmax + 1))} for r, (k, lab, _, blk) in enumerate(ROWS)],
        "cols": [{"c": c, "x0": colx[c][0], "x1": colx[c][1], "label": ("池" if c == 0 else f"{c} 重原子"),
                  "n": (len(POOL_KEYS) if c == 0 else sum(cnt.get((r, c), 0) for r in range(len(ROWS))))} for c in range(0, vmax + 1)],
        "cells": [dict(r=r, c=c, g0=cell_base[(r, c)], n=cnt[(r, c)], **geom((r, c))) for (r, c) in cell_list],
        "special": special, "alias": alias, "fam_counts": dict(Counter(FAMS[f] for f in rf.tolist())), "flows": flow_list,
        "max_level": int(lev.max()) if NS else 0, "sizes_MB": {k: round(v / 1e6, 1) for k, v in sizes.items()},
        "start_gid": start, "default_gid": start,
        "ui": {"title": f"EC 化学空间 · 重原子 ≤ {vmax}（含自由基、离子、Li⁺ 配位物）",
               "subtitle": ("每个点是一个物种。列为重原子数，行为物种类别。放大后分级加载标签与每条基元反应；"
                            "单击物种查看全部反应，以及从 EC + e⁻ / Li⁺ / H 出发的最短路径。" + ("（每格已抽样）" if sampled else "")),
               "class_cn": CLASS_CN, "class_col": CLASS_COL, "fam_cn": FAM_CN, "fam_col": FAM_COL, "fam_groups": FAM_GROUPS,
               "col_prefix": "重原子 ", "start_label": "EC"},
    }
    (vdir / "meta.js").write_text("window.__lodMeta = " + json.dumps(meta, ensure_ascii=False, separators=(",", ":")) + ";")
    src = Path(__file__).resolve().parent.parent / "viewer"
    shutil.copy(src / "index.html", vdir / "index.html")
    shutil.copy(src / "smiles-drawer.min.js", vdir / "smiles-drawer.min.js")
    shutil.copy(src / "smiles-drawer-LICENSE.md", vdir / "smiles-drawer-LICENSE.md")
    log(f"viewer: written {vdir}  ({nchunks} chunks, labels {sizes['l'] / 1e6:.1f} MB, reactions {sizes['x'] / 1e6:.1f} MB)")
