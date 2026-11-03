"""Multi-level (zoom-to-load) tiles for a generated network.

    python lod_build.py ../data/net_C16_full lod_full        # -> lod_full/meta.json, lod_full/l/*.bin, lod_full/x/*.bin

Species get "grid ids" (gid): core/pool first, then cell by cell (row, carbon number) in layout order,
so the page can place species procedurally from (cell, local index) without downloading coordinates.

  meta.json   rows, columns, cells (box, gid range, count), class/family lists, core positions,
              aggregated cell-to-cell reaction flows for the zoomed-out view, chunk table
  l/<k>.bin   gzip: per species of chunk k -> class, distance from CO2 (steps), log10(#shortest pathways), label
  x/<k>.bin   gzip: per species of chunk k -> every reaction it takes part in (family, its side, the other
              participants as zigzag-varint gid deltas with their side)

Distances and pathway counts are computed here with a level-synchronous hypergraph BFS from CO2 with
H2, H*, H+, H2O, [H], H. always available (same rules as the page used before).
"""
import gzip
import json
import math
import struct
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

T0 = time.time()
src = Path(sys.argv[1] if len(sys.argv) > 1 else "../data/net_C16_full")
dst = Path(sys.argv[2] if len(sys.argv) > 2 else "lod_full")
CH = int(sys.argv[3]) if len(sys.argv) > 3 else 40000
stats = json.loads((src / "stats.json").read_text())
CMAX, FULL = stats["cmax"], stats["full_upto"]


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


ROWS = [
    ("alc", ["alcohol"], "醇 / 烷氧基 RO*", "MET"), ("cbo", ["aldehyde", "ketone", "acid"], "醛（酰基 RCO*）/ 酮 / 羧酸", "MET"),
    ("ads", ["ads"], "吸附烷基 R*", "MET"), ("rad", ["rad"], "气相自由基 R•", "GAS"), ("alk", ["alkane"], "烷烃", "OX"),
    ("ole", ["alkene"], "烯烃", "OX"), ("ia", ["ion_alkyl"], "烷基碳正离子", "BA"), ("die", ["diene"], "二烯", "BA2"),
    ("ie", ["ion_alkenyl"], "烯基 / 烯丙基碳正离子", "BA"), ("cyc", ["naphthene", "cycloalkene"], "环烷 / 环烯", "BA2"),
    ("ic", ["ion_cyc"], "环状碳正离子", "BA"), ("cd", ["cyclodiene"], "环二烯", "BA2"), ("aro", ["aromatic"], "芳烃", "AR"),
]
CLASSES = ["alcohol", "aldehyde", "ketone", "acid", "ads", "rad", "alkane", "alkene", "ion_alkyl", "diene",
           "ion_alkenyl", "naphthene", "cycloalkene", "ion_cyc", "cyclodiene", "aromatic", "core", "pool"]
FAMS = ["Z:prot", "Z:hshift", "Z:alkylshift", "Z:pcp", "Z:beta", "Z:meth", "Z:cyclise", "Z:HT", "H:hyd",
        "M:betaH", "M:CH2", "M:CO", "M:hyd", "O:dehyd", "O:hyd", "O:ox", "R:abs", "R:beta",
        "A:meth", "A:iso", "A:FC", "A:hcp", "C1", "O:OH", "A:cyc"]
cls2row = {c: i for i, r in enumerate(ROWS) for c in r[1]}
CI = {c: i for i, c in enumerate(CLASSES)}

# ---------------------------------------------------------------- species -> gid
tid_cls, tid_c, tid_lab, tid_sk = [], [], [], []
with gzip.open(src / "species.tsv.gz", "rt") as fh:
    next(fh)
    for line in fh:
        i, c, n, s, sk = line.rstrip("\n").split("\t")
        tid_cls.append(c); tid_c.append(int(n)); tid_lab.append(s); tid_sk.append(sk)
NS = len(tid_cls)
log("species", NS)
C1_KEYS = {("alcohol", "CO"): "CH₃OH", ("aldehyde", "C=O"): "HCHO", ("acid", "O=CO"): "HCOOH",
           ("ads", "[CH3+]"): "CH₃*", ("rad", "[CH3+]"): "·CH₃", ("alkane", "C"): "CH₄"}
special = [t for t in range(NS) if tid_cls[t] in ("pool", "core") or tid_c[t] <= 1]
cells = defaultdict(list)
for t in range(NS):
    if tid_cls[t] in ("pool", "core") or tid_c[t] <= 1:
        continue
    cells[(cls2row[tid_cls[t]], tid_c[t])].append(t)
for (r, c), lst in cells.items():
    lst.sort(key=lambda t: (ROWS[r][1].index(tid_cls[t]), tid_sk[t].count("("), tid_sk[t], tid_lab[t]))
order = list(special)
cell_list = sorted(cells)
cell_base = {}
for k in cell_list:
    cell_base[k] = len(order)
    order.extend(cells[k])
gid_of = np.empty(NS, dtype=np.int64)
gid_of[np.array(order)] = np.arange(NS)
g_cls = np.array([CI[tid_cls[t]] for t in order], dtype=np.uint8)
g_c = np.array([tid_c[t] for t in order], dtype=np.uint8)


def label_of(t):
    c, s = tid_cls[t], tid_lab[t]
    if (c, s) in C1_KEYS:
        return C1_KEYS[(c, s)]
    if c == "rad":
        return s.replace("+]", "•]")
    if c == "ads":
        return s.replace("+]", "*]")
    return s


g_lab = [label_of(t) for t in order]
POOL = np.where(g_cls == CI["pool"])[0]

# ---------------------------------------------------------------- layout (same formula as web_export.py)
W = 100.0
cnt = {k: len(v) for k, v in cells.items()}
colmax = {c: max([cnt.get((r, c), 0) for r in range(len(ROWS))] + [1]) for c in range(2, CMAX + 1)}
rowmax = {r: max([cnt.get((r, c), 0) for c in range(2, CMAX + 1)] + [1]) for r in range(len(ROWS))}
C1W = 6.0
wraw = {c: max(1.6, math.sqrt(colmax[c])) for c in range(2, CMAX + 1)}
sc = (W - C1W) / sum(wraw.values())
colx, x = {1: (0.0, C1W)}, C1W
for c in range(2, CMAX + 1):
    colx[c] = (x, x + wraw[c] * sc)
    x += wraw[c] * sc
hraw = {r: max(18.0, math.sqrt(rowmax[r])) for r in range(len(ROWS))}
H = round(W * 0.78, 3)
sr = H / sum(hraw.values())
rowy, y = {}, 0.0
for r in range(len(ROWS)):
    rowy[r] = (y, y + hraw[r] * sr)
    y += hraw[r] * sr


def ry(key, f):
    r = [k for k, *_ in ROWS].index(key)
    return rowy[r][0] + f * (rowy[r][1] - rowy[r][0])


fx = lambda v: v * 0.78 / 3.3 * C1W  # noqa: E731
CORE_POS = {"CO₂": (fx(2.1), 0.012 * H), "CO₂^δ−*": (fx(1.0), ry("alc", .30)), "HCOO*": (fx(1.0), ry("alc", .62)),
            "H₂COO*": (fx(1.0), ry("alc", .90)), "H₂CO*": (fx(1.0), ry("cbo", .30)), "CH₃O*": (fx(1.0), ry("cbo", .62)),
            "CH₃OH": (fx(1.0), ry("cbo", .92)), "HCHO": (fx(2.1), ry("cbo", .45)), "HCOOH": (fx(2.1), ry("alc", .75)),
            "CO*": (fx(3.2), ry("alc", .40)), "HCO*": (fx(3.2), ry("cbo", .20)), "HCOH*": (fx(3.2), ry("cbo", .55)),
            "CH*/C*": (fx(3.2), ry("cbo", .88)), "CH₂*": (fx(3.2), ry("ads", .35)), "CH₃*": (fx(2.1), ry("ads", .65)),
            "·CH₃": (fx(2.1), ry("rad", .5)), "CH₄": (fx(2.1), ry("alk", .5)), "CH₃OH/SMS": (fx(1.5), ry("ia", .5)),
            "DME": (fx(3.0), ry("ia", .5)), "烃池 HCP": (fx(2.1), ry("ic", .5))}

# ---------------------------------------------------------------- reactions
log("reading reactions")
NONE = -1
rf_l, R4 = [], []
with gzip.open(src / "reactions.tsv.gz", "rt") as fh:
    next(fh)
    for line in fh:
        _, fam, r, p, _ = line.split("\t", 4)
        rs, ps = r.split(","), p.split(",")
        rf_l.append(FAMS.index(fam))
        R4.append((int(rs[0]), int(rs[1]) if len(rs) > 1 else -1, int(ps[0]), int(ps[1]) if len(ps) > 1 else -1))
rf = np.array(rf_l, dtype=np.uint8)
del rf_l
A = np.array(R4, dtype=np.int64)
del R4
NR = len(rf)
for j in range(4):
    m = A[:, j] >= 0
    A[m, j] = gid_of[A[m, j]]
log("reactions", NR)

# ---------------------------------------------------------------- levels + shortest-pathway counts
is_pool = np.zeros(NS, dtype=bool)
is_pool[POOL] = True
lev = np.full(NS, -1, dtype=np.int32)
co2 = g_lab.index("CO₂")
lev[POOL] = 0
lev[co2] = 0
one_way = rf == FAMS.index("A:hcp")
# directed edges: forward (in = reactants) and backward (in = products), backward of one-way families dropped
din = np.concatenate([A[:, 0:2], A[~one_way][:, 2:4]])
dout = np.concatenate([A[:, 2:4], A[~one_way][:, 0:2]])
dfam = np.concatenate([rf, rf[~one_way]])
elev = np.full(len(din), -1, dtype=np.int32)
L = 0
while True:
    l0 = np.where(din[:, 0] >= 0, lev[np.maximum(din[:, 0], 0)], 0)
    l1 = np.where(din[:, 1] >= 0, lev[np.maximum(din[:, 1], 0)], 0)
    ready = (elev < 0) & (l0 >= 0) & (l1 >= 0) & (np.maximum(l0, l1) == L)
    if not ready.any() and L > 0 and not (lev == L).any():
        break
    elev[ready] = L + 1
    outs = dout[ready].ravel()
    outs = outs[outs >= 0]
    newly = outs[lev[outs] < 0]
    lev[newly] = L + 1
    L += 1
    if L > 200:
        break
log("max level", lev.max(), "unreachable", int((lev < 0).sum()))
# pathway counts (float64, log10 stored): main inputs = non-pool inputs with level == elev-1
cnt_p = np.zeros(NS, dtype=np.float64)
cnt_p[co2] = 1.0
for Lv in range(1, int(lev.max()) + 1):
    sel = np.where(elev == Lv)[0]
    if not len(sel):
        continue
    for j in (0, 1):
        ins = din[sel, j]
        okin = (ins >= 0)
        ins_c = np.maximum(ins, 0)
        okin &= ~is_pool[ins_c] & (lev[ins_c] == Lv - 1)
        if j == 1:
            okin &= din[sel, 1] != din[sel, 0]
        for k in (0, 1):
            outs = dout[sel, k]
            m = okin & (outs >= 0)
            o = outs[m]
            m2 = (lev[o] == Lv) & ~is_pool[o]
            np.add.at(cnt_p, o[m2], cnt_p[ins_c[m][m2]])
lcount = np.where(cnt_p > 0, np.log10(np.maximum(cnt_p, 1e-300)), -1).astype(np.float32)
log("pathway counts done; max log10", float(lcount.max()))

# ---------------------------------------------------------------- aggregated cell flows (zoomed-out view)
cell_of = np.full(NS, -1, dtype=np.int32)
for ci, k in enumerate(cell_list):
    cell_of[cell_base[k]:cell_base[k] + cnt[k]] = ci
GROUP = {f: f.split(":")[0] for f in FAMS}
flows = Counter()
for ja in (0, 1):
    for jb in (2, 3):
        a, b = A[:, ja], A[:, jb]
        m = (a >= 0) & (b >= 0)
        ca, cb = cell_of[np.maximum(a, 0)], cell_of[np.maximum(b, 0)]
        m &= (ca >= 0) & (cb >= 0) & (ca != cb)
        key = (np.minimum(ca, cb)[m].astype(np.int64) * 4096 + np.maximum(ca, cb)[m]) * 32 + rf[m]
        u, n = np.unique(key, return_counts=True)
        for kk, nn in zip(u.tolist(), n.tolist()):
            flows[kk] += nn
flow_list = sorted(([k // 32 // 4096, k // 32 % 4096, k % 32, v] for k, v in flows.items()), key=lambda x: -x[3])
log("cell flows", len(flow_list))

# ---------------------------------------------------------------- incidence per species, chunked
log("building incidence")
inc_sp, inc_r, inc_side = [], [], []
for j, side in ((0, 0), (1, 0), (2, 1), (3, 1)):
    col = A[:, j]
    m = col >= 0
    if j in (1, 3):
        m &= col != A[:, j - 1]
    sp = col[m]
    keep = ~is_pool[sp]
    inc_sp.append(sp[keep]); inc_r.append(np.where(m)[0][keep]); inc_side.append(np.full(keep.sum(), side, dtype=np.uint8))
inc_sp = np.concatenate(inc_sp); inc_r = np.concatenate(inc_r); inc_side = np.concatenate(inc_side)
o = np.lexsort((inc_r, inc_sp))
inc_sp, inc_r, inc_side = inc_sp[o], inc_r[o], inc_side[o]
starts = np.searchsorted(inc_sp, np.arange(NS + 1))
log("incidences", len(inc_sp))


def uvar(v, out):
    while v >= 0x80:
        out.append((v & 0x7F) | 0x80)
        v >>= 7
    out.append(v)


(dst / "l").mkdir(parents=True, exist_ok=True)
(dst / "x").mkdir(parents=True, exist_ok=True)
nchunks = (NS + CH - 1) // CH
sizes = Counter()
Al = A.tolist()
rfl = rf.tolist()
for k in range(nchunks):
    g0, g1 = k * CH, min(NS, (k + 1) * CH)
    n = g1 - g0
    labs = "\n".join(g_lab[g0:g1]).encode()
    blob = struct.pack("<II", n, len(labs)) + g_cls[g0:g1].tobytes() + lev[g0:g1].astype(np.int16).tobytes()
    blob += b"\0" * ((-len(blob)) % 4) + lcount[g0:g1].tobytes() + labs
    lb = gzip.compress(blob, 9)
    (dst / "l" / f"{k}.bin").write_bytes(lb)
    out = bytearray(struct.pack("<I", n))
    for g in range(g0, g1):
        s0, s1 = starts[g], starts[g + 1]
        uvar(int(s1 - s0), out)
        for q in range(s0, s1):
            r = int(inc_r[q]); side = int(inc_side[q])
            a0, a1, p0, p1 = Al[r]
            others = [(x_, 0) for x_ in (a0, a1) if x_ >= 0] + [(x_, 1) for x_ in (p0, p1) if x_ >= 0]
            others.remove((g, side))
            out.append(rfl[r] | (side << 5))
            out.append(len(others))
            for x_, sd in others:
                d = x_ - g
                z = (d << 1) if d >= 0 else ((-d << 1) - 1)
                uvar((z << 1) | sd, out)
    xb = gzip.compress(bytes(out), 9)
    (dst / "x" / f"{k}.bin").write_bytes(xb)
    sizes["l"] += len(lb); sizes["x"] += len(xb)
    if k % 10 == 0:
        log(f"chunk {k}/{nchunks}  l {sizes['l'] / 1e6:.1f} MB  x {sizes['x'] / 1e6:.1f} MB")

def geom(key):
    r, c = key
    x0, x1 = colx[c]
    y0, y1 = rowy[r]
    pad = 0.012 * (x1 - x0) + 0.05
    top = 0.035 * (y1 - y0) + 0.18
    w, h = x1 - x0 - 2 * pad, y1 - y0 - top - pad
    n = cnt[key]
    ncol = max(1, round(math.sqrt(n * w / h)))
    nrow = math.ceil(n / ncol)
    ncol = math.ceil(n / nrow)
    return {"px": x0 + pad, "py": y0 + top, "dx": w / ncol, "dy": h / nrow, "ncol": ncol, "nrow": nrow}


def alias_gids():
    from rdkit import Chem
    sys.path.insert(0, str(Path(__file__).parent))
    names = [("乙烯", "C=C"), ("丙烯", "C=CC"), ("1-丁烯", "C=CCC"), ("2-丁烯", "CC=CC"), ("异丁烯", "C=C(C)C"),
             ("1,3-丁二烯", "C=CC=C"), ("乙烷", "CC"), ("丙烷", "CCC"), ("正丁烷", "CCCC"), ("异丁烷", "CC(C)C"),
             ("异辛烷", "CC(C)CC(C)(C)C"), ("正辛烷", "CCCCCCCC"), ("正癸烷", "CCCCCCCCCC"), ("正十二烷", "CCCCCCCCCCCC"),
             ("正十四烷", "CCCCCCCCCCCCCC"), ("正十六烷", "CCCCCCCCCCCCCCCC"), ("2-甲基十五烷", "CCCCCCCCCCCCCC(C)C"),
             ("乙醇", "CCO"), ("正丙醇", "CCCO"), ("异丙醇", "CC(C)O"), ("乙醛", "CC=O"), ("丙酮", "CC(C)=O"), ("乙酸", "CC(=O)O"),
             ("苯", "c1ccccc1"), ("甲苯", "Cc1ccccc1"), ("邻二甲苯", "Cc1ccccc1C"), ("间二甲苯", "Cc1cccc(C)c1"),
             ("对二甲苯", "Cc1ccc(C)cc1"), ("乙苯", "CCc1ccccc1"), ("均三甲苯", "Cc1cc(C)cc(C)c1"), ("偏三甲苯", "Cc1ccc(C)c(C)c1"),
             ("均四甲苯", "Cc1cc(C)c(C)cc1C"), ("五甲苯", "Cc1cc(C)c(C)c(C)c1C"), ("六甲苯", "Cc1c(C)c(C)c(C)c(C)c1C"),
             ("异丙苯", "CC(C)c1ccccc1"), ("萘", "c1ccc2ccccc2c1"), ("2-甲基萘", "Cc1ccc2ccccc2c1"), ("2,6-二甲基萘", "Cc1ccc2cc(C)ccc2c1"),
             ("茚满", "c1ccc2c(c1)CCC2"), ("四氢萘", "c1ccc2c(c1)CCCC2"), ("正癸基苯", "CCCCCCCCCCc1ccccc1"),
             ("环己烷", "C1CCCCC1"), ("甲基环己烷", "CC1CCCCC1"), ("环己烯", "C1=CCCCC1"), ("1,3-环己二烯", "C1=CCCC=C1"),
             ("异丙基正离子", "C[CH+]C"), ("叔丁基正离子", "C[C+](C)C"), ("甲基环戊基正离子", "C[C+]1CCCC1")]
    want = {Chem.MolToSmiles(Chem.MolFromSmiles(sm)): nm for nm, sm in names}
    out = {}
    for g, lab in enumerate(g_lab):
        if lab in want and CLASSES[g_cls[g]] not in ("rad", "ads"):
            out.setdefault(want[lab], g)
    return out


meta = {
    "cmax": CMAX, "full_upto": FULL, "W": W, "H": H, "chunk": CH, "nchunks": nchunks, "n_species": NS, "n_reactions": NR,
    "n_pool": int(len(POOL)), "classes": CLASSES, "fams": FAMS,
    "rows": [{"key": k, "label": lab, "block": blk, "cls": cl, "y0": rowy[r][0], "y1": rowy[r][1],
              "n": sum(cnt.get((r, c), 0) for c in range(2, CMAX + 1))} for r, (k, cl, lab, blk) in enumerate(ROWS)],
    "cols": [{"c": c, "x0": colx[c][0], "x1": colx[c][1], "n": sum(cnt.get((r, c), 0) for r in range(len(ROWS))) if c > 1 else len(special)}
             for c in range(1, CMAX + 1)],
    "cells": [dict(r=r, c=c, g0=cell_base[(r, c)], n=cnt[(r, c)], **geom((r, c))) for (r, c) in cell_list],
    "alias": alias_gids(),
    "special": [{"g": g, "label": g_lab[g], "cls": CLASSES[g_cls[g]], "xy": CORE_POS.get(g_lab[g])} for g in range(len(special))],
    "fam_counts": dict(Counter(FAMS[f] for f in rf.tolist())),
    "flows": flow_list[:6000],
    "max_level": int(lev.max()),
    "sizes_MB": {k: round(v / 1e6, 1) for k, v in sizes.items()},
}
(dst / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, separators=(",", ":")))
log("done", meta["sizes_MB"], "chunks", nchunks)
