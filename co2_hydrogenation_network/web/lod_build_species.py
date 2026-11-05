"""Zoom-to-load tiles for a species-level database (species_merge.py output): one node per species, no ions.

    python lod_build_species.py ../data/net_C16_species       out_full  [chunk=40000]
    python lod_build_species.py ../data/net_C16_shape_species out_shape

Writes <out>/meta.js, <out>/l/<k>.js, <out>/x/<k>.js (base64 gzip in <script> files: they load from file:// and
from static hosting), plus index.html, smiles-drawer.min.js and subnet_tpl.js (the sub-network viewer used by the
in-page extraction).  Chunk formats are those of lod_build.py with two extensions:
  l: class byte = row class (low 5 bits) | site mask << 5   (1 acid site R+/R-O-Z, 2 metal site R*, 4 gas radical R.)
  x: a reaction can have up to 3 species per side (others count <= 5)
Levels, shortest-pathway counts and path parents are computed from CO2 with H, H2, H2O always available.
"""
import base64
import gzip
import json
import math
import multiprocessing as mp
import shutil
import struct
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

T0 = time.time()
src = Path(sys.argv[1] if len(sys.argv) > 1 else "../data/net_C16_species")
dst = Path(sys.argv[2] if len(sys.argv) > 2 else "out_species")
CH = int(sys.argv[3]) if len(sys.argv) > 3 else 40000
HERE = Path(__file__).resolve().parent
W3 = 3


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


KIND2CLS = {"R 基团 · 烷基": "r_alkyl", "R 基团 · 烯基 / 烯丙基": "r_alkenyl", "R 基团 · 环状": "r_cyc", "烷烃": "alkane",
            "烯烃": "alkene", "二烯": "diene", "环烷": "naphthene", "环烯": "cycloalkene", "环二烯": "cyclodiene", "芳烃": "aromatic",
            "醇": "alcohol", "醚": "ether", "醛": "aldehyde", "酰基 RC(=O)": "acyl", "酮": "ketone", "羧酸": "acid",
            "C1 碎片": "core", "C1 分子": "core", "氢": "pool", "水": "pool"}
ROWS = [  # key, classes, label, band
    ("alc", ["alcohol", "ether"], "醇 / 醚", "MET"),
    ("cbo", ["aldehyde", "acyl", "ketone", "acid"], "醛 / 酰基 RC(=O) / 酮 / 羧酸", "MET"),
    ("ra", ["r_alkyl"], "R 基团 · 烷基（R⁺ / R* / R•）", "GAS"),
    ("alk", ["alkane"], "烷烃", "OX"),
    ("ole", ["alkene"], "烯烃", "OX"),
    ("re", ["r_alkenyl"], "R 基团 · 烯基 / 烯丙基", "BA"),
    ("die", ["diene"], "二烯", "BA2"),
    ("cyc", ["naphthene", "cycloalkene"], "环烷 / 环烯", "BA2"),
    ("rc", ["r_cyc"], "R 基团 · 环状", "BA"),
    ("cd", ["cyclodiene"], "环二烯", "BA2"),
    ("aro", ["aromatic"], "芳烃", "AR"),
]
CLASSES = ["alcohol", "ether", "aldehyde", "acyl", "ketone", "acid", "r_alkyl", "alkane", "alkene", "r_alkenyl", "diene",
           "naphthene", "cycloalkene", "r_cyc", "cyclodiene", "aromatic", "core", "pool"]
CLASS_CN = {"alcohol": "醇", "ether": "醚", "aldehyde": "醛", "acyl": "酰基 RC(=O)", "ketone": "酮", "acid": "羧酸",
            "r_alkyl": "R 基团 · 烷基", "alkane": "烷烃", "alkene": "烯烃", "r_alkenyl": "R 基团 · 烯基 / 烯丙基", "diene": "二烯",
            "naphthene": "环烷", "cycloalkene": "环烯", "r_cyc": "R 基团 · 环状", "cyclodiene": "环二烯", "aromatic": "芳烃",
            "core": "C1 物种 / 碎片", "pool": "H / H₂ / H₂O（始终可用）"}
CLASS_COL = {"alcohol": "#8E5BB5", "ether": "#B07CD8", "aldehyde": "#C2508A", "acyl": "#D9628F", "ketone": "#C77A2E", "acid": "#6A3D9A",
             "r_alkyl": "#D64545", "alkane": "#2F3338", "alkene": "#E07B00", "r_alkenyl": "#3D8FE0", "diene": "#C9A100",
             "naphthene": "#5C6B73", "cycloalkene": "#2F9C95", "r_cyc": "#6C4FD1", "cyclodiene": "#0E7C66", "aromatic": "#1E8C3A",
             "core": "#14171C", "pool": "#8C95A3"}
BASE_CN = {"Z:prot": "烯烃 + H ⇌ R（酸位质子化 / 去质子）", "Z:hshift": "1,2-H 迁移", "Z:alkylshift": "烷基迁移（A 型）",
           "Z:pcp": "PCP 支化（B 型）", "Z:beta": "β-断裂 / 齐聚（酸位）", "Z:meth": "甲基化（CH₃ + 烯烃）", "Z:cyclise": "烯基成环",
           "Z:HT": "氢转移（放出 H₂ 当量）", "H:hyd": "加氢 / 脱氢", "M:betaH": "β-H 消除（金属位）", "M:CH2": "CH₂ 插入（链增长）",
           "M:CO": "CO 插入 → 酰基", "M:acylH": "酰基加氢 → 醛", "M:hyd": "R + H → 烷烃（金属位）", "O:dehyd": "醇脱水",
           "O:hyd": "羰基加氢", "O:ox": "醛 → 酸", "O:OH": "R + H₂O → 醇 + H", "R:abs": "H 抽提（R + H₂ ⇌ RH + H）",
           "R:beta": "自由基 β-断裂", "A:meth": "芳环甲基化", "A:iso": "芳环 1,2-迁移", "A:FC": "侧链烷基化 / 脱烷基",
           "A:cyc": "侧链脱氢环化", "A:scm": "侧链甲基化（烃池）", "C1": "C1 核心"}
BASE_COL = {"Z:prot": "#5B9BE8", "Z:hshift": "#3A78D4", "Z:alkylshift": "#6A63DB", "Z:pcp": "#8C5AD8", "Z:beta": "#2C5DB8",
            "Z:meth": "#2FB0E0", "Z:cyclise": "#7B4CCB", "Z:HT": "#13A39A", "H:hyd": "#D99A00", "M:betaH": "#D64545", "M:CH2": "#E2603A",
            "M:CO": "#C2385E", "M:acylH": "#D9628F", "M:hyd": "#B04848", "O:dehyd": "#9058C0", "O:hyd": "#C2508A", "O:ox": "#C77A2E",
            "O:OH": "#B07CD8", "R:abs": "#8A9099", "R:beta": "#5F6670", "A:meth": "#1E9C4B", "A:iso": "#45B86D", "A:FC": "#2E8B57",
            "A:cyc": "#5FAE3E", "A:scm": "#7FCB97", "C1": "#14171C"}
GROUPS = [["分子筛 Brønsted 酸", ["Z:prot", "Z:hshift", "Z:alkylshift", "Z:pcp", "Z:beta", "Z:meth", "Z:cyclise"]],
          ["氢转移 / 加氢", ["Z:HT", "H:hyd"]], ["金属 / 碳化物", ["M:betaH", "M:CH2", "M:CO", "M:acylH", "M:hyd"]],
          ["含氧化合物", ["O:dehyd", "O:hyd", "O:ox", "O:OH"]], ["气相自由基", ["R:abs", "R:beta"]],
          ["芳烃（含烃池）", ["A:meth", "A:iso", "A:FC", "A:cyc", "A:scm"]], ["C1", ["C1"]]]
SITE_BITS = (("酸位", 1), ("金属位", 2), ("气相自由基", 4))
cls2row = {c: i for i, r in enumerate(ROWS) for c in r[1]}
CI = {c: i for i, c in enumerate(CLASSES)}

# ---------------------------------------------------------------- nodes -> gid
smi, lab, cls_, cn, sk, site = [], [], [], [], [], []
with gzip.open(src / "nodes.tsv.gz", "rt", encoding="utf-8") as fh:
    next(fh)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        smi.append(p[1]); lab.append(p[2]); cls_.append(KIND2CLS[p[3]]); cn.append(int(p[5])); sk.append(p[10])
        site.append(sum(b for name, b in SITE_BITS if name in p[8]) if p[3].startswith("R 基团") or p[1] == "[CH3]" else 0)
NS = len(smi)
CMAX = max(cn)
log("nodes", NS, "Cmax", CMAX)
special = [t for t in range(NS) if cn[t] <= 1]
special.sort(key=lambda t: (cls_[t] != "pool", t))
cells = defaultdict(list)
for t in range(NS):
    if cn[t] > 1:
        cells[(cls2row[cls_[t]], cn[t])].append(t)
for (r, c), lst in cells.items():
    lst.sort(key=lambda t: (ROWS[r][1].index(cls_[t]), sk[t].count("("), sk[t], smi[t]))
order = list(special)
cell_list = sorted(cells)
cell_base = {}
for k in cell_list:
    cell_base[k] = len(order)
    order.extend(cells[k])
gid_of = np.empty(NS, dtype=np.int64)
gid_of[np.array(order)] = np.arange(NS)
g_cls = np.array([CI[cls_[t]] for t in order], dtype=np.uint8)
g_site = np.array([site[t] for t in order], dtype=np.uint8)
SPECIAL_NAME = {"[H]": "H", "[H][H]": "H₂", "O": "H₂O", "O=C=O": "CO₂", "[C-]#[O+]": "CO", "C": "CH₄", "CO": "CH₃OH",
                "C=O": "HCHO", "O=CO": "HCOOH", "[CH3]": "CH₃（SMS / CH₃* / •CH₃）", "[CH2]": "CH₂", "[CH]": "CH / C",
                "C[O]": "CH₃O", "O[CH]": "HCOH", "[CH]=O": "HCO", "O=C[O]": "HCOO", "[O]C[O]": "H₂COO", "[O]CO": "H₂COOH"}
g_lab = [smi[t] for t in order]
spec_label = [SPECIAL_NAME.get(smi[t], lab[t] or smi[t]) for t in special]
POOL = np.array([i for i, t in enumerate(special) if cls_[t] == "pool"])
is_pool = np.zeros(NS, dtype=bool)
is_pool[POOL] = True

# ---------------------------------------------------------------- layout
Wd = 100.0
cnt = {k: len(v) for k, v in cells.items()}
colmax = {c: max([cnt.get((r, c), 0) for r in range(len(ROWS))] + [1]) for c in range(2, CMAX + 1)}
rowmax = {r: max([cnt.get((r, c), 0) for c in range(2, CMAX + 1)] + [1]) for r in range(len(ROWS))}
C1W = 6.0
EXP = 0.5 if NS < 1_000_000 else 1 / 3
wraw = {c: max(1.6 if EXP == 0.5 else 2.0, colmax[c] ** EXP) for c in range(2, CMAX + 1)}
sc = (Wd - C1W) / sum(wraw.values())
colx, x = {1: (0.0, C1W)}, C1W
for c in range(2, CMAX + 1):
    colx[c] = (x, x + wraw[c] * sc)
    x += wraw[c] * sc
hraw = {r: max(18.0 if EXP == 0.5 else 6.0, rowmax[r] ** EXP) for r in range(len(ROWS))}
Hd = round(Wd * 0.78, 3)
sr = Hd / sum(hraw.values())
rowy, y = {}, 0.0
for r in range(len(ROWS)):
    rowy[r] = (y, y + hraw[r] * sr)
    y += hraw[r] * sr


def ry(key, f):
    r = [k for k, *_ in ROWS].index(key)
    return rowy[r][0] + f * (rowy[r][1] - rowy[r][0])


fx = lambda v: v * 0.78 / 3.3 * C1W  # noqa: E731
CORE_POS = {"H": (fx(1.0), 0.012 * Hd + 0.6), "H₂": (fx(2.1), 0.012 * Hd + 0.6), "H₂O": (fx(3.2), 0.012 * Hd + 0.6),
            "CO₂": (fx(2.1), ry("alc", .25)), "HCOO": (fx(1.0), ry("alc", .55)), "H₂COO": (fx(1.0), ry("alc", .85)),
            "H₂COOH": (fx(1.0), ry("cbo", .2)), "HCOOH": (fx(2.1), ry("alc", .75)), "CO": (fx(3.2), ry("alc", .45)),
            "HCO": (fx(3.2), ry("cbo", .25)), "HCHO": (fx(2.1), ry("cbo", .45)), "CH₃O": (fx(1.0), ry("cbo", .55)),
            "CH₃OH": (fx(1.0), ry("cbo", .85)), "HCOH": (fx(3.2), ry("cbo", .6)), "CH / C": (fx(3.2), ry("cbo", .9)),
            "CH₂": (fx(3.2), ry("ra", .3)), "CH₃（SMS / CH₃* / •CH₃）": (fx(2.1), ry("ra", .6)), "CH₄": (fx(2.1), ry("alk", .5))}
CORE_SHORT = {"CH₃（SMS / CH₃* / •CH₃）": "CH₃"}

# ---------------------------------------------------------------- reactions
files = sorted(src.glob("reactions.part*.tsv.gz")) or [src / "reactions.tsv.gz"]
fam_combo = {}
rows_, blocks = [], []
for f in files:
    with gzip.open(f, "rt", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            p = line.split("\t", 5)
            r = [int(v) for v in p[1].split(",")]
            q = [int(v) for v in p[2].split(",")]
            rows_.append(r + [-1] * (W3 - len(r)) + q + [-1] * (W3 - len(q)) + [fam_combo.setdefault(p[3], len(fam_combo))])
            if len(rows_) >= 1_000_000:
                blocks.append(np.array(rows_, dtype=np.int64))
                rows_ = []
    log("read", f.name)
if rows_:
    blocks.append(np.array(rows_, dtype=np.int64))
Araw = np.concatenate(blocks)
del blocks, rows_
FAMS = [k for k, _ in sorted(fam_combo.items(), key=lambda kv: kv[1])]
assert len(FAMS) <= 32, len(FAMS)
rf = Araw[:, 6].astype(np.uint8)
A = Araw[:, :6].copy()
del Araw
m = A >= 0
A[m] = gid_of[A[m]]
NR = len(A)
log("reactions", NR, "family combinations", len(FAMS))

# ---------------------------------------------------------------- levels + pathway counts + parents (from CO2)
co2 = spec_label.index("CO₂")
lev = np.full(NS, -1, dtype=np.int32)
lev[POOL] = 0
lev[co2] = 0
din = np.concatenate([A[:, :3], A[:, 3:]])
dout = np.concatenate([A[:, 3:], A[:, :3]])
dfam = np.concatenate([rf, rf])
BIGL = 1 << 20
L = 0
while True:
    lv = np.where(din >= 0, lev[np.maximum(din, 0)], 0)
    lv = np.where(lv < 0, BIGL, lv)
    fin = lv.max(axis=1)
    ready = fin == L
    outs = dout[ready].ravel()
    outs = outs[outs >= 0]
    newly = outs[lev[outs] < 0]
    lev[newly] = L + 1
    L += 1
    if not (lev == L).any() or L > 300:
        break
lv = np.where(din >= 0, lev[np.maximum(din, 0)], 0)
lv = np.where(lv < 0, BIGL, lv)
elev = lv.max(axis=1) + 1
elev[elev > BIGL] = -1
log("max level", lev.max(), "unreachable", int((lev < 0).sum()))
cnt_p = np.zeros(NS, dtype=np.float64)
cnt_p[co2] = 1.0
for Lv in range(1, int(lev.max()) + 1):
    sel = np.flatnonzero(elev == Lv)
    if not len(sel):
        continue
    for j in range(W3):
        ins = din[sel, j]
        ic = np.maximum(ins, 0)
        okin = (ins >= 0) & ~is_pool[ic] & (lev[ic] == Lv - 1)
        for jj in range(j):
            okin &= din[sel, jj] != ins
        for k in range(W3):
            outs = dout[sel, k]
            mm = okin & (outs >= 0)
            o = outs[mm]
            m2 = (lev[o] == Lv) & ~is_pool[o]
            np.add.at(cnt_p, o[m2], cnt_p[ic[mm][m2]])
lcount = np.where(cnt_p > 0, np.log10(np.maximum(cnt_p, 1e-300)), -1).astype(np.float32)
RANK = ["C1", "M:CH2", "M:CO", "M:acylH", "M:betaH", "M:hyd", "Z:meth", "Z:prot", "Z:beta", "H:hyd", "Z:HT", "A:meth"]
frank = np.array([min([RANK.index(b) for b in f.split(";") if b in RANK] + [50]) for f in FAMS], dtype=np.int64)
po, pf, pi_ = [], [], []
for j in range(W3):
    ins = din[:, j]
    ic = np.maximum(ins, 0)
    okin = (ins >= 0) & (elev > 0) & ~is_pool[ic] & (lev[ic] == elev - 1)
    for jj in range(j):
        okin &= din[:, jj] != ins
    for k in range(W3):
        outs = dout[:, k]
        oc = np.maximum(outs, 0)
        mm = okin & (outs >= 0) & (lev[oc] == elev) & ~is_pool[oc]
        po.append(outs[mm]); pf.append(dfam[mm].astype(np.int64)); pi_.append(ins[mm])
po, pf, pi_ = np.concatenate(po), np.concatenate(pf), np.concatenate(pi_)
key = (po * 64 + pf) * NS + pi_
key, first = np.unique(key, return_index=True)
po, pf, pi_ = po[first], pf[first], pi_[first]
o = np.lexsort((pi_, frank[pf], po))
po, pf, pi_ = po[o], pf[o], pi_[o]
grp = np.searchsorted(po, po, side="left")
keep = (np.arange(len(po)) - grp) < 12
po, pf, pi_ = po[keep], pf[keep], pi_[keep]
par_start = np.searchsorted(po, np.arange(NS + 1))
del din, dout, dfam, lv
log("pathway counts and parents done")

# ---------------------------------------------------------------- cell flows
cell_of = np.full(NS, -1, dtype=np.int32)
for ci, k in enumerate(cell_list):
    cell_of[cell_base[k]:cell_base[k] + cnt[k]] = ci
flows = Counter()
for ja in range(3):
    for jb in range(3, 6):
        a, b = A[:, ja], A[:, jb]
        mm = (a >= 0) & (b >= 0)
        ca, cb = cell_of[np.maximum(a, 0)], cell_of[np.maximum(b, 0)]
        mm &= (ca >= 0) & (cb >= 0) & (ca != cb)
        kk = (np.minimum(ca, cb)[mm].astype(np.int64) * 4096 + np.maximum(ca, cb)[mm]) * 32 + rf[mm]
        u, n = np.unique(kk, return_counts=True)
        for k2, v in zip(u.tolist(), n.tolist()):
            flows[k2] += v
flow_list = sorted(([k // 32 // 4096, k // 32 % 4096, k % 32, v] for k, v in flows.items()), key=lambda z: -z[3])
log("cell flows", len(flow_list))

# ---------------------------------------------------------------- incidence (each species once per side)
inc_sp, inc_r, inc_side = [], [], []
for j in range(6):
    col = A[:, j]
    mm = col >= 0
    base = 0 if j < 3 else 3
    for jj in range(base, j):
        mm &= col != A[:, jj]
    sp = col[mm]
    kp = ~is_pool[sp]
    inc_sp.append(sp[kp]); inc_r.append(np.flatnonzero(mm)[kp]); inc_side.append(np.full(int(kp.sum()), 0 if j < 3 else 1, dtype=np.uint8))
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


if dst.exists():
    shutil.rmtree(dst)
(dst / "l").mkdir(parents=True)
(dst / "x").mkdir(parents=True)
nchunks = (NS + CH - 1) // CH
labels_all = list(g_lab)
for i, t in enumerate(special):
    labels_all[i] = spec_label[i]


def write_chunk(k):
    g0, g1 = k * CH, min(NS, (k + 1) * CH)
    n = g1 - g0
    labs = "\n".join(labels_all[g0:g1]).encode()
    clsb = (g_cls[g0:g1] | (g_site[g0:g1] << 5)).astype(np.uint8)
    blob = struct.pack("<II", n, len(labs)) + clsb.tobytes() + lev[g0:g1].astype(np.int16).tobytes()
    blob += b"\0" * ((-len(blob)) % 4) + lcount[g0:g1].tobytes() + labs
    pb = bytearray()
    for g in range(g0, g1):
        a0, a1 = int(par_start[g]), int(par_start[g + 1])
        uvar(a1 - a0, pb)
        for q in range(a0, a1):
            d = int(pi_[q]) - g
            pb.append(int(pf[q]))
            uvar((d << 1) if d >= 0 else ((-d << 1) - 1), pb)
    blob += bytes(pb)
    lb = gzip.compress(blob, 9)
    s0, s1 = int(starts[g0]), int(starts[g1])
    rr = inc_r[s0:s1]
    rows, fams, sds = A[rr].tolist(), rf[rr].tolist(), inc_side[s0:s1].tolist()
    counts = (starts[g0 + 1:g1 + 1] - starts[g0:g1]).tolist()
    out = bytearray(struct.pack("<I", n))
    q = 0
    for i in range(n):
        g = g0 + i
        uvar(counts[i], out)
        for _ in range(counts[i]):
            row = rows[q]
            side = sds[q]
            others = [(v, 0) for v in row[:3] if v >= 0] + [(v, 1) for v in row[3:] if v >= 0]
            others.remove((g, side))
            out.append(fams[q] | (side << 5))
            out.append(len(others))
            for v, sd in others:
                d = v - g
                z = (d << 1) if d >= 0 else ((-d << 1) - 1)
                uvar((z << 1) | sd, out)
            q += 1
    xb = gzip.compress(bytes(out), 9)
    for kind, data in (("l", lb), ("x", xb)):
        (dst / kind / f"{k}.js").write_text(f'window.__lodChunk("{kind}",{k},"{base64.b64encode(data).decode()}");')
    return k, len(lb), len(xb)


sizes = Counter()
with mp.get_context("fork").Pool(4) as pool:
    for k, nl, nx in pool.imap_unordered(write_chunk, range(nchunks)):
        sizes["l"] += nl; sizes["x"] += nx
        if k % 10 == 0:
            log(f"chunk {k}/{nchunks}  l {sizes['l'] / 1e6:.1f} MB  x {sizes['x'] / 1e6:.1f} MB")


def geom(key_):
    r, c = key_
    x0, x1 = colx[c]
    y0, y1 = rowy[r]
    pad = min(0.012 * (x1 - x0) + 0.05, 0.08 * (x1 - x0), 0.08 * (y1 - y0))
    top = min(0.035 * (y1 - y0) + 0.18, 0.18 * (y1 - y0))
    w, h = x1 - x0 - 2 * pad, y1 - y0 - top - pad
    n = cnt[key_]
    ncol = max(1, round(math.sqrt(n * w / h)))
    nrow = math.ceil(n / ncol)
    ncol = math.ceil(n / nrow)
    return {"px": x0 + pad, "py": y0 + top, "dx": w / ncol, "dy": h / nrow, "ncol": ncol, "nrow": nrow}


def alias_gids():
    from rdkit import Chem
    names = [("乙烯", "C=C"), ("丙烯", "C=CC"), ("1-丁烯", "C=CCC"), ("2-丁烯", "CC=CC"), ("异丁烯", "C=C(C)C"),
             ("1,3-丁二烯", "C=CC=C"), ("乙烷", "CC"), ("丙烷", "CCC"), ("正丁烷", "CCCC"), ("异丁烷", "CC(C)C"),
             ("异辛烷", "CC(C)CC(C)(C)C"), ("正辛烷", "CCCCCCCC"), ("正癸烷", "CCCCCCCCCC"), ("正十二烷", "CCCCCCCCCCCC"),
             ("正十四烷", "CCCCCCCCCCCCCC"), ("正十六烷", "CCCCCCCCCCCCCCCC"), ("2-甲基十五烷", "CCCCCCCCCCCCCC(C)C"),
             ("乙醇", "CCO"), ("正丙醇", "CCCO"), ("异丙醇", "CC(C)O"), ("乙醛", "CC=O"), ("丙酮", "CC(C)=O"), ("乙酸", "CC(=O)O"),
             ("二甲醚 DME", "COC"), ("苯", "c1ccccc1"), ("甲苯", "Cc1ccccc1"), ("邻二甲苯", "Cc1ccccc1C"), ("间二甲苯", "Cc1cccc(C)c1"),
             ("对二甲苯", "Cc1ccc(C)cc1"), ("乙苯", "CCc1ccccc1"), ("均三甲苯", "Cc1cc(C)cc(C)c1"), ("偏三甲苯", "Cc1ccc(C)c(C)c1"),
             ("均四甲苯", "Cc1cc(C)c(C)cc1C"), ("五甲苯", "Cc1cc(C)c(C)c(C)c1C"), ("六甲苯", "Cc1c(C)c(C)c(C)c(C)c1C"),
             ("异丙苯", "CC(C)c1ccccc1"), ("萘", "c1ccc2ccccc2c1"), ("2-甲基萘", "Cc1ccc2ccccc2c1"), ("2,6-二甲基萘", "Cc1ccc2cc(C)ccc2c1"),
             ("茚满", "c1ccc2c(c1)CCC2"), ("四氢萘", "c1ccc2c(c1)CCCC2"), ("正癸基苯", "CCCCCCCCCCc1ccccc1"),
             ("环己烷", "C1CCCCC1"), ("甲基环己烷", "CC1CCCCC1"), ("环己烯", "C1=CCCCC1"), ("1,3-环己二烯", "C1=CCCC=C1"),
             ("异丙基（R⁺ / R* / R•）", "C[CH]C"), ("叔丁基（R⁺ / R* / R•）", "C[C](C)C"), ("乙基（R* / R•）", "[CH2]C"),
             ("乙酰基", "C[C]=O")]
    want = {Chem.MolToSmiles(Chem.MolFromSmiles(s)): nm for nm, s in names}
    out = {}
    for g in range(len(special), NS):
        if g_lab[g] in want:
            out.setdefault(want[g_lab[g]], g)
    for i, nm in enumerate(spec_label):
        if nm in CORE_SHORT:
            out.setdefault("SMS 表面甲氧基 / CH₃* / •CH₃", i)
        out.setdefault(CORE_SHORT.get(nm, nm), i)
    return out


def combo_cn(f):
    return " ＋ ".join(BASE_CN.get(b, b) for b in f.split(";"))


fam_groups = []
for gname, bases in GROUPS:
    lst = [f for f in FAMS if f.split(";")[0] in bases]
    if lst:
        fam_groups.append([gname, lst])
alias = alias_gids()
pxy = alias.get("对二甲苯", -1)
total = sum(cnt.values()) + len(special)
title = (f"CO₂ 一步热催化加氢 · 物种级网络 C1–C{CMAX}（{'全量' if NS > 1_000_000 else '推荐版'}，一个物种一个节点，不含离子）")
meta = {
    "cmax": CMAX, "full_upto": CMAX, "W": Wd, "H": Hd, "chunk": CH, "nchunks": nchunks, "n_species": NS, "n_reactions": NR,
    "n_pool": int(len(POOL)), "classes": CLASSES, "fams": FAMS,
    "rows": [{"key": k, "label": lb, "block": blk, "cls": cl, "y0": rowy[r][0], "y1": rowy[r][1],
              "n": sum(cnt.get((r, c), 0) for c in range(2, CMAX + 1))} for r, (k, cl, lb, blk) in enumerate(ROWS)],
    "cols": [{"c": c, "x0": colx[c][0], "x1": colx[c][1], "label": ("C0–C1" if c == 1 else None),
              "n": sum(cnt.get((r, c), 0) for r in range(len(ROWS))) if c > 1 else len(special)} for c in range(1, CMAX + 1)],
    "cells": [dict(r=r, c=c, g0=cell_base[(r, c)], n=cnt[(r, c)], **geom((r, c))) for (r, c) in cell_list],
    "alias": alias,
    "special": [{"g": i, "label": CORE_SHORT.get(spec_label[i], spec_label[i]), "cls": cls_[t], "c": cn[t],
                 "smiles": smi[t], "xy": CORE_POS.get(spec_label[i])} for i, t in enumerate(special)],
    "pool_gids": [int(v) for v in POOL], "co2_gid": co2,
    "fam_counts": dict(Counter(FAMS[f] for f in rf.tolist())),
    "flows": flow_list[:6000], "max_level": int(lev.max()),
    "sizes_MB": {k: round(v / 1e6, 1) for k, v in sizes.items()},
    "start_gid": co2, "default_gid": pxy,
    "ui": {"title": title,
           "subtitle": ("R⁺（酸位 R–O–Z）、R*（金属位）、R•（气相）合并为同一个 R 基团节点，共振式合并，所有反应元素守恒。"
                        "缩放分级加载；单击物种看全部反应和从 CO₂ 出发的最短路径；左侧“子网络提取”可指定起点、终点和多走的步数并导出。"),
           "class_cn": CLASS_CN, "class_col": CLASS_COL, "fam_cn": {f: combo_cn(f) for f in FAMS},
           "fam_col": {f: BASE_COL.get(f.split(";")[0], "#888") for f in FAMS}, "fam_groups": fam_groups,
           "block_cn": {"MET": "含氧物种", "GAS": "R 基团", "OX": "烷烃 / 烯烃", "BA": "R 基团（烯基 / 环状）", "BA2": "多烯 / 环", "AR": "芳烃"},
           "col_prefix": "C", "start_label": "CO₂"},
}
(dst / "meta.js").write_text("window.__lodMeta = " + json.dumps(meta, ensure_ascii=False, separators=(",", ":")) + ";")
page = (HERE / "species_index.html").read_text(encoding="utf-8")
page = page.replace("<title>CO₂ 加氢物种级网络</title>", f"<title>CO₂ 加氢物种级网络（{'全量 C16' if NS > 1_000_000 else '推荐版'}）</title>")
(dst / "index.html").write_text(page, encoding="utf-8")
shutil.copy(HERE / "vendor" / "smiles-drawer.min.js", dst / "smiles-drawer.min.js")
shutil.copy(HERE / "vendor" / "smiles-drawer-LICENSE.md", dst / "smiles-drawer-LICENSE.md")
tpl = (HERE.parent / "scripts" / "subnet_view.html").read_text(encoding="utf-8")
tpl = tpl.replace("/*__SMILES_DRAWER__*/", (HERE / "vendor" / "smiles-drawer.min.js").read_text(encoding="utf-8"))
(dst / "subnet_tpl.js").write_text("window.__SUBNET_TPL = " + json.dumps(tpl, ensure_ascii=False) + ";", encoding="utf-8")
log("done", meta["sizes_MB"], "chunks", nchunks, "total species", total)
