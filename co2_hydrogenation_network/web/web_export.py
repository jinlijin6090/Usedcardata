"""Grid layout + compact JSON export of a generated network for the interactive page.

    python web_export.py ../data/net_C16_shape .        # -> nodes.json, reactions.json

Input: species.tsv.gz / reactions.tsv.gz / stats.json written by scripts/network_gen.py.
Species ids are kept (array index = id in species.tsv.gz), so the page and the tables agree.
Layout: columns = carbon number (C1 = hand-placed surface core), rows = species class.
Arrays are base64 little-endian typed arrays.
"""
import base64
import gzip
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem

src = Path(sys.argv[1] if len(sys.argv) > 1 else "../data/net_C16_shape")
dst = Path(sys.argv[2] if len(sys.argv) > 2 else ".")
stats = json.loads((src / "stats.json").read_text())
CMAX, FULL = stats["cmax"], stats["full_upto"]
W = 100.0

ROWS = [
    ("alc", ["alcohol"], "醇 / 烷氧基 RO*", "MET"),
    ("cbo", ["aldehyde", "ketone", "acid"], "醛（酰基 RCO*）/ 酮 / 羧酸", "MET"),
    ("ads", ["ads"], "吸附烷基 R*", "MET"),
    ("rad", ["rad"], "气相自由基 R•", "GAS"),
    ("alk", ["alkane"], "烷烃", "OX"),
    ("ole", ["alkene"], "烯烃", "OX"),
    ("ia", ["ion_alkyl"], "烷基碳正离子", "BA"),
    ("die", ["diene"], "二烯", "BA2"),
    ("ie", ["ion_alkenyl"], "烯基 / 烯丙基碳正离子", "BA"),
    ("cyc", ["naphthene", "cycloalkene"], "环烷 / 环烯", "BA2"),
    ("ic", ["ion_cyc"], "环状碳正离子", "BA"),
    ("cd", ["cyclodiene"], "环二烯", "BA2"),
    ("aro", ["aromatic"], "芳烃", "AR"),
]
CLASSES = ["alcohol", "aldehyde", "ketone", "acid", "ads", "rad", "alkane", "alkene", "ion_alkyl", "diene",
           "ion_alkenyl", "naphthene", "cycloalkene", "ion_cyc", "cyclodiene", "aromatic", "core", "pool"]
cls2row = {c: i for i, r in enumerate(ROWS) for c in r[1]}

# ---- species
sp = []  # (id, cls, C, smiles, skeleton)
with gzip.open(src / "species.tsv.gz", "rt") as fh:
    next(fh)
    for line in fh:
        i, c, n, s, sk = line.rstrip("\n").split("\t")
        sp.append((int(i), c, int(n), s, sk))
assert [x[0] for x in sp] == list(range(len(sp)))
N = len(sp)
key2id = {f"{c}|{s}": i for i, c, n, s, sk in sp if c not in ("core", "pool")}

cell = defaultdict(list)
for i, c, n, s, sk in sp:
    if c in ("core", "pool") or n <= 1:
        continue
    cell[(cls2row[c], n)].append(i)
cnt = {k: len(v) for k, v in cell.items()}
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

xy = np.full((N, 2), -1.0, dtype=np.float32)
sz = np.zeros(N, dtype=np.float32)
for (r, c), ids in cell.items():
    x0, x1 = colx[c]
    y0, y1 = rowy[r]
    pad = 0.012 * (x1 - x0) + 0.05
    top = 0.035 * (y1 - y0) + 0.18
    w, h = x1 - x0 - 2 * pad, y1 - y0 - top - pad
    n = len(ids)
    ncol = max(1, round(math.sqrt(n * w / h)))
    nrow = math.ceil(n / ncol)
    ncol = math.ceil(n / nrow)
    dx, dy = w / ncol, h / nrow
    ids.sort(key=lambda i: (ROWS[r][1].index(sp[i][1]), sp[i][4].count("("), sp[i][4], sp[i][3]))
    for k, i in enumerate(ids):
        rr, cc = divmod(k, ncol)
        xy[i] = (x0 + pad + (cc + 0.5) * dx, y0 + top + (rr + 0.5) * dy)
        sz[i] = min(dx, dy)


def ry(key, f):
    r = [k for k, *_ in ROWS].index(key)
    return rowy[r][0] + f * (rowy[r][1] - rowy[r][0])


fx = lambda v: v * 0.78 / 3.3 * C1W  # noqa: E731
CORE_POS = {
    "core|CO2": (fx(2.1), 0.012 * H), "core|CO2d": (fx(1.0), ry("alc", .30)), "core|HCOO": (fx(1.0), ry("alc", .62)),
    "core|H2COO": (fx(1.0), ry("alc", .90)), "core|H2CO": (fx(1.0), ry("cbo", .30)), "core|CH3O": (fx(1.0), ry("cbo", .62)),
    "alcohol|CO": (fx(1.0), ry("cbo", .92)), "aldehyde|C=O": (fx(2.1), ry("cbo", .45)), "acid|O=CO": (fx(2.1), ry("alc", .75)),
    "core|CO*": (fx(3.2), ry("alc", .40)), "core|HCO": (fx(3.2), ry("cbo", .20)), "core|HCOH": (fx(3.2), ry("cbo", .55)),
    "core|CH": (fx(3.2), ry("cbo", .88)), "core|CH2*": (fx(3.2), ry("ads", .35)), "ads|[CH3+]": (fx(2.1), ry("ads", .65)),
    "rad|[CH3+]": (fx(2.1), ry("rad", .5)), "alkane|C": (fx(2.1), ry("alk", .5)), "core|SMS": (fx(1.5), ry("ia", .5)),
    "core|DME": (fx(3.0), ry("ia", .5)), "core|HCP": (fx(2.1), ry("ic", .5)),
}
CORE_LABEL = {"alcohol|CO": "CH₃OH", "aldehyde|C=O": "HCHO", "acid|O=CO": "HCOOH", "ads|[CH3+]": "CH₃*",
              "rad|[CH3+]": "·CH₃", "alkane|C": "CH₄"}
labels = []
CORE_KEY_BY_LABEL = {"CO₂": "core|CO2", "CO₂^δ−*": "core|CO2d", "HCOO*": "core|HCOO", "H₂COO*": "core|H2COO",
                     "H₂CO*": "core|H2CO", "CH₃O*": "core|CH3O", "CO*": "core|CO*", "HCO*": "core|HCO", "HCOH*": "core|HCOH",
                     "CH*/C*": "core|CH", "CH₂*": "core|CH2*", "CH₃OH/SMS": "core|SMS", "DME": "core|DME", "烃池 HCP": "core|HCP"}
for i, c, n, s, sk in sp:
    if c in ("core", "pool"):
        labels.append(s)                      # species.tsv already stores the readable label for core/pool
        if c == "core":
            kk = CORE_KEY_BY_LABEL[s]
            xy[i] = CORE_POS[kk]
            sz[i] = 1.2
        continue
    k = f"{c}|{s}"
    if k in CORE_POS:
        xy[i] = CORE_POS[k]
        sz[i] = 1.2
    lab = CORE_LABEL.get(k, s)
    if c == "rad" and k not in CORE_LABEL:
        lab = s.replace("+]", "•]")
    if c == "ads" and k not in CORE_LABEL:
        lab = s.replace("+]", "*]")
    labels.append(lab)
assert len(labels) == N

# ---- reactions
FAMS = ["Z:prot", "Z:hshift", "Z:alkylshift", "Z:pcp", "Z:beta", "Z:meth", "Z:cyclise", "Z:HT", "H:hyd",
        "M:betaH", "M:CH2", "M:CO", "M:hyd", "O:dehyd", "O:hyd", "O:ox", "R:abs", "R:beta",
        "A:meth", "A:iso", "A:FC", "A:hcp", "C1", "O:OH", "A:cyc"]
NONE = 0xFFFFFFFF
rf, ra, rb, pa, pb = [], [], [], [], []
with gzip.open(src / "reactions.tsv.gz", "rt") as fh:
    next(fh)
    for line in fh:
        _, fam, r, p, _ = line.rstrip("\n").split("\t")
        rs, ps = [int(v) for v in r.split(",")], [int(v) for v in p.split(",")]
        assert len(rs) <= 2 and len(ps) <= 2
        rf.append(FAMS.index(fam))
        ra.append(rs[0]); rb.append(rs[1] if len(rs) > 1 else NONE)
        pa.append(ps[0]); pb.append(ps[1] if len(ps) > 1 else NONE)
R = len(rf)

ALIAS = [("乙烯", "C=C"), ("丙烯", "C=CC"), ("1-丁烯", "C=CCC"), ("2-丁烯", "CC=CC"), ("异丁烯", "C=C(C)C"),
         ("1,3-丁二烯", "C=CC=C"), ("异戊二烯", "C=CC(C)=C"), ("乙烷", "CC"), ("丙烷", "CCC"), ("正丁烷", "CCCC"),
         ("异丁烷", "CC(C)C"), ("异戊烷", "CCC(C)C"), ("新戊烷", "CC(C)(C)C"), ("异辛烷", "CC(C)CC(C)(C)C"),
         ("正辛烷", "CCCCCCCC"), ("正癸烷", "CCCCCCCCCC"), ("正十二烷", "CCCCCCCCCCCC"), ("正十四烷", "CCCCCCCCCCCCCC"),
         ("正十六烷", "CCCCCCCCCCCCCCCC"), ("2-甲基十五烷", "CCCCCCCCCCCCCC(C)C"),
         ("乙醇", "CCO"), ("正丙醇", "CCCO"), ("异丙醇", "CC(C)O"), ("正丁醇", "CCCCO"), ("乙醛", "CC=O"), ("丙醛", "CCC=O"),
         ("丙酮", "CC(C)=O"), ("乙酸", "CC(=O)O"), ("丙酸", "CCC(=O)O"),
         ("苯", "c1ccccc1"), ("甲苯", "Cc1ccccc1"), ("邻二甲苯", "Cc1ccccc1C"), ("间二甲苯", "Cc1cccc(C)c1"),
         ("对二甲苯", "Cc1ccc(C)cc1"), ("乙苯", "CCc1ccccc1"), ("均三甲苯", "Cc1cc(C)cc(C)c1"), ("偏三甲苯", "Cc1ccc(C)c(C)c1"),
         ("连三甲苯", "Cc1cccc(C)c1C"), ("均四甲苯", "Cc1cc(C)c(C)cc1C"), ("五甲苯", "Cc1cc(C)c(C)c(C)c1C"),
         ("六甲苯", "Cc1c(C)c(C)c(C)c(C)c1C"), ("异丙苯", "CC(C)c1ccccc1"), ("正丙苯", "CCCc1ccccc1"),
         ("对甲基异丙苯", "Cc1ccc(C(C)C)cc1"), ("萘", "c1ccc2ccccc2c1"), ("1-甲基萘", "Cc1cccc2ccccc12"),
         ("2-甲基萘", "Cc1ccc2ccccc2c1"), ("2,6-二甲基萘", "Cc1ccc2cc(C)ccc2c1"), ("茚满", "c1ccc2c(c1)CCC2"),
         ("四氢萘", "c1ccc2c(c1)CCCC2"), ("正己基苯", "CCCCCCc1ccccc1"), ("正癸基苯", "CCCCCCCCCCc1ccccc1"),
         ("环戊烷", "C1CCCC1"), ("环己烷", "C1CCCCC1"), ("甲基环己烷", "CC1CCCCC1"), ("环己烯", "C1=CCCCC1"),
         ("1,3-环己二烯", "C1=CCCC=C1"), ("1,4-环己二烯", "C1=CCC=CC1"), ("环戊二烯", "C1=CCC=C1"),
         ("异丙基正离子", "C[CH+]C"), ("叔丁基正离子", "C[C+](C)C"), ("仲丁基正离子", "CC[CH+]C"),
         ("叔戊基正离子", "CC[C+](C)C"), ("环己基正离子", "[CH+]1CCCCC1"),
         ("甲基环戊基正离子", "C[C+]1CCCC1"), ("甲基环己基正离子", "C[C+]1CCCCC1")]
alias, missing = {}, []
for name, smi in ALIAS:
    c = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    found = any(f"{k}|{c}" in key2id for k in ("alkane", "alkene", "diene", "aromatic", "naphthene", "cycloalkene",
                                                 "cyclodiene", "alcohol", "aldehyde", "ketone", "acid", "ion_alkyl", "ion_cyc"))
    if found:
        alias[c] = name
    else:
        missing.append(name)

b64 = lambda a: base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()  # noqa: E731
meta = {
    "cmax": CMAX, "full_upto": FULL, "W": W, "H": H,
    "rows": [{"key": k, "label": lab, "block": blk, "cls": cl, "y0": rowy[r][0], "y1": rowy[r][1],
              "n": sum(cnt.get((r, c), 0) for c in range(2, CMAX + 1))} for r, (k, cl, lab, blk) in enumerate(ROWS)],
    "cols": [{"c": c, "x0": colx[c][0], "x1": colx[c][1],
              "n": sum(cnt.get((r, c), 0) for r in range(len(ROWS))) if c > 1 else sum(1 for x in sp if x[1] == "core" or (x[2] == 1 and x[1] != "pool"))}
             for c in range(1, CMAX + 1)],
    "cells": [[r, c, n] for (r, c), n in sorted(cnt.items())],
    "classes": CLASSES, "fams": FAMS,
    "fam_counts": Counter(FAMS[f] for f in rf),
    "class_counts": Counter(x[1] for x in sp),
    "n_nodes": N, "n_reactions": R, "alias": alias,
    "label": "\n".join(labels),
    "cls": b64(np.array([CLASSES.index(x[1]) for x in sp], dtype=np.uint8)),
    "cn": b64(np.array([x[2] for x in sp], dtype=np.uint8)),
    "xy": b64(xy.ravel()), "sz": b64(sz),
}
rx = {"n_reactions": R, "rf": b64(np.array(rf, dtype=np.uint8)), "ra": b64(np.array(ra, dtype=np.uint32)),
      "rb": b64(np.array(rb, dtype=np.uint32)), "pa": b64(np.array(pa, dtype=np.uint32)), "pb": b64(np.array(pb, dtype=np.uint32))}
dst.mkdir(parents=True, exist_ok=True)
for name, obj in (("nodes.json", meta), ("reactions.json", rx)):
    (dst / name).write_text(json.dumps(obj, ensure_ascii=False, separators=(",", ":")))
    print(name, round((dst / name).stat().st_size / 1e6, 2), "MB")
print("species", N, "reactions", R, "aliases", len(alias), "missing", missing)
