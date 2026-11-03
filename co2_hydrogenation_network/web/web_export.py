"""Grid layout + compact JSON export of the isomer-resolved network for the interactive page.

python web_export.py ../data/isomer_graph_C12.json.gz network.json

Layout: columns = carbon number (C1 = hand-placed surface core), rows = species class
(same grid as deck/c10_layout.py). Inside a cell, species sorted by carbon skeleton.
Arrays are base64 little-endian typed arrays to keep the file small and fast to parse.
"""
import base64
import gzip
import json
import math
import sys
from collections import Counter, defaultdict

import numpy as np

src = sys.argv[1] if len(sys.argv) > 1 else "../data/isomer_graph_C12.json.gz"
dst = sys.argv[2] if len(sys.argv) > 2 else "network.json"
G = json.load(gzip.open(src, "rt"))
CMAX = G["cmax"]
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
           "ion_alkenyl", "naphthene", "cycloalkene", "ion_cyc", "cyclodiene", "aromatic", "core", "hub"]
FAMS = ["Z:prot", "Z:hshift", "Z:alkylshift", "Z:pcp", "Z:beta", "Z:meth", "Z:cyclise", "Z:HT", "H:hyd",
        "M:betaH", "M:CH2", "M:CO", "M:hyd", "M:hub", "O:dehyd", "O:hyd", "O:ox", "R:abs", "R:beta",
        "A:meth", "A:iso", "A:FC", "A:hub", "C1"]
cls2row = {c: i for i, r in enumerate(ROWS) for c in r[1]}

nodes = {n[0]: n for n in G["nodes"]}
cell = defaultdict(list)
for nid, cls, nc, sk in G["nodes"]:
    if cls != "hub" and nc > 1:
        cell[(cls2row[cls], nc)].append(nid)
cnt = {k: len(v) for k, v in cell.items()}
colmax = {c: max([cnt.get((r, c), 0) for r in range(len(ROWS))] + [1]) for c in range(2, CMAX + 1)}
rowmax = {r: max([cnt.get((r, c), 0) for c in range(2, CMAX + 1)] + [1]) for r in range(len(ROWS))}
C1W = 6.5
wraw = {c: max(1.6, math.sqrt(colmax[c])) for c in range(2, CMAX + 1)}
sc = (W - C1W) / sum(wraw.values())
colx, x = {1: (0.0, C1W)}, C1W
for c in range(2, CMAX + 1):
    colx[c] = (x, x + wraw[c] * sc)
    x += wraw[c] * sc
# row heights: same unit area as columns so cells stay roughly square-packed
hraw = {r: max(18.0, math.sqrt(rowmax[r])) for r in range(len(ROWS))}
H = round(W * 0.78, 3)
sr = H / sum(hraw.values())
rowy, y = {}, 0.0
for r in range(len(ROWS)):
    rowy[r] = (y, y + hraw[r] * sr)
    y += hraw[r] * sr

pos, size = {}, {}
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
    ids.sort(key=lambda i: (ROWS[r][1].index(nodes[i][1]), nodes[i][3].count("("), nodes[i][3], i))
    for k, i in enumerate(ids):
        rr, cc = divmod(k, ncol)
        pos[i] = (x0 + pad + (cc + 0.5) * dx, y0 + top + (rr + 0.5) * dy)
        size[i] = min(dx, dy)


def ry(key, f):
    r = [k for k, *_ in ROWS].index(key)
    return rowy[r][0] + f * (rowy[r][1] - rowy[r][0])


fx = lambda v: v * 0.78 / 3.3 * C1W  # noqa: E731  (fractions of the slide's C1 panel)
CORE = {
    "core|CO2": (fx(2.1), 0.012 * H, "CO₂"), "core|CO2d": (fx(1.0), ry("alc", .30), "CO₂^δ−*"),
    "core|HCOO": (fx(1.0), ry("alc", .62), "HCOO*"), "core|H2COO": (fx(1.0), ry("alc", .90), "H₂COO*"),
    "core|H2CO": (fx(1.0), ry("cbo", .30), "H₂CO*"), "core|CH3O": (fx(1.0), ry("cbo", .62), "CH₃O*"),
    "alcohol|CO": (fx(1.0), ry("cbo", .92), "CH₃OH"), "aldehyde|C=O": (fx(2.1), ry("cbo", .45), "HCHO"),
    "acid|O=CO": (fx(2.1), ry("alc", .75), "HCOOH"), "core|CO": (fx(3.2), ry("alc", .30), "CO*"),
    "hub|CO*": (fx(3.2), ry("alc", .62), "CO* 插入"), "core|HCO": (fx(3.2), ry("cbo", .20), "HCO*"),
    "core|HCOH": (fx(3.2), ry("cbo", .55), "HCOH*"), "core|CH": (fx(3.2), ry("cbo", .88), "CH*/C*"),
    "hub|CH2*": (fx(3.2), ry("ads", .35), "CH₂* 插入"), "ads|[CH3+]": (fx(2.1), ry("ads", .65), "CH₃*"),
    "rad|[CH3+]": (fx(2.1), ry("rad", .5), "·CH₃"), "alkane|C": (fx(2.1), ry("alk", .5), "CH₄"),
    "hub|CH3OH/SMS": (fx(1.5), ry("ia", .5), "CH₃OH/SMS"), "core|DME": (fx(3.0), ry("ia", .5), "DME"),
    "core|HCP": (fx(2.1), ry("ic", .5), "烃池 HCP"),
}
CORE_E = [("core|CO2", "core|CO2d"), ("core|CO2d", "core|HCOO"), ("core|HCOO", "core|H2COO"), ("core|H2COO", "core|H2CO"),
          ("core|H2CO", "core|CH3O"), ("core|CH3O", "alcohol|CO"), ("core|HCOO", "acid|O=CO"), ("core|H2CO", "aldehyde|C=O"),
          ("core|CO2", "core|CO"), ("core|CO", "hub|CO*"), ("core|CO", "core|HCO"), ("core|HCO", "core|HCOH"),
          ("core|HCOH", "core|CH"), ("core|CH", "hub|CH2*"), ("hub|CH2*", "ads|[CH3+]"), ("ads|[CH3+]", "alkane|C"),
          ("rad|[CH3+]", "alkane|C"), ("core|HCO", "aldehyde|C=O"), ("alcohol|CO", "hub|CH3OH/SMS"),
          ("hub|CH3OH/SMS", "core|DME"), ("hub|CH3OH/SMS", "core|HCP")]
labels = {}
for k, (x_, y_, lab) in CORE.items():
    pos[k] = (x_, y_)
    size[k] = 1.2
    labels[k] = lab
    if k not in nodes:
        nodes[k] = [k, "core", 1, k]

ids = sorted(pos, key=lambda i: (nodes[i][2], i))
idx = {k: i for i, k in enumerate(ids)}


def label(i):
    if i in labels:
        return labels[i]
    cls, s = nodes[i][1], i.split("|", 1)[1]
    if cls in ("rad", "ads"):
        s = s.replace("+]", "•]" if cls == "rad" else "*]")
    return s


edges = [(idx[u], idx[v], FAMS.index(f)) for u, v, f in G["edges"] if u in idx and v in idx]
edges += [(idx[u], idx[v], FAMS.index("C1")) for u, v in CORE_E]
E = np.array(edges, dtype=np.uint32)
b64 = lambda a: base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()  # noqa: E731
out = {
    "cmax": CMAX, "W": W, "H": H,
    "rows": [{"key": k, "label": lab, "block": blk, "cls": cl, "y0": rowy[r][0], "y1": rowy[r][1],
              "n": sum(cnt.get((r, c), 0) for c in range(2, CMAX + 1))} for r, (k, cl, lab, blk) in enumerate(ROWS)],
    "cols": [{"c": c, "x0": colx[c][0], "x1": colx[c][1],
              "n": sum(cnt.get((r, c), 0) for r in range(len(ROWS))) if c > 1 else len(CORE)} for c in range(1, CMAX + 1)],
    "cells": [[r, c, n] for (r, c), n in sorted(cnt.items())],
    "classes": CLASSES, "fams": FAMS,
    "fam_counts": Counter(FAMS[f] for *_, f in edges),
    "class_counts": Counter(nodes[i][1] for i in ids),
    "n_nodes": len(ids), "n_edges": len(edges),
    "label": "\n".join(label(i) for i in ids),
    "cls": b64(np.array([CLASSES.index(nodes[i][1]) for i in ids], dtype=np.uint8)),
    "cn": b64(np.array([nodes[i][2] for i in ids], dtype=np.uint8)),
    "xy": b64(np.array([pos[i] for i in ids], dtype=np.float32).ravel()),
    "sz": b64(np.array([size[i] for i in ids], dtype=np.float32)),
    "eu": b64(E[:, 0]), "ev": b64(E[:, 1]), "ef": b64(E[:, 2].astype(np.uint8)),
}
ALIAS = [("乙烯", "C=C"), ("丙烯", "C=CC"), ("1-丁烯", "C=CCC"), ("2-丁烯", "CC=CC"), ("异丁烯", "C=C(C)C"),
         ("1,3-丁二烯", "C=CC=C"), ("异戊二烯", "C=CC(C)=C"), ("乙烷", "CC"), ("丙烷", "CCC"), ("正丁烷", "CCCC"),
         ("异丁烷", "CC(C)C"), ("异戊烷", "CCC(C)C"), ("新戊烷", "CC(C)(C)C"), ("异辛烷", "CC(C)CC(C)(C)C"),
         ("正辛烷", "CCCCCCCC"), ("正癸烷", "CCCCCCCCCC"), ("正十二烷", "CCCCCCCCCCCC"),
         ("乙醇", "CCO"), ("正丙醇", "CCCO"), ("异丙醇", "CC(C)O"), ("正丁醇", "CCCCO"), ("乙醛", "CC=O"), ("丙醛", "CCC=O"),
         ("丙酮", "CC(C)=O"), ("乙酸", "CC(=O)O"), ("丙酸", "CCC(=O)O"),
         ("苯", "c1ccccc1"), ("甲苯", "Cc1ccccc1"), ("邻二甲苯", "Cc1ccccc1C"), ("间二甲苯", "Cc1cccc(C)c1"),
         ("对二甲苯", "Cc1ccc(C)cc1"), ("乙苯", "CCc1ccccc1"), ("均三甲苯", "Cc1cc(C)cc(C)c1"), ("偏三甲苯", "Cc1ccc(C)c(C)c1"),
         ("连三甲苯", "Cc1cccc(C)c1C"), ("均四甲苯", "Cc1cc(C)c(C)cc1C"), ("五甲苯", "Cc1cc(C)c(C)c(C)c1C"),
         ("六甲苯", "Cc1c(C)c(C)c(C)c(C)c1C"), ("异丙苯", "CC(C)c1ccccc1"), ("正丙苯", "CCCc1ccccc1"),
         ("对甲基异丙苯", "Cc1ccc(C(C)C)cc1"), ("萘", "c1ccc2ccccc2c1"), ("1-甲基萘", "Cc1cccc2ccccc12"),
         ("2-甲基萘", "Cc1ccc2ccccc2c1"), ("茚满", "c1ccc2c(c1)CCC2"), ("四氢萘", "c1ccc2c(c1)CCCC2"),
         ("环戊烷", "C1CCCC1"), ("环己烷", "C1CCCCC1"), ("甲基环己烷", "CC1CCCCC1"), ("环己烯", "C1=CCCCC1"),
         ("1,3-环己二烯", "C1=CCCC=C1"), ("1,4-环己二烯", "C1=CCC=CC1"), ("环戊二烯", "C1=CCC=C1"),
         ("异丙基正离子", "C[CH+]C"), ("叔丁基正离子", "C[C+](C)C"), ("仲丁基正离子", "CC[CH+]C"),
         ("叔戊基正离子", "CC[C+](C)C"), ("环己基正离子", "[CH+]1CCCCC1"),
         ("甲基环戊基正离子", "C[C+]1CCCC1"), ("甲基环己基正离子", "C[C+]1CCCCC1")]
from rdkit import Chem  # noqa: E402
lab_set = {}
for i in ids:
    lab_set.setdefault(label(i), i)
alias = {}
for name, smi in ALIAS:
    c = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    if c in lab_set:
        alias[c] = name
out["alias"] = alias
print("aliases matched", len(alias), "of", len(ALIAS), [n for n, s in ALIAS if Chem.MolToSmiles(Chem.MolFromSmiles(s)) not in lab_set])
json.dump(out, open(dst, "w"), ensure_ascii=False, separators=(",", ":"))
print(dst, "nodes", len(ids), "edges", len(edges), "bytes", len(json.dumps(out, ensure_ascii=False).encode()))
