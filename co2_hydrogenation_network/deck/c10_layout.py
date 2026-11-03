"""Grid layout + render of the full C1-C10 isomer network (one slide).

columns = carbon number (C1 column = hand-placed C1 surface core), rows = species class.
Writes c10_net.png / c10_net.pdf (plot area only) and c10_layout.json (boxes for native labels).
"""
import gzip
import json
import math
import sys
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

G = json.load(gzip.open(sys.argv[1] if len(sys.argv) > 1 else "../g10.json.gz", "rt"))
W, H = 41.1, 22.0          # plot area (in), placed at slide (2.6, 1.55)
DPI = int(sys.argv[2]) if len(sys.argv) > 2 else 200

ROWS = [  # key, classes, label, block
    ("alc", ["alcohol"], "醇 ROH / 烷氧基 RO*", "MET"),
    ("cbo", ["aldehyde", "ketone", "acid"], "醛（酰基 RCO*）/ 酮 / 羧酸", "MET"),
    ("ads", ["ads"], "吸附烷基 R*（金属/碳化物）", "MET"),
    ("rad", ["rad"], "气相烷基自由基 R•", "GAS"),
    ("alk", ["alkane"], "烷烃", "OX"),
    ("ole", ["alkene"], "烯烃", "OX"),
    ("ia", ["ion_alkyl"], "烷基碳正离子 R⁺", "BA"),
    ("die", ["diene"], "二烯", "BA2"),
    ("ie", ["ion_alkenyl"], "烯基 / 烯丙基碳正离子", "BA"),
    ("cyc", ["naphthene", "cycloalkene"], "环烷 / 环烯", "BA2"),
    ("ic", ["ion_cyc"], "环状碳正离子", "BA"),
    ("cd", ["cyclodiene"], "环二烯", "BA2"),
    ("aro", ["aromatic"], "芳烃", "AR"),
]
BLOCK = {"MET": "#FBE3E1", "GAS": "#ECECEA", "OX": "#FBF1DC", "BA": "#DDE9F8", "BA2": "#EAF1FB", "AR": "#DFF1E5"}
NODE_COL = {
    "alcohol": "#8E5BB5", "aldehyde": "#B5487F", "ketone": "#C77A2E", "acid": "#6A3D9A",
    "ads": "#D64545", "rad": "#6E6E6E", "alkane": "#3A3A3A", "alkene": "#E07B00",
    "ion_alkyl": "#1F5FBF", "diene": "#C9A100", "ion_alkenyl": "#3D8FE0", "naphthene": "#5C6B73",
    "cycloalkene": "#2F9C95", "ion_cyc": "#6C4FD1", "cyclodiene": "#0E7C66", "aromatic": "#1E8C3A",
    "hub": "#111111", "core": "#111111",
}
EDGE_COL = {"Z": "#2A6FD6", "H": "#D99A00", "M": "#D64545", "O": "#8E5BB5", "R": "#777777", "A": "#1E9C4B"}
EDGE_SPECIAL = {"Z:HT": "#13A39A", "Z:prot": "#5B9BE8"}

cls2row = {c: i for i, r in enumerate(ROWS) for c in r[1]}
nodes = {n[0]: n for n in G["nodes"]}
cell = defaultdict(list)
for nid, cls, nc, sk in G["nodes"]:
    if cls == "hub" or nc == 1:
        continue
    cell[(cls2row[cls], nc)].append(nid)
NCOL = G["cmax"]
# --- column widths / row heights
colmax = {c: max([len(cell[(r, c)]) for r in range(len(ROWS))] + [1]) for c in range(2, NCOL + 1)}
rowmax = {r: max([len(cell[(r, c)]) for c in range(2, NCOL + 1)] + [1]) for r in range(len(ROWS))}
C1W = 3.3
wraw = {c: max(1.6, math.sqrt(colmax[c])) for c in range(2, NCOL + 1)}
s = (W - C1W) / sum(wraw.values())
colx, x = {1: (0.0, C1W)}, C1W
for c in range(2, NCOL + 1):
    colx[c] = (x, x + wraw[c] * s)
    x += wraw[c] * s
hraw = {r: max(9.0, math.sqrt(rowmax[r])) for r in range(len(ROWS))}
s = H / sum(hraw.values())
rowy, y = {}, 0.0
for r in range(len(ROWS)):
    rowy[r] = (y, y + hraw[r] * s)
    y += hraw[r] * s

# --- place nodes inside cells
pos, size = {}, {}
nbranch = lambda sk: sk.count("(")  # noqa: E731
for (r, c), ids in cell.items():
    if not ids:
        continue
    x0, x1 = colx[c]
    y0, y1 = rowy[r]
    px, py = 0.10, 0.16
    w, h = x1 - x0 - 2 * px, y1 - y0 - py - 0.08
    n = len(ids)
    nc_ = max(1, round(math.sqrt(n * w / h)))
    nr_ = math.ceil(n / nc_)
    nc_ = math.ceil(n / nr_)
    dx, dy = w / nc_, h / nr_
    ids.sort(key=lambda i: (ROWS[r][1].index(nodes[i][1]), nbranch(nodes[i][3]), nodes[i][3], i))
    for k, i in enumerate(ids):
        rr, cc = divmod(k, nc_)
        pos[i] = (x0 + px + (cc + 0.5) * dx, y0 + py + (rr + 0.5) * dy)
        size[i] = min(dx, dy)

# --- C1 core panel (hand placed), relative to column 1
def ry(rkey, f=0.5):
    r = [k for k, *_ in ROWS].index(rkey)
    return rowy[r][0] + f * (rowy[r][1] - rowy[r][0])

CORE = {  # id: (x, y, label)
    "core|CO2": (2.1, 0.25, r"$\mathrm{CO_2}$"),
    "core|CO2d": (1.0, ry("alc", 0.30), r"$\mathrm{CO_2^{\delta-}*}$"),
    "core|HCOO": (1.0, ry("alc", 0.62), r"$\mathrm{HCOO*}$"),
    "core|H2COO": (1.0, ry("alc", 0.90), r"$\mathrm{H_2COO*}$"),
    "core|H2CO": (1.0, ry("cbo", 0.30), r"$\mathrm{H_2CO*}$"),
    "core|CH3O": (1.0, ry("cbo", 0.62), r"$\mathrm{CH_3O*}$"),
    "alcohol|CO": (1.0, ry("cbo", 0.92), r"$\mathrm{CH_3OH}$"),
    "aldehyde|C=O": (2.1, ry("cbo", 0.45), r"$\mathrm{HCHO}$"),
    "acid|O=CO": (2.1, ry("alc", 0.75), r"$\mathrm{HCOOH}$"),
    "core|CO": (3.2, ry("alc", 0.30), r"$\mathrm{CO*}$"),
    "hub|CO*": (3.2, ry("alc", 0.62), r"$\mathrm{CO*}$ 插入"),
    "core|HCO": (3.2, ry("cbo", 0.20), r"$\mathrm{HCO*}$"),
    "core|HCOH": (3.2, ry("cbo", 0.55), r"$\mathrm{HCOH*}$"),
    "core|CH": (3.2, ry("cbo", 0.88), r"$\mathrm{CH*/C*}$"),
    "hub|CH2*": (3.2, ry("ads", 0.35), r"$\mathrm{CH_2*}$"),
    "ads|[CH3+]": (2.1, ry("ads", 0.6), r"$\mathrm{CH_3*}$"),
    "rad|[CH3+]": (2.1, ry("rad", 0.5), r"$\mathrm{\cdot CH_3}$"),
    "alkane|C": (2.1, ry("alk", 0.5), r"$\mathrm{CH_4}$"),
    "hub|CH3OH/SMS": (1.4, ry("ia", 0.5), r"$\mathrm{CH_3OH / SMS}$"),
    "core|DME": (3.0, ry("ia", 0.5), r"DME"),
    "core|HCP": (2.1, ry("ic", 0.5), "烃池 HCP"),
}
CORE_E = [("core|CO2", "core|CO2d", "O"), ("core|CO2d", "core|HCOO", "O"), ("core|HCOO", "core|H2COO", "O"),
          ("core|H2COO", "core|H2CO", "O"), ("core|H2CO", "core|CH3O", "O"), ("core|CH3O", "alcohol|CO", "O"),
          ("core|HCOO", "acid|O=CO", "O"), ("core|H2CO", "aldehyde|C=O", "O"), ("core|CO2", "core|CO", "M"),
          ("core|CO", "hub|CO*", "M"), ("core|CO", "core|HCO", "M"), ("core|HCO", "core|HCOH", "M"),
          ("core|HCOH", "core|CH", "M"), ("core|CH", "hub|CH2*", "M"), ("hub|CH2*", "ads|[CH3+]", "M"),
          ("ads|[CH3+]", "alkane|C", "M"), ("rad|[CH3+]", "alkane|C", "R"), ("core|HCO", "aldehyde|C=O", "M"),
          ("alcohol|CO", "hub|CH3OH/SMS", "Z"), ("hub|CH3OH/SMS", "core|DME", "Z"), ("hub|CH3OH/SMS", "core|HCP", "Z"),
          ("core|HCP", "hub|CH3OH/SMS", "Z")]
for k, (x_, y_, _) in CORE.items():
    pos[k] = (x_ * 0.78, y_)
    size[k] = 0.5
# unplaced C1 graph nodes (e.g. C1 ketone/rad duplicates) -> drop
EDGES = [e for e in G["edges"] if e[0] in pos and e[1] in pos]

# --- draw
fig = plt.figure(figsize=(W, H), dpi=DPI)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.axis("off")
plt.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]
plt.rcParams["mathtext.fontset"] = "dejavusans"
# row background blocks
for r, (key, cl, lab, blk) in enumerate(ROWS):
    y0, y1 = rowy[r]
    ax.add_patch(Rectangle((C1W, y0 + 0.02), W - C1W, y1 - y0 - 0.04, fc=BLOCK[blk], ec="none", zorder=0))
ax.add_patch(FancyBboxPatch((0.06, 0.04), C1W - 0.14, H - 0.08, boxstyle="round,pad=0,rounding_size=0.15",
                            fc="#F7F6F2", ec="#BDBAB2", lw=1.2, zorder=0))
for c in range(3, NCOL + 1):
    ax.plot([colx[c][0]] * 2, [0, H], color="white", lw=2.2, zorder=0.5)
# cell counts
for (r, c), ids in cell.items():
    if not ids:
        continue
    ax.text(colx[c][0] + 0.08, rowy[r][0] + 0.05, f"{len(ids)}", fontsize=7.5, color="#55534E",
            va="top", ha="left", zorder=5)
# edges
segs = defaultdict(list)
for u, v, fam in EDGES:
    segs[fam].append((pos[u], pos[v]))
for u, v, cat in CORE_E:
    segs["CORE:" + cat].append((pos[u], pos[v]))
for fam, ss in segs.items():
    cat = fam.split(":")[0] if not fam.startswith("CORE") else fam.split(":")[1]
    col = EDGE_SPECIAL.get(fam, EDGE_COL[cat])
    hubfam = fam in ("M:hub", "A:hub", "Z:meth")
    alpha = 0.9 if fam.startswith("CORE") else 0.035 if hubfam else 0.22
    lw = 1.4 if fam.startswith("CORE") else 0.35
    ax.add_collection(LineCollection(ss, colors=col, linewidths=lw, alpha=alpha, zorder=1 if not fam.startswith("CORE") else 3))
# nodes
xs, ys, ss_, cs = [], [], [], []
for i, (x_, y_) in pos.items():
    if i in CORE:
        continue
    d = min(0.16, max(0.035, 0.5 * size[i]))
    xs.append(x_)
    ys.append(y_)
    ss_.append((d * 72) ** 2)
    cs.append(NODE_COL[nodes[i][1]])
ax.scatter(xs, ys, s=ss_, c=cs, linewidths=0, zorder=4)
# SMILES labels (tiny; zoom in the PDF)
for i, (x_, y_) in pos.items():
    if i in CORE:
        continue
    fs = min(6.5, max(1.2, size[i] * 72 * 0.17))
    smi = i.split("|", 1)[1]
    if nodes[i][1] in ("rad", "ads"):
        smi = smi.replace("+]", "•]" if nodes[i][1] == "rad" else "*]")
    ax.text(x_, y_ + min(0.16, max(0.035, 0.5 * size[i])) * 0.62, smi, fontsize=fs, ha="center", va="top",
            color="#4A4945", zorder=4, family="DejaVu Sans Mono")
# core nodes
for k, (x_, y_, lab) in CORE.items():
    ishub = k.startswith("hub")
    ax.text(x_ * 0.78, y_, lab, fontsize=11 if ishub else 10, ha="center", va="center", color="white" if ishub else "#111",
            zorder=6, bbox=dict(boxstyle="round,pad=0.3,rounding_size=0.25", fc="#2B2B2B" if ishub else "white",
                                ec="#2B2B2B", lw=1.1))
fig.savefig("c10_net.png", dpi=DPI, facecolor="white")
from PIL import Image  # noqa: E402
Image.MAX_IMAGE_PIXELS = None
Image.open("c10_net.png").convert("RGB").save("c10_net.jpg", quality=90, optimize=True)
if "--pdf" in sys.argv:
    fig.savefig("c10_net.pdf", facecolor="white")
json.dump({"W": W, "H": H, "colx": colx, "rowy": {str(k): v for k, v in rowy.items()},
           "rows": [[k, lab, blk, sum(len(cell[(r, c)]) for c in range(2, NCOL + 1))] for r, (k, cl, lab, blk) in enumerate(ROWS)],
           "cols": {c: sum(len(cell[(r, c)]) for r in range(len(ROWS))) for c in range(2, NCOL + 1)},
           "n_nodes": len(pos), "n_edges": len(EDGES) + len(CORE_E),
           "edge_fams": Counter(e[2] for e in EDGES)},
          open("c10_layout.json", "w"), ensure_ascii=False, indent=1)
print("nodes", len(pos), "edges", len(EDGES) + len(CORE_E))
