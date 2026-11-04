"""Figure: DRM -> syngas -> FT -> upgrading vs direct CO2 hydrogenation share one reaction network.

    python route_figure.py [--pdf]   # -> ../figures/route_shared_network.png (/.pdf)

Inputs: ../web/nodes.json (layout of the recommended C16 network), ../data/net_C16_shape/reactions.tsv.gz,
../data/route_codes_shape.npy (+ _species.npy) and ../data/route_compare_{shape,full}.json from route_compare.py.
Fonts: Noto Sans SC / Noto Serif SC (SIL OFL) from FONT_DIR, DejaVu Sans as glyph fallback (subscripts).
"""
import base64
import gzip
import json
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
NET = HERE.parent
FONT_DIR = Path(os.environ.get("FONT_DIR", str(NET / "assets" / "fonts")))
for f in FONT_DIR.glob("*.ttf"):
    fm.fontManager.addfont(str(f))
SANS = ["Noto Sans SC", "DejaVu Sans"]
SERIF = ["Noto Serif SC", "Noto Sans SC", "DejaVu Sans"]
plt.rcParams["font.family"] = SANS

INK, INK2, MUTED, LINE, PAPER = "#16191F", "#3F4651", "#6B7380", "#D5DAE1", "#FFFFFF"
SHARED = "#2F6DB5"      # chemistry used by both routes
BONLY = "#E8590C"       # used only by direct CO2 hydrogenation (methanol entry + methylation)
AROUTE = "#13906F"      # route A accent (DRM / syngas / FT)
BAND = {"MET": "#FBEDEC", "GAS": "#F1F2F4", "OX": "#FBF5E8", "BA": "#EAF1FA", "BA2": "#F2F6FB", "AR": "#EAF5EE"}

# ------------------------------------------------------------------ data
D = json.loads((NET / "web" / "nodes.json").read_text())
b64 = lambda s, t: np.frombuffer(base64.b64decode(s), dtype=t)  # noqa: E731
xy = b64(D["xy"], np.float32).reshape(-1, 2)
cls = b64(D["cls"], np.uint8)
classes = D["classes"]
labels = D["label"].split("\n")
POOL, CORE = classes.index("pool"), classes.index("core")
code = np.load(NET / "data" / "route_codes_shape.npy")
sp_code = np.load(NET / "data" / "route_codes_shape_species.npy")
RS = json.loads((NET / "data" / "route_compare_shape.json").read_text())
RFp = NET / "data" / "route_compare_full.json"
RF = json.loads(RFp.read_text()) if RFp.exists() else None
ST = RF or RS

segs = {0: [], 1: []}
core_segs = {0: [], 1: []}
with gzip.open(NET / "data" / "net_C16_shape" / "reactions.tsv.gz", "rt") as fh:
    next(fh)
    for k, line in enumerate(fh):
        c = int(code[k])
        if c > 1:
            continue
        _, f, r, p, _ = line.split("\t", 4)
        for a in map(int, r.split(",")):
            for q in map(int, p.split(",")):
                if POOL in (cls[a], cls[q]) or a == q:
                    continue
                if f == "C1":
                    core_segs[c].append((a, q))
                elif CORE not in (cls[a], cls[q]):        # hub fans (CH2*, CO*, SMS) are shown by the C1 panel instead
                    segs[c].append((a, q))

# ------------------------------------------------------------------ canvas (inches, y down)
FW, FH = 20.0, 11.6
fig = plt.figure(figsize=(FW, FH), dpi=200)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, FW)
ax.set_ylim(FH, 0)
ax.axis("off")
fig.patch.set_facecolor(PAPER)

NX0, NY0, NW, NH = 5.25, 1.78, 14.30, 5.50
C1W = D["cols"][0]["x1"]
C1PX = 2.15                              # widened C1 column so the entry chemistry is readable
ky = NH / D["H"]
COLS = D["cols"][1:]
orig = np.array([c["x1"] - c["x0"] for c in COLS])
wt = np.maximum(orig * (NW - C1PX) / orig.sum(), 0.34)            # every carbon number gets a visible column
big = wt > 0.34
wt[big] *= ((NW - C1PX) - wt[~big].sum()) / wt[big].sum()
X0S = NX0 + C1PX + np.concatenate([[0], np.cumsum(wt)[:-1]])


def P(x, y):
    if x <= C1W:
        X = NX0 + x / C1W * C1PX
    else:
        j = min(len(COLS) - 1, int(np.searchsorted([c["x1"] for c in COLS], x)))
        c = COLS[j]
        X = X0S[j] + (x - c["x0"]) / (c["x1"] - c["x0"]) * wt[j]
    return X, NY0 + y * ky


def text(x, y, s, size=11, color=INK, weight="normal", ha="left", va="baseline", family=None, **kw):
    return ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va, fontfamily=family or SANS, **kw)


def card(x, y, w, h, fc="#FFFFFF", ec=LINE, lw=1.0, r=0.12, z=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, zorder=z))


def chip(x, y, s, color, size=8.4, z=4):
    t = text(x, y, s, size=size, color=color, va="center", zorder=z + 1)
    fig.canvas.draw_idle()
    w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi + 0.16
    ax.add_patch(FancyBboxPatch((x - 0.08, y - 0.12), w, 0.24, boxstyle="round,pad=0,rounding_size=0.11",
                                fc=color + "14", ec=color + "60", lw=0.7, zorder=z))
    return w


def badge(x, y, n, color):
    ax.add_patch(Circle((x, y), 0.16, fc=color, ec="none", zorder=4))
    text(x, y + 0.01, str(n), size=9, color="white", weight="bold", ha="center", va="center", zorder=5)


# ------------------------------------------------------------------ title
text(0.45, 0.66, "两条航煤路线，同一张反应网络", size=27, weight="bold", family=SERIF)
text(0.47, 1.12, "甲烷干重整 → 合成气费托 → 加氢提质，与 CO₂ 直接加氢制航煤：从 C2 起走的是同一批中间体，差别只在 C1 入口和链增长的单体",
     size=12.5, color=INK2)

# ------------------------------------------------------------------ network panel
for r in D["rows"]:
    x0, y0 = P(C1W, r["y0"])
    x1, y1 = P(D["W"], r["y1"])
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=BAND[r["block"]], ec="none", zorder=0))
ax.add_patch(FancyBboxPatch((NX0 + 0.03, NY0 + 0.02), C1PX - 0.1, NH - 0.04, boxstyle="round,pad=0,rounding_size=0.1",
                            fc="#F6F7F9", ec=LINE, lw=0.8, zorder=0))
for c in D["cols"][2:]:
    x, _ = P(c["x0"], 0)
    ax.plot([x, x], [NY0, NY0 + NH], color=PAPER, lw=1.1, zorder=0.5)
for c in D["cols"]:
    x0, _ = P(c["x0"], 0)
    x1, _ = P(c["x1"], 0)
    if x1 - x0 > 0.26:
        text((x0 + x1) / 2, NY0 - 0.07, "C1 入口" if c["c"] == 1 else f"C{c['c']}", size=8.6 if c["c"] > 1 else 9, color=MUTED if c["c"] > 1 else INK2,
             ha="center", weight="normal" if c["c"] > 1 else "bold")
jx0, _ = P(D["cols"][7]["x0"], 0)
jx1, _ = P(D["W"], 0)
ax.add_patch(Rectangle((jx0, NY0 - 0.30), jx1 - jx0, 0.055, fc=AROUTE, ec="none", zorder=1))
text(jx1, NY0 - 0.38, "航煤馏分 C8–C16", size=10.5, color=AROUTE, weight="bold", ha="right")

# hand-placed C1 entry chemistry: formate / methanol lane (left), FT lane (right)
lx_, rx_, mx_ = NX0 + 0.52, NX0 + 1.62, NX0 + 1.07
yy_ = [NY0 + 0.28 + i * 0.5 for i in range(11)]
CUSTOM = {"CO₂": (lx_, yy_[0]), "CO₂^δ−*": (lx_, yy_[1]), "HCOO*": (lx_, yy_[2]), "H₂COO*": (lx_, yy_[3]), "H₂CO*": (lx_, yy_[4]),
          "CH₃O*": (lx_, yy_[5]), "CH₃OH": (lx_, yy_[6]), "CH₃OH/SMS": (lx_ + 0.12, yy_[7] + 0.15), "DME": (mx_ + 0.42, yy_[8] + 0.05),
          "烃池 HCP": (lx_ + 0.12, yy_[9] + 0.15), "CO*": (rx_, yy_[1]), "HCO*": (rx_, yy_[2]), "HCOH*": (rx_, yy_[3]),
          "CH*/C*": (rx_, yy_[4]), "CH₂*": (rx_, yy_[5]), "CH₃*": (rx_, yy_[6]), "CH₄": (rx_, yy_[7]), "·CH₃": (mx_ + 0.02, yy_[7] - 0.18),
          "HCOOH": (mx_, yy_[3] - 0.25), "HCHO": (mx_, yy_[4] + 0.25)}
LAB_I = {}
for i, l in enumerate(labels):
    if l in CUSTOM and cls[i] != POOL and (cls[i] == CORE or xy[i, 0] <= C1W):
        LAB_I[i] = l


def Q(i):
    return CUSTOM[LAB_I[i]] if i in LAB_I else P(*xy[i])


tr = lambda pairs: [(Q(a), Q(q)) for a, q in pairs]  # noqa: E731
ax.add_collection(LineCollection(tr(segs[0]), colors=SHARED, linewidths=0.2, alpha=0.05, zorder=1.2, rasterized=True))
ax.add_collection(LineCollection(tr(segs[1]), colors=BONLY, linewidths=0.18, alpha=0.045, zorder=1.3, rasterized=True))
vis = (cls != POOL) & (cls != CORE) & (xy[:, 0] > C1W)
pts = np.array([P(*v) for v in xy[vis]])
ax.scatter(pts[:, 0], pts[:, 1], s=0.3, c="#4C5B70", linewidths=0, alpha=0.5, zorder=1.6, rasterized=True)
short = {"alc": "醇", "cbo": "醛 / 酮 / 酸", "ads": "吸附烷基 R*", "rad": "自由基 R•", "alk": "烷烃", "ole": "烯烃", "ia": "烷基碳正离子",
         "die": "二烯", "ie": "烯基碳正离子", "cyc": "环烷 / 环烯", "ic": "环状碳正离子", "cd": "环二烯", "aro": "芳烃"}
for r in D["rows"]:
    _, y0 = P(0, r["y0"])
    _, y1 = P(0, r["y1"])
    t = text(NX0 + NW - 0.06, (y0 + y1) / 2, short[r["key"]], size=8.2, color=INK2, ha="right", va="center", zorder=5)
    t.set_bbox(dict(boxstyle="round,pad=0.16,rounding_size=0.25", fc="white", ec="none", alpha=0.85))

# C1 panel: entry chemistry drawn explicitly
ax.add_collection(LineCollection(tr(core_segs[0]), colors=SHARED, linewidths=1.5, alpha=0.9, zorder=6))
ax.add_collection(LineCollection(tr(core_segs[1]), colors=BONLY, linewidths=1.5, alpha=0.9, zorder=6))
POS = {l: CUSTOM[l] for l in LAB_I.values()}
text(lx_, NY0 + NH - 0.12, "甲酸盐–甲醇", size=7.8, color=BONLY, ha="center")
text(rx_, NY0 + NH - 0.12, "费托 C1", size=7.8, color=SHARED, ha="center")
ENTRY = {"CO₂": BONLY, "CO*": AROUTE, "CH₃OH/SMS": BONLY}
B_SET = {"CO₂", "CO₂^δ−*", "HCOO*", "H₂COO*", "H₂CO*", "CH₃O*", "CH₃OH/SMS", "DME", "烃池 HCP"}
NICE = {"CO₂^δ−*": "CO₂δ⁻*", "CH₃OH/SMS": "CH₃OH / SMS", "烃池 HCP": "芳烃池"}
for l, (x, y) in POS.items():
    fill = ENTRY.get(l)
    col = BONLY if l in B_SET else SHARED
    t = text(x, y, NICE.get(l, l), size=7.6 if not fill else 8.4, color="white" if fill else col, ha="center", va="center",
             weight="bold" if fill else "normal", zorder=8)
    t.set_bbox(dict(boxstyle="round,pad=0.22,rounding_size=0.45", fc=fill or "white", ec=fill or col, lw=0.8))

# legend under the network
LY = NY0 + NH + 0.30
xx = NX0
for key, name in [("MET", "金属 / 碳化物"), ("GAS", "气相"), ("OX", "加氢"), ("BA", "分子筛酸位"), ("AR", "芳烃池")]:
    ax.add_patch(Rectangle((xx, LY - 0.12), 0.22, 0.16, fc=BAND[key], ec="#C9CED6", lw=0.5, zorder=2))
    t = text(xx + 0.3, LY + 0.01, name, size=8.6, color=INK2)
    xx += 0.3 + t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi + 0.3
for col, name in [(SHARED, "两条路线共用的基元反应 / 物种"), (BONLY, "仅 CO₂ 直接加氢使用（甲酸盐–甲醇入口、SMS 甲基化）")]:
    ax.plot([xx, xx + 0.4], [LY - 0.04] * 2, color=col, lw=2.4)
    t = text(xx + 0.48, LY + 0.01, name, size=8.6, color=INK2)
    xx += 0.48 + t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi + 0.35

# ------------------------------------------------------------------ left column: the two process chains
LX, LW = 0.45, 4.35


def step(y, h, n, color, title, sub, eq, chips):
    card(LX, y, LW, h)
    badge(LX + 0.3, y + 0.3, n, color)
    t = text(LX + 0.56, y + 0.37, title, size=11.8, weight="bold")
    w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
    text(LX + 0.56 + w + 0.12, y + 0.37, sub, size=8.6, color=MUTED)
    if eq:
        text(LX + 0.56, y + 0.68, eq, size=9, color=INK2)
    cx = LX + 0.56
    for s_ in chips:
        cx += chip(cx, y + h - 0.22, s_.lstrip("※"), BONLY if s_.startswith("※") else SHARED) + 0.07


text(LX, 1.60, "路线 A", size=11.5, color=AROUTE, weight="bold")
text(LX + 0.72, 1.60, "甲烷干重整 → 合成气费托 → 加氢提质", size=11, color=INK)
step(1.74, 0.80, 1, AROUTE, "甲烷干重整 DRM", "Ni / Co · 700–900 °C", "CH₄ + CO₂ ⇌ 2 CO + 2 H₂（另有积碳、RWGS）", [])
step(2.64, 1.06, 2, AROUTE, "费托合成 FTS", "Fe / Co · 200–340 °C", "CO* → HCO* → CHₓ* → 逐个碳增长（ASF 分布）",
     ["CH₂* 插入", "β-H 消除", "CO* 插入", "表面加氢"])
step(3.80, 1.06, 3, AROUTE, "加氢提质", "Pt / Pd–分子筛 · 加氢裂化 / 异构化 / 齐聚", "长链蜡、烯烃 → 支链 C8–C16",
     ["质子化", "PCP 支化", "β-断裂 / 齐聚", "加氢"])
text(LX, 5.22, "路线 B", size=11.5, color=BONLY, weight="bold")
text(LX + 0.72, 5.22, "CO₂ 直接加氢", size=11, color=INK)
card(LX, 5.36, LW, 1.06)
badge(LX + 0.3, 5.66, 1, BONLY)
text(LX + 0.56, 5.73, "C1 活化", size=11.8, weight="bold")
text(LX + 1.43, 5.73, "两种入口", size=8.6, color=MUTED)
text(LX + 0.56, 6.02, "① RWGS + 费托（Fe–K）：CO₂ → CO*，此后同路线 A", size=8.9, color=INK2)
text(LX + 0.56, 6.28, "② OX-ZEO（In₂O₃、ZnZrOₓ）：CO₂ → HCOO* → CH₃OH / SMS", size=8.9, color=INK2)
step(6.52, 1.06, 2, BONLY, "分子筛酸位", "HZSM-5 · 300–400 °C", "烯烃循环：SMS 甲基化 ⇄ β-断裂；氢转移、芳构化",
     ["※SMS 甲基化", "质子化", "β-断裂 / 齐聚", "氢转移"])


def arrow(p0, p1, color, rad=0.0, lw=1.6, z=9):
    ax.add_patch(FancyArrowPatch(p0, p1, connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>,head_length=5.5,head_width=3.2",
                                 color=color, lw=lw, zorder=z, shrinkA=1, shrinkB=10))


arrow((LX + LW, 2.12), POS["CO*"], AROUTE, rad=-0.15)
text(LX + LW + 0.06, 1.96, "CO + H₂", size=8.6, color=AROUTE, weight="bold")
arrow((LX + LW, 5.98), POS["CO₂"], BONLY, rad=-0.22, lw=1.3)
arrow((LX + LW, 6.24), POS["CH₃OH/SMS"], BONLY, rad=0.12)
arrow((LX + LW, 3.25), (NX0 + C1PX + 0.15, 2.95), AROUTE, rad=0.0, lw=1.1)
arrow((LX + LW, 4.33), (NX0 + C1PX + 0.15, 4.7), AROUTE, rad=0.0, lw=1.1)
arrow((LX + LW, 7.05), (NX0 + C1PX + 0.15, 5.25), BONLY, rad=0.0, lw=1.1)

# ------------------------------------------------------------------ bottom statistics (full network)
BY, BH = 7.92, 3.25
pct_sp = 100 * ST["species"]["both"] / ST["species"]["B"]
pct_rx = 100 * ST["reactions"]["both"] / ST["reactions"]["B"]
# (1) intermediates
card(LX, BY, LW, BH)
text(LX + 0.3, BY + 0.42, "可达中间体", size=12.5, weight="bold")
text(LX + 0.3, BY + 0.72, f"全量 C1–C16 网络，{ST['n_species']:,} 个物种", size=8.6, color=MUTED)
cx_, cy_ = LX + 1.2, BY + 1.72
ax.add_patch(Circle((cx_, cy_), 0.76, fc=BONLY + "16", ec=BONLY, lw=1.6, zorder=3))
ax.add_patch(Circle((cx_, cy_), 0.70, fc=SHARED + "22", ec=AROUTE, lw=1.6, ls=(0, (4, 2)), zorder=4))
text(cx_, cy_ - 0.08, "100%", size=19, weight="bold", color=SHARED, ha="center", va="center", zorder=6)
text(cx_, cy_ + 0.27, "C2 以上共享", size=8.6, color=INK2, ha="center", va="center", zorder=6)
X1 = LX + 2.35
text(X1, BY + 1.05, "路线 A 可达", size=9, color=AROUTE, weight="bold")
text(X1, BY + 1.33, f"{ST['species']['A']:,}", size=13.5, weight="bold")
text(X1, BY + 1.75, "路线 B 可达", size=9, color=BONLY, weight="bold")
text(X1, BY + 2.03, f"{ST['species']['B']:,}", size=13.5, weight="bold")
text(X1, BY + 2.38, f"共享 {pct_sp:.4f}%", size=9, color=SHARED, weight="bold")
text(LX + 0.3, BY + 2.66, f"只有 {ST['species']['B_only']} 个物种为路线 B 独有，全是 C1 入口：CO₂、CO₂δ⁻*、HCOO*、\nH₂COO*、H₂CO*、CH₃O*、CH₃OH / SMS、DME、芳烃池。", size=8.1, color=MUTED, linespacing=1.45, va="top")

# (2) reactions by family
X2 = LX + LW + 0.25
W2 = 7.05
card(X2, BY, W2, BH)
text(X2 + 0.3, BY + 0.42, "基元反应：按反应族", size=12.5, weight="bold")
text(X2 + 0.3, BY + 0.72, f"共享 {ST['reactions']['both']:,} / {ST['reactions']['B']:,} 条（{pct_rx:.1f}%）；路线 B 多出的几乎都是 SMS 甲基化",
     size=8.6, color=MUTED)
FAMN = {"Z:prot": "质子化", "Z:pcp": "PCP 支化", "Z:HT": "氢转移", "Z:hshift": "1,2-H 迁移", "Z:alkylshift": "烷基迁移",
        "Z:meth": "SMS 甲基化", "Z:beta": "β-断裂 / 齐聚", "H:hyd": "加氢 / 脱氢", "Z:cyclise": "成环", "M:betaH": "β-H 消除",
        "O:dehyd": "醇脱水", "M:hyd": "表面加氢", "O:OH": "R* + OH* 成醇"}
fams = sorted(ST["by_family"].items(), key=lambda kv: -(kv[1]["both"] + kv[1]["B_only"]))
fams = [kv for kv in fams if kv[0] in FAMN][:11]
mx = max(v["both"] + v["B_only"] for _, v in fams)
bx0, bw = X2 + 1.45, W2 - 2.35
for i, (f, v) in enumerate(fams):
    y = BY + 0.98 + i * 0.196
    tot = v["both"] + v["B_only"]
    text(bx0 - 0.1, y + 0.075, FAMN[f], size=8.3, color=BONLY if v["B_only"] else INK2, ha="right", va="center")
    ax.add_patch(Rectangle((bx0, y), bw * v["both"] / mx, 0.15, fc=SHARED, ec="none", zorder=3))
    ax.add_patch(Rectangle((bx0 + bw * v["both"] / mx, y), bw * v["B_only"] / mx, 0.15, fc=BONLY, ec="none", zorder=3))
    text(bx0 + bw * tot / mx + 0.06, y + 0.075, f"{tot / 1e6:.2f} M", size=7.8, color=MUTED, va="center")

# (3) what differs
X3 = X2 + W2 + 0.25
card(X3, BY, FW - X3 - 0.45, BH)
text(X3 + 0.3, BY + 0.42, "差别只在两处", size=12.5, weight="bold")
items = [
    ("C1 入口", AROUTE, "A：DRM 直接给出 CO + H₂，以 CO* 进入网络。B：CO₂ 经 RWGS 也得到 CO*；\n或者经甲酸盐 → 甲醇 → 表面甲氧基 SMS 进入分子筛。"),
    ("链增长的单体", BONLY, "A：金属 / 碳化物表面 CH₂* 插入，每次加一个碳。B（OX-ZEO）：SMS 给烯烃加甲基\n（烯烃循环）。两种方式长出的是同一批烷基、烯烃和碳正离子。"),
    ("其余完全相同", SHARED, f"异构化、裂化、齐聚、氢转移、成环芳构化共用同一批中间体。\nC8–C16 航煤区间可达异构体：A {ST['jet_C8_C16']['A']:,}，B {ST['jet_C8_C16']['B']:,}，完全重合。"),
]
yy = BY + 0.80
for h, col, body in items:
    ax.add_patch(Rectangle((X3 + 0.3, yy - 0.15), 0.05, 0.72, fc=col, ec="none", zorder=3))
    text(X3 + 0.45, yy, h, size=10.2, weight="bold", color=INK)
    text(X3 + 0.45, yy + 0.12, body, size=8.6, color=INK2, linespacing=1.5, va="top")
    yy += 0.80

foot = (f"数据：显式枚举的同分异构体反应网络。中间图为推荐版（C1–C12 全部异构体 + C13–C16 ZSM-5 择形，{RS['n_species']:,} 个物种、{RS['n_reactions']:,} 条基元反应），"
        f"下方统计为全量 C1–C16（{ST['n_species']:,} 个物种、{ST['n_reactions']:,} 条）。可达性：双向超图 BFS，H₂、H*、H⁺、H₂O 视为可用；"
        "路线 A 从 CO* 出发、不允许甲醇相关反应族，路线 B 从 CO₂ 出发、允许全部反应族。")
text(0.45, FH - 0.2, foot, size=7.8, color=MUTED)

out = NET / "figures" / "route_shared_network"
fig.savefig(str(out) + ".png", dpi=200, facecolor=PAPER)
if "--pdf" in sys.argv:
    fig.savefig(str(out) + ".pdf", facecolor=PAPER, dpi=300)
print("saved", out)
