#!/usr/bin/env python3
"""把 C1–C8 的完整网络画成一张节点-箭头图: 258 个物种 = 258 个节点, 578 个基元步骤 = 578 条箭头。

布局: 每个碳数一个区块, 区块沿阶梯从左上排到右下。区块内每个链增长片段(次烷基/亚烷基/烷基)
对应一条"带", 带内中间一行(列)是烃类中间体, 两侧是 CO 插入得到的含氧中间体。
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import build_network as bn  # noqa: E402

NMAX = 8
bn.build(NMAX)
RX = list(csv.DictReader(open(ROOT / "data/c8/reactions.csv", encoding="utf-8-sig")))
SPEC = {s.name: s for s in bn.SPECIES.values()}

plt.rcParams["font.sans-serif"] = ["Liberation Sans", "WenQuanYi Zen Hei", "DejaVu Sans"]
W, H = 56.0, 31.5
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.axis("off")

# 颜色
K_H, K_OH, K_CO, K_CC, K_INS = "#3B6FB6", "#7A6FD0", "#C0392B", "#2E8B57", "#8E44AD"
K_ENOL, K_ALK, K_DES, K_WGS = "#008B8B", "#8B5A2B", "#D9580D", "#5F6B73"
NH = 0.25          # 节点高
FS = 6.0

# ------------------------------------------------------------------ 几何
PU, SV, BH = 1.30, 0.50, 2.00       # H 型: 列距, 子行距, 带高
BW, UV = 3.45 * PU, 0.80            # V 型: 带宽, 行距
LANE = 2.62
C2X0, C2Y0 = 6.9, 3.6
STEPX = 2 * BW + LANE * PU + 1.45
ORI = {}                            # n -> (type, x0, y0)
ORI[2] = ("H", C2X0, C2Y0)
x, y = C2X0, C2Y0
for n in range(3, NMAX + 1):
    if n % 2:                        # V 型放在上一个 H 型正下方
        prev_b = 4 if n == 3 else 3
        y = y + (prev_b - 1) * BH + LANE * SV + 1.15
        ORI[n] = ("V", x, y)
    else:                            # H 型放在上一个 V 型右侧
        x = x + STEPX
        ORI[n] = ("H", x, y)


def P(n, u, band, v):
    t, x0, y0 = ORI[n]
    if t == "H":
        return (x0 + u * PU, y0 + band * BH + v * SV)
    return (x0 + band * BW + v * PU, y0 + u * UV)


POS, STYLE = {}, {}


def put(name, xy, style):
    assert name in SPEC and name not in POS, name
    POS[name], STYLE[name] = xy, style


def width(name):
    return min(0.98, max(0.42, 0.052 * len(name) + 0.12))


# C1 区
XA, XB, XC, XD = 1.9, 3.15, 4.4, 5.6
row = lambda z: C2Y0 + z * BH + SV            # 与 C2 各带的烃行对齐
for p, nm in enumerate(["CO*", "HCO*", "CH2O*", "CH3O*"]):
    put(nm, (XB, row(p)), "oxy" if p else "small")
for p, nm in enumerate(["COH*", "HCOH*", "CH2OH*"]):
    put(nm, (XC, row(p)), "oxy")
for z in range(4):
    put(bn.c1(z).name, (XD, row(z)), "frag")
put("CO(g)", (XB, row(0) - 1.0), "gasin")
put("COOH*", (XA, row(0)), "small")
put("CO2*", (XA, row(0) + 1.0), "small")
put("HCOO*", (XA, row(1)), "small")
put("CO2(g)", (0.85, row(0) + 1.0), "gas")
put("CH4*", (XD, row(3) + 0.85), "hc")
put("CH4(g)", (XD, row(3) + 1.7), "gas")
YS = row(3) + 3.0
put("H2(g)", (XA, YS), "gasin")
put("H*", (XB, YS), "small")
put("O*", (XA, YS + 1.0), "small")
put("OH*", (XB, YS + 1.0), "small")
put("H2O*", (XC, YS + 1.0), "small")
put("H2O(g)", (XD, YS + 1.0), "gas")

# C2 … C8
FRAG = {1: [bn.c1(z) for z in range(4)]}
for n in range(2, NMAX + 1):
    FRAG[n] = [f for f, _ in bn.frags(n)]
for n in range(2, NMAX + 1):
    if n == 2:
        for a in range(4):
            for b in range(a + 1):
                st = "frag" if (a == 3 and b < 3) else "hc"
                put(bn.hc(2, a, b).name, P(2, b, a, 1), st)
        put(bn.alkane_g(2).name, P(2, 3.8, 3, 1), "gas")
        put(bn.alkene_g(2).name, P(2, 2.8, 2, 1), "gas")
    else:
        for a in range(3):
            for b in range(4):
                st = "frag" if (a == 2 and b < 3) else "hc"
                put(bn.hc(n, a, b).name, P(n, b, a, 1), st)
        put(bn.allyl(n).name, P(n, 3.55, 1, 2), "allyl")
        put(bn.alkane_g(n).name, P(n, 3.8, 2, 1), "gas")
        put(bn.alkene_g(n).name, P(n, 3.62, 1, 0), "gas")
    for a, f in enumerate(FRAG[n - 1]):
        for p in range(3):
            for q in range(2):
                put(bn.oxy(f, p, q).name, P(n, p + 0.5, a, 0 if q == 0 else 2), "oxy")

assert set(POS) == set(SPEC), set(SPEC) - set(POS)

# ------------------------------------------------------------------ 画节点
FILL = {"hc": ("#FFFFFF", "#34495E"), "frag": ("#DFF3E6", "#2E8B57"), "oxy": ("#F7ECF7", "#9B6BA6"),
        "small": ("#EEF3F4", "#5F6B73"), "gas": ("#FFE3CC", "#D9580D"), "gasin": ("#E4ECF7", "#3B6FB6"),
        "allyl": ("#F3E9DC", "#8B5A2B")}
for nm, (x, y) in POS.items():
    fc, ec = FILL[STYLE[nm]]
    if SPEC[nm].sites and STYLE[nm] == "hc" and nm != "CH4*" and 2 * SPEC[nm].C + 2 == SPEC[nm].H:
        fc = "#FFF3E6"                         # 吸附的烷烃
    w = width(nm)
    ax.add_patch(FancyBboxPatch((x - w / 2, y - NH / 2), w, NH, boxstyle="round,pad=0,rounding_size=0.05",
                                fc=fc, ec=ec, lw=0.7, zorder=3))
    ax.text(x, y, nm, ha="center", va="center", fontsize=FS, zorder=4,
            fontweight="bold" if STYLE[nm] in ("gas", "frag") else "normal", color="#1B2A33")


# ------------------------------------------------------------------ 画边
DRAWN = []


def edge_pt(nm, toward):
    """节点边框上朝向 toward 的点。"""
    x, y = POS[nm]
    w, h = width(nm) / 2 + 0.02, NH / 2 + 0.02
    dx, dy = toward[0] - x, toward[1] - y
    if dx == 0 and dy == 0:
        return (x, y)
    s = min(w / abs(dx) if dx else 1e9, h / abs(dy) if dy else 1e9)
    return (x + dx * s, y + dy * s)


def arrow(p, q, color, ls="-", rad=0.0, lw=0.75, z=2):
    ax.annotate("", xy=q, xytext=p, zorder=z,
                arrowprops=dict(arrowstyle="-|>,head_length=0.35,head_width=0.16", color=color, lw=lw, linestyle=ls,
                                shrinkA=0, shrinkB=0, connectionstyle=f"arc3,rad={rad}", mutation_scale=9))


def rid_label(xy, rid, color):
    ax.text(xy[0], xy[1], rid[1:].lstrip("0"), fontsize=3.4, color=color, ha="center", va="center", zorder=5,
            bbox=dict(fc="white", ec="none", pad=0.25, alpha=0.9))


def straight(r, s, d, color, ls="-", rad=0.0, lab=0.5):
    if rad:
        p, q = edge_pt(s, POS[d]), edge_pt(d, POS[s])
    else:
        p, q = edge_pt(s, POS[d]), edge_pt(d, POS[s])
    arrow(p, q, color, ls, rad)
    mx, my = p[0] + (q[0] - p[0]) * lab, p[1] + (q[1] - p[1]) * lab
    if rad:                                    # 弧线中点偏移
        mx -= rad * (q[1] - p[1]) * 0.5
        my += rad * (q[0] - p[0]) * 0.5
    rid_label((mx, my), r["id"], color)
    DRAWN.append(r["id"])


def poly(pts, color, lw=0.75, head=True):
    xs, ys = zip(*pts)
    ax.plot(xs[:-1] if head else xs, ys[:-1] if head else ys, color=color, lw=lw, zorder=1.5, solid_joinstyle="round")
    if head:
        arrow(pts[-2], pts[-1], color, lw=lw, z=1.5)


def parse(side):
    out = []
    for t in side.split(" + "):
        t = t.strip()
        k = 1
        if " " in t:
            k, t = t.split(" ")
            k = int(k)
        if t != "*":
            out += [t] * k
    return out


SMALL = {"H*", "O*", "OH*"}


def endpoints(r):
    L, R = (parse(s) for s in r["equation"].split(" ⇌ "))
    cat = r["category"][:2]
    eq = r["equation"]
    if cat == "1 ":
        return L[0], R[0]
    if cat == "3a":
        return {"O* + H*": ("O*", "OH*"), "OH* + H*": ("OH*", "H2O*"), "2 OH*": ("OH*", "H2O*"),
                "H2O*": ("H2O*", "H2O(g)")}[eq.split(" ⇌ ")[0]]
    if eq.startswith("2 CO*"):
        return "CO*", "C*"
    cl = [s for s in L if s not in SMALL]
    cr = [s for s in R if s not in SMALL and s != "H2O*"]
    if len(set(cl)) > 1 and cat not in ("3b", "2a"):
        cl = [s for s in cl if s != "CO*"] or cl
    if len(set(cl)) > 1 and cat == "5e":
        cl = [s for s in cl if s != "HCOH*"]
    if cat == "3b" and len(set(cl)) > 1:
        cl = [s for s in cl if s not in ("H2O*",)]
    s = max(cl, key=lambda nm: (SPEC[nm].C, SPEC[nm].H, nm != "CO*"))
    d = max(cr, key=lambda nm: (SPEC[nm].C, SPEC[nm].H))
    return s, d


def block_of(nm):
    return SPEC[nm].C


# --- 5a / 5b: 跨区块的链增长, 走"干线 + 分支"
by_prod = {}
for r in RX:
    if r["category"][:2] in ("5a", "5b"):
        s, d = endpoints(r)
        by_prod[(r["category"][:2], s, d)] = r

for n in range(2, NMAX + 1):
    t, x0, y0 = ORI[n]
    srcs = FRAG[n - 1]
    lines = []                                   # (目标坐标, kind, band, frag)
    for a, f in enumerate(srcs):
        lines.append(("P", a, f))
        lines.append(("G", a, f))
    for i, (kind, a, f) in enumerate(lines):
        sx, sy = POS[f.name]
        sw = width(f.name) / 2
        off = 0.045 if kind == "P" else 0.115
        if t == "H":                             # 从左侧进入
            start = (sx + sw, sy + (off - 0.08 if n == 2 else off))
            chx = x0 - 0.55 - 0.15 - 0.1 * i if n > 2 else x0 - 0.55 - 0.2 - 0.12 * (kind == "G")
            if kind == "P":
                tgt = bn.oxy(f, 0, 0).name
                tx, ty = POS[tgt]
                pts = [start, (chx, start[1]), (chx, ty), (tx - width(tgt) / 2, ty)]
                r = by_prod[("5b", f.name, tgt)]
                poly(pts, K_INS)
                rid_label((tx - width(tgt) / 2 - 0.22, ty), r["id"], K_INS)
                DRAWN.append(r["id"])
            else:
                ly = y0 + a * BH + LANE * SV
                targets = [bn.hc(n, a, z) for z in (range(a + 1) if n == 2 else range(4))]
                xend = max(POS[h.name][0] for h in targets) - 0.1
                poly([start, (chx, start[1]), (chx, ly), (xend, ly)], K_CC, head=False)
                for h in targets:
                    hx, hy = POS[h.name]
                    arrow((hx - 0.1, ly), (hx - 0.1, hy + NH / 2), K_CC, z=1.5)
                    r = by_prod[("5a", f.name, h.name)]
                    rid_label((hx - 0.1, ly - 0.2), r["id"], K_CC)
                    DRAWN.append(r["id"])
        else:                                    # 从上方进入
            start = (sx + off, sy + NH / 2)
            chy = y0 - 0.35 - 0.1 * i
            if kind == "P":
                tgt = bn.oxy(f, 0, 0).name
                tx, ty = POS[tgt]
                pts = [start, (start[0], chy), (tx, chy), (tx, ty - NH / 2)]
                r = by_prod[("5b", f.name, tgt)]
                poly(pts, K_INS)
                rid_label((tx, ty - NH / 2 - 0.13), r["id"], K_INS)
                DRAWN.append(r["id"])
            else:
                lx = x0 + a * BW + LANE * PU
                targets = [bn.hc(n, a, z) for z in range(4)]
                yend = max(POS[h.name][1] for h in targets) - 0.075
                poly([start, (start[0], chy), (lx, chy), (lx, yend)], K_CC, head=False)
                for h in targets:
                    hx, hy = POS[h.name]
                    arrow((lx, hy - 0.075), (hx + width(h.name) / 2, hy - 0.075), K_CC, z=1.5)
                    r = by_prod[("5a", f.name, h.name)]
                    rid_label((lx - 0.2, hy - 0.075), r["id"], K_CC)
                    DRAWN.append(r["id"])

# --- 其余反应: 每条一根箭头
pair_count = {}
for r in RX:
    cat = r["category"][:2]
    if cat in ("5a", "5b"):
        continue
    s, d = endpoints(r)
    k = pair_count.get((s, d), 0)
    pair_count[(s, d)] = k + 1
    (sx, sy), (dx, dy) = POS[s], POS[d]
    eq = r["equation"]
    if cat in ("1 ", "5h") or d.endswith("(g)"):
        straight(r, s, d, K_DES if d.endswith("(g)") else K_H)
    elif cat == "3a":
        straight(r, s, d, K_WGS, rad=0.5 * k)
    elif cat == "3b":
        straight(r, s, d, K_WGS, rad=[0, 0.45, -0.45][k])
    elif cat == "2a":
        straight(r, s, d, K_CO, rad=-0.22 - 0.16 * k, ls="-" if k == 0 else ":")
    elif cat in ("2c", "5d"):
        ls = "--" if "OH*" in parse(eq.split(" ⇌ ")[1]) else "-"
        if cat == "2c" and s in ("HCO*", "CH2O*"):
            straight(r, s, d, K_CO, ls, rad=-0.22)
        elif block_of(d) == 2 and abs(dy - sy) > 1.2:      # C2: 落到别的带
            straight(r, s, d, K_CO, ls, rad=0.12, lab=0.35)
        else:
            straight(r, s, d, K_CO, ls)
    elif cat == "5e":
        straight(r, s, d, K_ENOL, "--", rad=0.18, lab=0.5)
    elif cat == "5f":
        cross = block_of(s) != block_of(d)
        straight(r, s, d, K_ALK, rad=0.15 if cross else 0.0, lab=0.5 if cross else 0.45)
    else:                                                  # 加氢: 2b, 4, 5c, 5g
        ohyd = SPEC[d].O and ("OH" in d and "OH" not in s)
        col = K_OH if ohyd else K_H
        same_blk = block_of(s) == block_of(d) and cat == "5g" and block_of(s) >= 2
        t = ORI[block_of(s)][0] if same_blk else None
        if same_blk and t == "H" and abs(dy - sy) > 1.0:           # β-碳加氢 (跨带, 竖直)
            p, q = (sx + 0.09, sy + NH / 2), (sx + 0.09, dy - NH / 2)
            arrow(p, q, col)
            rid_label((p[0], sy + 0.75 * SV), r["id"], col)
            DRAWN.append(r["id"])
        elif same_blk and t == "V" and abs(dx - sx) > 2.0:         # β-碳加氢 (跨带, 水平)
            p, q = (sx + width(s) / 2, sy + 0.075), (dx - width(d) / 2, sy + 0.075)
            arrow(p, q, col)
            rid_label((sx + width(s) / 2 + 0.18, sy + 0.075), r["id"], col)
            DRAWN.append(r["id"])
        else:
            straight(r, s, d, col)

assert sorted(DRAWN) == sorted(r["id"] for r in RX), (len(DRAWN), len(set(DRAWN)))

# ------------------------------------------------------------------ 区块标签
plt.rcParams["font.sans-serif"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]   # 以下为中文文字
ax.text(XB + 0.6, C2Y0 - 1.75, "C1", fontsize=40, fontweight="bold", color="#0F4C5C", alpha=0.9, ha="center")
ax.text(XB + 0.6, C2Y0 - 1.35, "CO 活化 → CHx 单体", fontsize=11, color="#0F4C5C", ha="center")
ax.text(XB + 0.6, YS - 0.55, "氢解离与除氧", fontsize=11, color="#0F4C5C", ha="center")
for n in range(2, NMAX + 1):
    t, x0, y0 = ORI[n]
    if t == "H":
        ax.text(x0 + 1.5 * PU, y0 - 0.95, f"C{n}", fontsize=40, fontweight="bold", color="#0F4C5C", ha="center", va="center")
        for a, f in enumerate(FRAG[n - 1]):
            ax.text(x0 + 3.15 * PU, y0 + a * BH - 0.36, f"来自 {f.name}", fontsize=7, color="#8A9AA3", ha="right", va="center")
    else:
        ax.text(x0 - 1.05, y0 + 1.5 * UV, f"C{n}", fontsize=40, fontweight="bold", color="#0F4C5C", ha="center", va="center")
        for a, f in enumerate(FRAG[n - 1]):
            ax.text(x0 + a * BW + PU, y0 + 3.75 * UV, f"来自 {f.name}", fontsize=7, color="#8A9AA3", ha="center", va="center")

# ------------------------------------------------------------------ 图例
LX, LY = 22.5, 3.3
T = "#1B2A33"
ax.text(LX, LY, "节点", fontsize=24, fontweight="bold", color="#0F4C5C", va="center")
items = [("hc", "烃类表面中间体（自由基 / 卡宾 / 卡拜型）"), ("frag", "链增长片段：再接 CHx 或 CO"), ("oxy", "含氧中间体（CO 插入 / 氢助路径）"),
         ("allyl", "η3-烯丙基"), ("small", "小分子吸附物种"), ("gas", "气相产物"), ("gasin", "气相反应物")]
for i, (st, txt) in enumerate(items):
    yy = LY + 0.9 + i * 0.72
    fc, ec = FILL[st]
    ax.add_patch(FancyBboxPatch((LX, yy - 0.2), 1.1, 0.4, boxstyle="round,pad=0,rounding_size=0.07", fc=fc, ec=ec, lw=1.4))
    ax.text(LX + 1.4, yy, txt, fontsize=18, va="center", color=T)
EX = LX + 11.5
ax.text(EX, LY, "箭头（每条 = 一个可逆基元步骤）", fontsize=24, fontweight="bold", color="#0F4C5C", va="center")
edges = [(K_H, "-", "加氢 / 脱氢（± H*）"), (K_OH, "-", "O 端加氢"), (K_CO, "-", "C–O 断裂，生成 O*"), (K_CO, "--", "C–O 断裂，生成 OH*"),
         (K_CC, "-", "C–C 偶联：片段 + CHx*（碳化物机理）"), (K_INS, "-", "CO 插入：片段 + CO*"),
         (K_ENOL, "--", "烯醇缩合：+ HCOH*，− H2O*"), (K_ALK, "-", "烯基机理：烯基 + CH2* → 烯丙基"),
         (K_DES, "-", "脱附 → 气相产物"), (K_WGS, "-", "除氧 / 水煤气变换")]
for i, (c, ls, txt) in enumerate(edges):
    col, rw = divmod(i, 5)
    yy = LY + 0.9 + rw * 0.72
    xx = EX + col * 10.3
    arrow((xx, yy), (xx + 1.5, yy), c, ls, lw=2.2)
    ax.text(xx + 1.8, yy, txt, fontsize=18, va="center", color=T)
NX, NY = 1.0, 20.2
ax.text(NX, NY, "怎么读这张图", fontsize=24, fontweight="bold", color="#0F4C5C", va="center")
notes = ["每个节点是一个物种，共 258 个；每条箭头是一个可逆基元步骤，共 578 条，",
         "箭头指向书写的正方向，旁边的小数字是反应编号（1–578）。",
         "箭头只连接含碳主体；共反应物（H*、O*、OH*、CHx*、CO*、空位 *）",
         "由箭头颜色表示。",
         "绿色干线：上一碳数的片段 R–CHa* 进入对应的带，分支到 R–CHa–CHb*",
         "（b = 0…3），即分别与 C*、CH*、CH2*、CH3* 偶联。",
         "紫色线：同一片段 + CO* 得到 R–CHa–CO*，再经加氢和 C–O 断裂",
         "落回该带的烃类中间体。",
         "带内相邻的烃类节点差一个 α-氢，跨带相邻的差一个 β-氢。",
         "终点：R–CH2–CH3* 脱附为正构烷烃，R–CH=CH2* 脱附为直链 α-烯烃。"]
for i, tline in enumerate(notes):
    ax.text(NX, NY + 1.0 + i * 0.68, tline, fontsize=18, va="center", color=T)

out = ROOT / "figures"
fig.savefig(out / "ft_network_C1-C8_full.png", dpi=200, facecolor="white")
fig.savefig(out / "ft_network_C1-C8_full.pdf", facecolor="white")
print(f"节点 {len(POS)}, 箭头 {len(DRAWN)}")
