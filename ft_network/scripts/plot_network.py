#!/usr/bin/env python3
"""绘制费托反应网络示意图: (a) C1 网络; (b) C_n 中间体格子与链增长/终止。"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

plt.rcParams["font.sans-serif"] = ["WenQuanYi Zen Hei", "Noto Sans CJK SC", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

C_H, C_CO, C_OH, C_CC, C_OUT = "#3b6fb6", "#c0392b", "#8e6bbf", "#2e8b57", "#d98200"


def node(ax, xy, text, fc="white", ec="#444", fs=9, bold=False):
    ax.text(*xy, text, ha="center", va="center", fontsize=fs, zorder=3,
            fontweight="bold" if bold else "normal",
            bbox=dict(boxstyle="round,pad=0.32", fc=fc, ec=ec, lw=0.9))


def arrow(ax, p, q, color, style="<|-|>", rad=0.0, lw=1.2, shrink=17, ls="-"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=9, color=color, lw=lw,
                                 shrinkA=shrink, shrinkB=shrink, zorder=1, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}"))


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 7.4), gridspec_kw=dict(width_ratios=[1, 1.12]))

# ------------------------------------------------------------------ (a) C1
ax = ax1
X = [0, 1.5, 3.0, 4.5]
top = {0: "CO*", 1: "HCO*", 2: "CH2O*", 3: "CH3O*"}
mid = {0: "COH*", 1: "HCOH*", 2: "CH2OH*"}
bot = {0: "C*", 1: "CH*", 2: "CH2*", 3: "CH3*"}
yT, yM, yB = 3.0, 1.9, 0.0
for p, t in top.items():
    node(ax, (X[p], yT), t, fc="#fdecea")
for p, t in mid.items():
    node(ax, (X[p], yM), t, fc="#f3ecfa")
for p, t in bot.items():
    node(ax, (X[p], yB), t, fc="#e8f1fb", bold=True)
for p in range(3):
    arrow(ax, (X[p], yT), (X[p + 1], yT), C_H)
    arrow(ax, (X[p], yB), (X[p + 1], yB), C_H)
    arrow(ax, (X[p], yT), (X[p], yM), C_OH)                      # O 端加氢
for p in range(2):
    arrow(ax, (X[p], yM), (X[p + 1], yM), C_H)
for p in range(4):                                               # C–O 断裂 (→ CHp + O)
    arrow(ax, (X[p], yT), (X[p], yB), C_CO, style="-|>", rad=-0.33)
for p in range(3):                                               # → CHp + OH
    arrow(ax, (X[p], yM), (X[p], yB), C_CO, style="-|>", ls="--")
node(ax, (-1.5, yT), "CO(g)", fc="#eee")
arrow(ax, (-1.5, yT), (X[0], yT), "#444")
node(ax, (-1.5, yB), "H2(g) ⇌ 2H*", fc="#eee")
node(ax, (6.0, yB), "CH4(g)", fc="#fff3d6", ec=C_OUT)
arrow(ax, (X[3], yB), (6.0, yB), C_OUT, style="-|>")
ax.text(5.3, yB + 0.22, "+H*", fontsize=8, color=C_OUT, ha="center")
# 除氧
node(ax, (0.6, -1.6), "O*", fc="#f5f5f5")
node(ax, (2.1, -1.6), "OH*", fc="#f5f5f5")
node(ax, (3.6, -1.6), "H2O*", fc="#f5f5f5")
node(ax, (5.3, -1.6), "H2O(g)", fc="#fff3d6", ec=C_OUT)
arrow(ax, (0.6, -1.6), (2.1, -1.6), C_H)
arrow(ax, (2.1, -1.6), (3.6, -1.6), C_H)
arrow(ax, (3.6, -1.6), (5.3, -1.6), C_OUT, style="-|>")
node(ax, (0.6, -2.7), "CO2*", fc="#f5f5f5")
node(ax, (2.1, -2.7), "COOH*", fc="#f5f5f5")
node(ax, (3.6, -2.7), "HCOO*", fc="#f5f5f5")
node(ax, (5.3, -2.7), "CO2(g)", fc="#fff3d6", ec=C_OUT)
arrow(ax, (0.6, -1.6), (0.6, -2.7), "#777")
arrow(ax, (2.1, -2.7), (0.6, -2.7), "#777")
arrow(ax, (3.6, -2.7), (0.6, -2.7), "#777", rad=0.35)
arrow(ax, (0.6, -2.7), (5.3, -2.7), C_OUT, style="-|>", rad=-0.42)
ax.text(-1.5, -1.6, "除氧 → H2O", fontsize=9, ha="center")
ax.text(-1.5, -2.7, "CO+O / CO+OH /\nHCO+O → CO2\n(水煤气变换)", fontsize=8.5, ha="center", va="center")
ax.text(-1.5, yM, "氢助路径\n含氧中间体", fontsize=9, ha="center", va="center", color=C_OH)
ax.text(2.25, 3.75, "(a) C1 网络：CO 活化 → CHx 单体 → CH4；除氧", fontsize=12, ha="center", fontweight="bold")
ax.text(2.25, 0.85, "C–O 断裂：实线 → CHp* + O*，虚线 → CHp* + OH*", fontsize=8.5, ha="center", color=C_CO,
        bbox=dict(fc="white", ec="none", pad=1.5), zorder=2)
ax.set_xlim(-2.5, 6.8)
ax.set_ylim(-3.4, 4.1)
ax.axis("off")

# ------------------------------------------------------------------ (b) Cn 格子
ax = ax2
XB = [0, 1.7, 3.4, 5.1]          # b = α-碳 H 数
YA = {0: 3.0, 1: 1.5, 2: 0.0}    # a = β-碳 H 数
lab = {
    (0, 0): "R–C≡C*\n炔基", (0, 1): "R–C≡CH*\n1-炔", (0, 2): "R–C(*)=CH2\n1-烯-2-基", (0, 3): "R–C(*)–CH3\n2-亚烷基",
    (1, 0): "R–CH=C*\n亚乙烯基型", (1, 1): "R–CH=CH*\n烯基", (1, 2): "R–CH=CH2*\nα-烯烃(吸附)", (1, 3): "R–CH(*)–CH3\n2-烷基",
    (2, 0): "R–CH2–C*\n次烷基", (2, 1): "R–CH2–CH*\n亚烷基", (2, 2): "R–CH2–CH2*\n烷基", (2, 3): "R–CH2–CH3*\n正构烷烃",
}
for (a, b), t in lab.items():
    fc = "#e3f4ea" if (a == 2 and b < 3) else ("#fff3d6" if (a, b) in ((2, 3), (1, 2)) else "white")
    node(ax, (XB[b], YA[a]), t, fc=fc, fs=8)
for a in range(3):
    for b in range(3):
        arrow(ax, (XB[b], YA[a]), (XB[b + 1], YA[a]), C_H, shrink=30)
for a in range(2):
    for b in range(4):
        arrow(ax, (XB[b], YA[a]), (XB[b], YA[a + 1]), C_H, shrink=19)
# 产物
node(ax, (6.9, YA[2]), "CnH2n+2(g)\n正构烷烃", fc="#fff3d6", ec=C_OUT, fs=8.5, bold=True)
arrow(ax, (XB[3], YA[2]), (6.9, YA[2]), C_OUT, shrink=33)
node(ax, (6.9, YA[1] + 0.75), "CnH2n(g)\n直链 α-烯烃", fc="#fff3d6", ec=C_OUT, fs=8.5, bold=True)
arrow(ax, (XB[2], YA[1]), (6.9, YA[1] + 0.75), C_OUT, rad=-0.2, shrink=30)
# 链增长入口/出口
ax.text(2.55, 4.15,
        "入口：C(n-1) 片段 R'–CH_a* (a=0,1,2) + CH_b* (b=0…3) → 格点 (a, b)\n"
        "或 + CO* → R'–CH_a–CH_pOH_q* → C–O 断裂 → 格点 (a, p)", fontsize=8.8, ha="center", color=C_CC,
        bbox=dict(boxstyle="round,pad=0.4", fc="#f2faf5", ec=C_CC, lw=0.9))
for b in range(4):
    arrow(ax, (XB[b], 3.85), (XB[b], YA[0]), C_CC, style="-|>", shrink=14, lw=1.0)
ax.text(2.55, -1.55,
        "出口（链增长到 C(n+1)）：绿色三格 次烷基/亚烷基/烷基 + CH_z*、+ CO*\n"
        "烯基 + CH2* → η3-烯丙基 → 高一级烯基（烯基机理）", fontsize=8.8, ha="center", color=C_CC,
        bbox=dict(boxstyle="round,pad=0.4", fc="#f2faf5", ec=C_CC, lw=0.9))
for b in range(3):
    arrow(ax, (XB[b], YA[2]), (XB[b], -1.2), C_CC, style="-|>", shrink=14, lw=1.0)
ax.text(-1.25, 1.5, "β-碳加氢 ↓\n(逆向: β-H 消除)", fontsize=8.5, ha="center", va="center", color=C_H, rotation=90)
ax.text(2.55, 3.45, "α-碳加氢 →", fontsize=8.5, ha="center", color=C_H,
        bbox=dict(fc="white", ec="none", pad=1), zorder=2)
ax.text(3.0, 4.95, "(b) C_n (n ≥ 3) 中间体格子：链增长、加氢/脱氢与链终止", fontsize=12, ha="center", fontweight="bold")
ax.set_xlim(-1.7, 8.0)
ax.set_ylim(-2.2, 5.3)
ax.axis("off")

fig.text(0.5, 0.015, "蓝: 加氢/脱氢 (±H*)   紫: O 端加氢   红: C–O 断裂   绿: C–C 成键(链增长)   橙: 脱附(产物)   "
         "所有步骤可逆；R = CH3(CH2)n-3–", ha="center", fontsize=9.5)
out = Path(__file__).resolve().parent.parent / "figures" / "ft_network_scheme.png"
out.parent.mkdir(exist_ok=True)
fig.tight_layout(rect=(0, 0.035, 1, 1))
fig.savefig(out, dpi=200)
print(out)
