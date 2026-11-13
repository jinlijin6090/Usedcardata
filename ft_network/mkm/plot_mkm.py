#!/usr/bin/env python3
"""MKM 结果图: (a) 基元反应网络与能垒; (b) 产物分布 (ASF 图); (c) 各指标与实验对比; (d) 参数调整量与速率控制度。"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from params_lit import DFT_STEP, TARGETS  # noqa: E402

BEST = sys.argv[1] if len(sys.argv) > 1 else "D"
VARIANTS = [v for v in ("A", "B", "C", "D") if (HERE / f"fit_result_{v}.json").exists()]
FIT = {v: json.loads((HERE / f"fit_result_{v}.json").read_text()) for v in VARIANTS}
AN = json.loads((HERE / "analysis.json").read_text())

plt.rcParams["font.sans-serif"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
INK, MUTE, GRID, SURF = "#1B2A33", "#6B7880", "#E3E7EA", "#FCFCFB"
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]       # 已用 validate_palette.js 校验
COL = {"DFT": PAL[1], "A": PAL[3], "B": PAL[2], "C": PAL[4], "D": PAL[0]}
LAB = {"DFT": "纯 DFT 参数", "A": "拟合 A：单位点", "B": "拟合 B：+ θCO 约束", "C": "拟合 C：+ CO 解离需两个空位", "D": "拟合 D：+ H 独立位点"}


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTE)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTE, labelsize=10.5)


def title(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=14, fontweight="bold", color=INK, pad=24)
    ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=10, color=MUTE, va="bottom")


fig = plt.figure(figsize=(17, 13.2), facecolor=SURF)
gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.17, left=0.12, right=0.985, top=0.885, bottom=0.05,
                      width_ratios=[1.15, 1])

# ------------------------------------------------------------------ (a) 网络
ax = fig.add_subplot(gs[0, 0])
ax.set_xlim(0, 12)
ax.set_ylim(8.4, -0.3)
ax.axis("off")
title(ax, "(a) MKM 用到的基元反应网络", f"边上数字：活化能 Ea（eV），DFT 文献值 → 拟合 {BEST} 后的值；粗体为速率或选择性控制步")
fit = FIT[BEST]["fit"]
dft = FIT[BEST]["dft"]
N = {"CO(g)": (0.9, 0.3), "CO*": (0.9, 1.7), "HCO*": (3.0, 1.7), "C*": (0.9, 3.5), "CH*": (3.0, 3.5), "CH2*": (5.1, 3.5),
     "CH3* = R1": (7.4, 3.5), "CH4(g)": (7.4, 1.7), "O*": (0.9, 5.6), "OH*": (3.0, 5.6), "H2O(g)": (5.1, 5.6),
     "CO2(g)": (3.0, 7.4), "H2(g) ⇌ 2H*": (5.1, 0.3),
     "Rn*（烷基）": (10.0, 3.5), "R(n+1)*": (10.0, 1.4), "CnH2n+2(g)": (11.0, 5.7), "En*（烯烃）": (8.9, 5.7),
     "CnH2n(g)": (8.9, 7.6)}
GAS = {"CO(g)", "CH4(g)", "H2O(g)", "CO2(g)", "CnH2n+2(g)", "CnH2n(g)", "H2(g) ⇌ 2H*"}
for k, (x, y) in N.items():
    w = 0.16 * len(k) + 0.35 if "（" in k or "⇌" in k or "=" in k else 0.13 * len(k) + 0.45
    fc, ec = ("#FFE3CC", PAL[1]) if k in GAS else ("#E4ECF7", PAL[0])
    if k.startswith("R") or k.startswith("E"):
        fc, ec = "#DFF3E6", "#2E8B57"
    ax.add_patch(FancyBboxPatch((x - w / 2, y - 0.27), w, 0.54, boxstyle="round,pad=0,rounding_size=0.1", fc=fc, ec=ec, lw=1.2, zorder=3))
    ax.text(x, y, k, ha="center", va="center", fontsize=11, color=INK, zorder=4)
KEY = {"dir", "ch4", "par", "grow"}


def edge(a, b, step, dx=0.0, dy=0.0, rad=0.0, label_pos=0.5, two=False):
    (x1, y1), (x2, y2) = N[a], N[b]
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), zorder=2,
                arrowprops=dict(arrowstyle="<|-|>" if two else "-|>", color=MUTE, lw=1.3, shrinkA=17, shrinkB=17,
                                connectionstyle=f"arc3,rad={rad}"))
    if step is None:
        return
    mx, my = x1 + (x2 - x1) * label_pos + dx, y1 + (y2 - y1) * label_pos + dy
    txt = f"{dft[step]:.2f} → {fit[step]:.2f}"
    ax.text(mx, my, txt, ha="center", va="center", fontsize=9.5, color=INK, zorder=5,
            fontweight="bold" if step in KEY else "normal", bbox=dict(fc=SURF, ec="none", pad=1.2))


edge("CO(g)", "CO*", None, two=True)
ax.text(1.15, 1.0, f"吸附 {dft['E_CO']:.2f} → {fit['E_CO']:.2f}", fontsize=9.5, color=INK, va="center")
edge("CO*", "C*", "dir", dx=-0.62)
edge("CO*", "HCO*", "hco", dy=-0.28, two=True)
edge("HCO*", "CH*", "hcod", dx=0.6)
edge("C*", "CH*", "ch", dy=-0.28, two=True)
edge("CH*", "CH2*", "ch2", dy=-0.28, two=True)
edge("CH2*", "CH3* = R1", "ch3", dy=-0.28, two=True)
edge("CH3* = R1", "CH4(g)", "ch4", dx=0.62)
edge("CO*", "O*", None, rad=0.35)
edge("O*", "OH*", "oh", dy=-0.28, two=True)
edge("OH*", "H2O(g)", "h2o", dy=-0.28)
ax.text(4.05, 6.15, f"2OH → H2O + O：{dft['ohoh']:.2f} → {fit['ohoh']:.2f}", fontsize=9.5, color=INK, ha="center")
edge("O*", "CO2(g)", "co2", dx=-0.25, dy=0.25, rad=-0.15)
edge("CH3* = R1", "Rn*（烷基）", None)
edge("Rn*（烷基）", "R(n+1)*", "grow", dx=0.95)
ax.text(9.85, 2.45, "+ CH2*", fontsize=10, color=MUTE, ha="right", va="center")
edge("Rn*（烷基）", "CnH2n+2(g)", "par", dx=0.62)
edge("Rn*（烷基）", "En*（烯烃）", "beta", dx=-0.68, two=True)
edge("En*（烯烃）", "CnH2n(g)", None)
ax.text(8.9 - 0.15, 6.65, f"脱附 {dft['E_des0']:.2f} → {fit['E_des0']:.2f}\n每个碳 +{fit['c_des']:.3f}；乙烯 +{fit['d_C2']:.2f}",
        fontsize=9.5, color=INK, ha="right", va="center")
ax.text(0.1, 8.25, "R_n 为烷基链（n = 1…80），链增长、终止、烯烃再加氢对每个 n 重复；H* 为共反应物，未画出。",
        fontsize=10, color=MUTE)

# ------------------------------------------------------------------ (b) 产物分布
ax = fig.add_subplot(gs[0, 1])
style(ax)
n = np.arange(1, 41)
alpha, sch4 = 0.90, 0.075
ref = (1 - alpha) * alpha ** (n - 1)
ax.plot(n, np.log10(ref), color=MUTE, lw=1.5, ls="--", label="ASF 参照线 α = 0.90（实验 α 0.88–0.92）")
ax.axhspan(math.log10(0.05), math.log10(0.10), xmin=0, xmax=0.06, color=GRID, zorder=0)
for v in ["DFT"] + VARIANTS:
    mol = np.array(AN[v]["mol"])[1:41]
    ok = mol > 0
    ax.plot(n[ok], np.log10(mol[ok]), color=COL[v], lw=2, marker="o", ms=3.8, mec=SURF, mew=0.8, label=LAB[v])
ax.annotate("实验 CH4 选择性 5–10 %", (1, math.log10(0.07)), xytext=(28, 18), textcoords="offset points", fontsize=10.5,
            color=INK, arrowprops=dict(arrowstyle="-", color=MUTE, lw=0.8))
ax.set_xlabel("碳数 n", color=MUTE, fontsize=11)
ax.set_ylabel("log10（摩尔分数）", color=MUTE, fontsize=11)
ax.set_ylim(-6, 0.15)
ax.set_xlim(0, 40.5)
ax.legend(loc="upper right", frameon=False, fontsize=10, labelcolor=INK)
title(ax, "(b) 产物分布：纯 DFT 几乎只出甲烷，拟合后接近 α≈0.9",
      "210 °C，20 bar，H2/CO = 2.1；直线越平 α 越大；拟合曲线约 C10 后变平（双 α）")

# ------------------------------------------------------------------ (c) 指标对比
ax = fig.add_subplot(gs[1, 0])
style(ax)
names = {"log10_TOF": "log10 TOF", "S_CH4": "CH4 选择性", "alpha": "α", "Ea": "表观活化能", "n_H2": "H2 级数",
         "n_CO": "CO 级数", "S_C5": "C5+ 选择性", "log10_OP3": "log10 C3 烯/烷"}
keys = list(names)
y = np.arange(len(keys))
ax.axvspan(-1, 1, color="#E8F3EC", zorder=0)
ax.axvline(0, color="#2E8B57", lw=1)
off = {"DFT": -0.3, "A": -0.15, "B": 0.0, "C": 0.15, "D": 0.3}
src = FIT[VARIANTS[0]]
tg = {k: t for k, t, *_ in TARGETS}
tol = {k: s for k, _, s, *_ in TARGETS}
for v in ["DFT"] + VARIANTS:
    pred = src["pred_dft"] if v == "DFT" else FIT[v]["pred_fit"]
    z = [np.clip((pred[k] - tg[k]) / tol[k], -8, 8) for k in keys]
    ax.scatter(z, y + off[v], s=55, color=COL[v], edgecolor=SURF, linewidth=1, zorder=3, label=LAB[v])
ax.set_yticks(y)
ax.set_yticklabels([f"{names[k]}（实验 {tg[k]:g}）" for k in keys], fontsize=10.5, color=INK)
ax.invert_yaxis()
ax.set_xlim(-8.5, 8.5)
ax.set_xlabel("（预测 − 实验）/ 容差；绿色带 = 在容差内，超出 ±8 截断显示", color=MUTE, fontsize=10.5)
ax.legend(loc="lower right", frameon=False, fontsize=9.5, labelcolor=INK, ncol=1)
title(ax, "(c) 各指标与实验的偏差", "纯 DFT 参数多数指标偏离很远；拟合的好坏以落在绿色带内的数目衡量")

# ------------------------------------------------------------------ (d) 参数调整 + DRC
ax = fig.add_subplot(gs[1, 1])
style(ax)
pn = [k for k, *_ in DFT_STEP]
sig = FIT[BEST]["sigma"]
shift = [(fit[k] - dft[k]) / sig[k] for k in pn]
yy = np.arange(len(pn))
ax.barh(yy, shift, color=[PAL[0] if abs(s) <= 1 else PAL[1] for s in shift], height=0.62)
ax.axvline(0, color=MUTE, lw=1)
for s in (-1, 1):
    ax.axvline(s, color=MUTE, lw=0.8, ls=":")
drc = AN.get(BEST, {}).get("drc", {})
for i, k in enumerate(pn):
    d = drc.get(k)
    if d and (abs(d["TOF"]) > 0.05 or abs(d["S_CH4"]) > 0.05):
        ax.text(2.35, i, f"DRC：TOF {d['TOF']:+.2f}，CH4 {d['S_CH4']:+.2f}", va="center", fontsize=9.5, color=INK)
ax.set_yticks(yy)
ax.set_yticklabels(pn, fontsize=10, color=INK, fontfamily="DejaVu Sans Mono")
ax.invert_yaxis()
ax.set_xlim(-2.5, 4.6)
ax.set_xlabel("拟合值 − DFT 值，以先验不确定度 σ 为单位（蓝 ≤ 1σ，橙 > 1σ）", color=MUTE, fontsize=10.5)
title(ax, f"(d) 拟合 {BEST} 需要把参数调多少，哪些步骤控制结果",
      "DRC = 速率控制度 ∂ln y / ∂ln k；只标出 |DRC| > 0.05 的步骤")

fig.suptitle("钴基费托 MKM：文献 DFT 参数能否复现实验？", x=0.05, y=0.975, ha="left", fontsize=19, fontweight="bold", color=INK)
fig.text(0.05, 0.935, "平均场集总模型（台阶位），21 个可调参数以 DFT 文献值为先验；实验基准为 210 °C、20 bar、H2/CO = 2.1、"
         "低转化率、Co 颗粒 > 10 nm。文献数值来自摘要，未核对原文。", fontsize=11.5, color=MUTE)
out = HERE.parent / "figures" / "co_mkm_fit.png"
fig.savefig(out, dpi=170, facecolor=SURF)
print(out)
