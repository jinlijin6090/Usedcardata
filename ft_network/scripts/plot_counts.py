#!/usr/bin/env python3
"""随碳数增加, 表面中间体数和基元步骤数怎么变 (C1–C12)。
上排: 每个碳数新增的数量及其组成; 下排: 累计数量。"""
import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
NMAX = 12
SP = list(csv.DictReader(open(ROOT / f"data/c{NMAX}/species.csv", encoding="utf-8-sig")))
RX = list(csv.DictReader(open(ROOT / f"data/c{NMAX}/reactions.csv", encoding="utf-8-sig")))
plt.rcParams["font.sans-serif"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

INK, MUTE, GRID, SURF = "#1B2A33", "#6B7880", "#E3E7EA", "#FCFCFB"
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]     # 已用 validate_palette.js 校验
stage = lambda n: max(1, int(n))
NS = list(range(1, NMAX + 1))

SMALL = {"H*", "O*", "OH*", "H2O*", "CO*", "CO2*", "COOH*", "HCOO*"}
SP_TYPES = ["烃类中间体（CHx、CnHm）", "含氧中间体", "η3-烯丙基", "小分子吸附物种"]


def sp_type(r):
    if r["species"] in SMALL:
        return SP_TYPES[3]
    if r["type"] == "η3-烯丙基":
        return SP_TYPES[2]
    return SP_TYPES[1] if int(r["O"]) else SP_TYPES[0]


sp = Counter((sp_type(r), stage(r["carbon_number"])) for r in SP if r["phase"] == "surface")
RX_GROUPS = [("中间体加氢 / 脱氢", ["4 ", "5g"]), ("CO 插入中间体的加氢与 C–O 断裂", ["5c", "5d"]),
             ("C–C 成键（偶联、CO 插入、烯醇、烯基）", ["5a", "5b", "5e", "5f"]),
             ("CO 活化与除氧（仅 C1）", ["2a", "2b", "2c", "3a", "3b"]), ("吸附 / 脱附", ["1 ", "5h"])]
rx = Counter()
for r in RX:
    g = next(name for name, keys in RX_GROUPS if r["category"][:2] in keys)
    if r["equation"].startswith("CH4* ⇌"):
        g = "吸附 / 脱附"
    rx[(g, stage(r["n"]))] += 1

fig, axes = plt.subplots(2, 2, figsize=(15, 10.2), facecolor=SURF, gridspec_kw=dict(height_ratios=[1.25, 1], hspace=0.42, wspace=0.16))


def style(ax, ylab):
    ax.set_facecolor(SURF)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTE)
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTE, length=0, labelsize=11)
    ax.set_xticks(NS)
    ax.set_xticklabels([f"C{n}" for n in NS])
    ax.set_ylabel(ylab, color=MUTE, fontsize=11)


def stacked(ax, cats, data, title, ylab):
    style(ax, ylab)
    bottom = [0] * len(NS)
    for i, c in enumerate(cats):
        vals = [data[(c, n)] for n in NS]
        ax.bar(NS, vals, 0.72, bottom=bottom, color=PAL[i], edgecolor=SURF, linewidth=2, label=c)
        for n, v, b in zip(NS, vals, bottom):
            if v >= max(3, 0.045 * ax_max[ax]):
                ax.text(n, b + v / 2, str(v), ha="center", va="center", fontsize=9.5, color=INK)
        bottom = [b + v for b, v in zip(bottom, vals)]
    for n, t in zip(NS, bottom):
        ax.text(n, t + 0.012 * ax_max[ax], str(t), ha="center", va="bottom", fontsize=11, fontweight="bold", color=INK)
    ax.set_ylim(0, ax_max[ax])
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", color=INK, pad=54)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.01, 1.0), ncol=3 if len(cats) > 4 else 4, frameon=False, fontsize=9.5,
              handlelength=1.0, handleheight=1.0, columnspacing=1.2, labelcolor=INK)
    return bottom


def cumulative(ax, per, title, ylab, color, slope):
    style(ax, ylab)
    cum = [sum(per[:i + 1]) for i in range(len(per))]
    ax.plot(NS, cum, color=color, lw=2, marker="o", ms=6, mec=SURF, mew=1.5)
    for n, v in zip(NS, cum):
        if n in (1, 2, 3, 6, 9, 12):
            ax.text(n, v + 0.035 * cum[-1], str(v), ha="center", va="bottom", fontsize=11, color=INK,
                    fontweight="bold" if n == NMAX else "normal")
    ax.set_ylim(0, cum[-1] * 1.15)
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", color=INK, pad=12)
    ax.text(0.98, 0.08, f"C3 起是直线：每多一个碳 +{slope}", transform=ax.transAxes, ha="right", fontsize=11.5, color=MUTE)
    ax.set_xlabel("链增长到的碳数", color=MUTE, fontsize=11)
    return cum


ax_max = {axes[0, 0]: 42, axes[0, 1]: 95}
tot_sp = stacked(axes[0, 0], SP_TYPES, sp, "每个碳数新增的表面中间体", "新增中间体数")
tot_rx = stacked(axes[0, 1], [g for g, _ in RX_GROUPS], rx, "每个碳数新增的基元步骤", "新增基元步骤数")
cum_sp = cumulative(axes[1, 0], tot_sp, "累计表面中间体数", "累计中间体数", PAL[0], tot_sp[-1])
cum_rx = cumulative(axes[1, 1], tot_rx, "累计基元步骤数", "累计基元步骤数", PAL[0], tot_rx[-1])
assert sum(tot_rx) == len(RX) and sum(tot_sp) == sum(1 for r in SP if r["phase"] == "surface")

fig.suptitle("费托直链烃网络：C2 之后每个碳数新增的中间体和反应数固定不变", x=0.065, y=0.985, ha="left",
             fontsize=18, fontweight="bold", color=INK)
fig.text(0.065, 0.93, f"C1 = 通用小分子 + C1 物种；只计表面中间体，不含 {len(SP) - sum(tot_sp)} 个气相分子。"
         f"到 C{NMAX} 共 {cum_sp[-1]} 个表面中间体、{cum_rx[-1]} 个可逆基元步骤。", fontsize=11.5, color=MUTE)
fig.subplots_adjust(left=0.065, right=0.985, top=0.80, bottom=0.07)
out = ROOT / "figures" / f"ft_counts_vs_carbon_C1-C{NMAX}.png"
fig.savefig(out, dpi=200, facecolor=SURF)
print(out, tot_sp, cum_sp[-1], tot_rx, cum_rx[-1])
