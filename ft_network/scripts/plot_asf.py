#!/usr/bin/env python3
"""ASF (Anderson–Schulz–Flory) 分布示意图。
(a) 不同 α 下的产物质量分布; (b) ASF 图: 理想直线与实际常见偏差(示意); (c) 各馏分质量分数随 α 的变化。"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
plt.rcParams["font.sans-serif"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
INK, MUTE, GRID, SURF = "#1B2A33", "#6B7880", "#E3E7EA", "#FCFCFB"
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]     # 已用 validate_palette.js 校验


def W(n, a):
    """链长 n 的质量分数。"""
    return n * (1 - a) ** 2 * a ** (n - 1)


def style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTE)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTE, labelsize=11)


def title(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=14, fontweight="bold", color=INK, pad=26)
    ax.text(0, 1.025, sub, transform=ax.transAxes, fontsize=10.5, color=MUTE, va="bottom")


fig = plt.figure(figsize=(16.5, 10.5), facecolor=SURF)
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.42, wspace=0.2, left=0.06, right=0.985, top=0.83, bottom=0.07)
ax1, ax2, ax3 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

# ------------------------------------------------------------------ (a) 质量分布
style(ax1)
n = np.arange(1, 41)
for a, c in zip((0.70, 0.85, 0.92), PAL):
    w = W(n, a) * 100
    ax1.plot(n, w, color=c, lw=2, marker="o", ms=4.5, mec=SURF, mew=1)
    k = {0.70: 2, 0.85: 5, 0.92: 24}[a]
    ax1.annotate(f"α = {a:.2f}", (n[k], w[k]), xytext=(14, 10), textcoords="offset points", fontsize=11.5,
                 color=INK, fontweight="bold", arrowprops=dict(arrowstyle="-", color=MUTE, lw=0.8))
ax1.set_xlabel("碳数 n", color=MUTE, fontsize=11.5)
ax1.set_ylabel("质量分数 Wn（%）", color=MUTE, fontsize=11.5)
ax1.set_xlim(0, 40.5)
ax1.set_ylim(0, 15)
title(ax1, "(a) α 越大，产物越重、分布越宽", "Wn = n (1−α)² α^(n−1)；峰值位置约在 n ≈ −1/ln α")

# ------------------------------------------------------------------ (b) ASF 图
style(ax2)
n = np.arange(1, 31)
a1, a2, nb = 0.72, 0.88, 10
ideal = np.log(W(n, a1) / n)
real = ideal.copy()
real[n > nb] = ideal[n == nb] + (n[n > nb] - nb) * np.log(a2)          # 长链段斜率变缓
real[0] += 0.9                                                       # CH4 偏多
real[1] -= 0.75                                                      # C2 偏少
ax2.plot(n, ideal, color=PAL[0], lw=2, label="理想 ASF：一条直线")
ax2.plot(n, real, color=PAL[1], lw=0, marker="o", ms=6.5, mec=SURF, mew=1.2, label="常见的实际数据（示意）")
ax2.annotate("CH4 偏多", (1, real[0]), xytext=(28, 4), textcoords="offset points", fontsize=11.5, color=INK,
             arrowprops=dict(arrowstyle="-", color=MUTE, lw=0.8))
ax2.annotate("C2 偏少\n（乙烯再吸附）", (2, real[1]), xytext=(18, -40), textcoords="offset points", fontsize=11.5, color=INK,
             arrowprops=dict(arrowstyle="-", color=MUTE, lw=0.8))
ax2.annotate(f"约 C{nb} 后斜率变缓：\n长链 α 更大（烯烃再吸附 / 双活性位）", (20, real[n == 20][0]), xytext=(-10, 34),
             textcoords="offset points", fontsize=11.5, color=INK, ha="center",
             arrowprops=dict(arrowstyle="-", color=MUTE, lw=0.8))
ax2.text(23, ideal[n == 23][0] - 0.25, f"斜率 = ln α = ln {a1}", fontsize=11, color=PAL[0], ha="center", va="top",
         rotation=np.degrees(np.arctan(np.log(a1) * 0.62)), rotation_mode="anchor")
ax2.set_xlabel("碳数 n", color=MUTE, fontsize=11.5)
ax2.set_ylabel("ln(Wn / n)", color=MUTE, fontsize=11.5)
ax2.legend(loc="lower left", frameon=False, fontsize=11, labelcolor=INK)
ax2.set_xlim(0, 30.5)
title(ax2, "(b) ASF 图：用直线斜率求 α", "ln(Wn/n) = n ln α + ln[(1−α)²/α]；橙点为示意，不是某一组实验数据")

# ------------------------------------------------------------------ (c) 馏分 vs α
style(ax3)
a = np.linspace(0.30, 0.98, 1500)
cuts = [("CH4", 1, 1), ("C2–C4 低碳烃", 2, 4), ("C5–C11 汽油段", 5, 11), ("C12–C20 柴油段", 12, 20), ("C21+ 蜡", 21, None)]
for (lab, lo, hi), c in zip(cuts, PAL):
    if hi is None:
        f = 1 - sum(W(k, a) for k in range(1, lo))
    else:
        f = sum(W(k, a) for k in range(lo, hi + 1))
    f *= 100
    ax3.plot(a, f, color=c, lw=2.2, label=lab)
    if hi is not None and lo > 1:
        i = int(np.argmax(f))
        ax3.plot(a[i], f[i], "o", color=c, ms=7, mec=SURF, mew=1.5, zorder=5)
        dx = -70 if lo == 12 else 0
        ax3.annotate(f"最高 {f[i]:.0f}%（α = {a[i]:.2f}）", (a[i], f[i]), xytext=(dx, 9), textcoords="offset points",
                     ha="center", fontsize=11, color=INK)
    continue
    pos = {"CH4": 0.36, "C2–C4 低碳烃": 0.30, "C5–C11 汽油段": 0.55, "C12–C20 柴油段": 0.95, "C21+ 蜡": 0.955}[lab]
    j = int(np.argmin(abs(a - pos)))
    off = {"CH4": (6, 4), "C2–C4 低碳烃": (4, 8), "C5–C11 汽油段": (-6, 6), "C12–C20 柴油段": (-10, 6),
           "C21+ 蜡": (-12, -18)}[lab]
    ax3.annotate(lab, (a[j], f[j]), xytext=off, textcoords="offset points", fontsize=12, color=c, fontweight="bold",
                 ha="right" if off[0] < 0 else "left")
ax3.legend(loc="upper left", frameon=False, fontsize=12, labelcolor=INK, ncol=5, handlelength=1.6, columnspacing=1.6)
ax3.set_xlabel("链增长概率 α", color=MUTE, fontsize=11.5)
ax3.set_ylabel("质量分数（%）", color=MUTE, fontsize=11.5)
ax3.set_xlim(0.30, 0.985)
ax3.set_ylim(0, 112)
title(ax3, "(c) 只有一个参数 α：除了 CH4，任何一段馏分都做不到高选择性",
      "C5–C11 最多约 48%、C12–C20 最多约 30%；工业上常把 α 推高多产蜡，再加氢裂化成柴油")

fig.suptitle("ASF（Anderson–Schulz–Flory）分布：每一步以概率 α 增长、以 1−α 终止", x=0.06, y=0.975, ha="left",
             fontsize=19, fontweight="bold", color=INK)
fig.text(0.06, 0.925, "假设 α 与链长无关：链长 n 的摩尔分数 xn = (1−α) α^(n−1)，质量分数 Wn = n (1−α)² α^(n−1)。"
         "(a)(c) 由公式计算，(b) 中的偏差为示意。", fontsize=12, color=MUTE)
out = ROOT / "figures" / "asf_distribution.png"
fig.savefig(out, dpi=200, facecolor=SURF)
fig.savefig(out.with_suffix(".pdf"), facecolor=SURF)
print(out)
