"""画网络规模随重原子上限的变化：节点数、反应数（同一对数纵轴），以及平均每个节点的反应数。

    python plot_scale.py results/subnet_scale_n12.csv -o results/scale_vs_max_heavy.png
"""

import argparse
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter  # noqa: E402

# 参考配色（dataviz 默认调色板，浅色模式）
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
WASH = "#f0efec"
NODES = "#2a78d6"      # 分类色 1
REACTIONS = "#eb6834"  # 分类色 2
RATIO = "#1baf7a"      # 分类色 3（对比度 < 3:1，按规范配直接标注）


def pick_cjk_font():
    for name in ("Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Zen Hei", "PingFang SC",
                 "Microsoft YaHei", "SimHei"):
        if any(f.name == name for f in font_manager.fontManager.ttflist):
            return name
    return None


def compact(v, _pos=None):
    for div, suf in ((1e9, "B"), (1e6, "M"), (1e3, "K")):
        if v >= div:
            x = v / div
            return f"{x:.0f}{suf}" if x >= 10 or x == int(x) else f"{x:.1f}{suf}"
    return f"{v:.0f}"


def style_axes(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(1)
    ax.tick_params(which="both", colors=MUTED, labelcolor=INK_2, labelsize=10, length=0, pad=6)
    ax.grid(True, which="major", axis="y", color=GRID, linewidth=1, linestyle="-")
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-o", "--out", default="scale_vs_max_heavy.png")
    ap.add_argument("--min-heavy", type=int, default=1)
    ap.add_argument("--full-max", type=int, default=5, help="小分子段上限（图中阴影）")
    args = ap.parse_args()

    rows = [r for r in csv.DictReader(open(args.csv)) if int(r["max_heavy"]) >= args.min_heavy]
    n = [int(r["max_heavy"]) for r in rows]
    nodes = [int(r["nodes"]) for r in rows]
    rxns = [int(r["reactions"]) for r in rows]
    ratio = [r / v for r, v in zip(rxns, nodes)]

    font = pick_cjk_font()
    if font:
        plt.rcParams["font.family"] = [font, "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.6), dpi=200,
                                   gridspec_kw={"width_ratios": [1.55, 1], "wspace": 0.28})
    fig.patch.set_facecolor(SURFACE)
    line = dict(linewidth=2, solid_capstyle="round", solid_joinstyle="round")
    marker = dict(marker="o", markersize=7, markeredgecolor=SURFACE, markeredgewidth=2)

    # ---------------- 左：节点数与反应数（同一对数轴）
    style_axes(ax1)
    ax1.axvspan(min(n) - 0.5, args.full_max + 0.5, color=WASH, zorder=0, linewidth=0)
    ax1.text(min(n) - 0.3, 3e8, f"≤ {args.full_max} 个重原子：\n任意两个物种都配对",
             color=INK_2, fontsize=9.5, va="top", ha="left")
    ax1.plot(n, rxns, color=REACTIONS, label="反应数", zorder=3, **line, **marker)
    ax1.plot(n, nodes, color=NODES, label="节点数（物种）", zorder=3, **line, **marker)
    ax1.set_yscale("log")
    ax1.set_ylim(5, 1e9)
    ax1.yaxis.set_major_locator(LogLocator(base=10, numticks=12))
    ax1.yaxis.set_minor_formatter(NullFormatter())
    ax1.yaxis.set_major_formatter(FuncFormatter(compact))
    ax1.set_xticks(n)
    ax1.set_xlim(min(n) - 0.5, max(n) + 1.9)
    ax1.set_xlabel("重原子数上限 N（C + O）", color=INK_2, fontsize=10.5)
    ax1.set_ylabel("数量（对数坐标）", color=INK_2, fontsize=10.5)
    # 直接标注：只标末端
    ax1.annotate(f"{rxns[-1]:,}", (n[-1], rxns[-1]), xytext=(9, 0), textcoords="offset points",
                 va="center", color=INK, fontsize=10)
    ax1.annotate(f"{nodes[-1]:,}", (n[-1], nodes[-1]), xytext=(9, 0), textcoords="offset points",
                 va="center", color=INK, fontsize=10)
    leg = ax1.legend(loc="lower right", frameon=False, fontsize=10, handlelength=1.6)
    for t in leg.get_texts():
        t.set_color(INK)

    # ---------------- 右：平均每个节点的反应数
    style_axes(ax2)
    ax2.axvspan(min(n) - 0.5, args.full_max + 0.5, color=WASH, zorder=0, linewidth=0)
    ax2.plot(n, ratio, color=RATIO, zorder=3, **line, **marker)
    ax2.set_yscale("log")
    ax2.set_ylim(3, 2e4)
    ax2.yaxis.set_major_locator(LogLocator(base=10, numticks=8))
    ax2.yaxis.set_minor_formatter(NullFormatter())
    ax2.yaxis.set_major_formatter(FuncFormatter(compact))
    ax2.set_xticks(n)
    ax2.set_xlim(min(n) - 0.5, max(n) + 1.2)
    ax2.set_xlabel("重原子数上限 N（C + O）", color=INK_2, fontsize=10.5)
    ax2.set_ylabel("反应数 / 节点数（对数坐标）", color=INK_2, fontsize=10.5)
    k = max(range(len(ratio)), key=ratio.__getitem__)
    ax2.annotate(f"峰值 {ratio[k]:,.0f}（N = {n[k]}）", (n[k], ratio[k]), xytext=(8, 8),
                 textcoords="offset points", color=INK, fontsize=10)
    ax2.annotate(f"{ratio[-1]:,.0f}", (n[-1], ratio[-1]), xytext=(9, 0), textcoords="offset points",
                 va="center", color=INK, fontsize=10)
    ax2.annotate(f"{ratio[0]:,.0f}", (n[0], ratio[0]), xytext=(0, -16), textcoords="offset points",
                 ha="center", color=INK, fontsize=10)

    fig.text(0.06, 0.965, "费托反应网络规模随重原子上限的变化", color=INK, fontsize=15, weight="bold",
             va="top")
    fig.text(0.06, 0.905, f"默认配置（≤ 2 个支链）；由 N = {max(n)} 的完整网络按重原子上限截取"
             "（反应中所有物种都 ≤ N）。左：节点数与反应数；右：平均每个节点的反应数",
             color=INK_2, fontsize=10.5, va="top")
    fig.subplots_adjust(left=0.07, right=0.97, top=0.80, bottom=0.12)
    fig.savefig(args.out, facecolor=SURFACE)
    print("saved", args.out)


if __name__ == "__main__":
    main()
