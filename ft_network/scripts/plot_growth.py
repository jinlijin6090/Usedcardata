#!/usr/bin/env python3
"""网络组成随碳数的变化: C1 → 加入 C2 → 加入 C3。
三个面板累积画网络, 每个面板只给新增的节点和箭头上色, 之前已有的变灰; 下方是组成统计表。"""
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

import plot_full_network as net

W, H = 43.0, 27.3
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.axis("off")

XOFF = [0.0, 7.7, 21.7]
YOFF = 1.2
for k, xo in enumerate(XOFF, 1):
    net.draw_network(ax, nmax=k, hl=k, xoff=xo, yoff=YOFF, legend=False)

INK, MUTE, HEAD = "#1B2A33", "#6B7880", "#0F4C5C"
plt.rcParams["font.sans-serif"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]

# ------------------------------------------------------------------ 统计
GROUPS = [("反应物吸附", ["1 "]), ("CO 活化（直接解离 + 氢助）", ["2a", "2b", "2c"]), ("除氧 / 水煤气变换", ["3a", "3b"]),
          ("CHx 与 Cn 中间体加氢/脱氢", ["4 ", "5g"]), ("C–C 偶联（碳化物机理）", ["5a"]), ("CO 插入", ["5b"]),
          ("CO 插入中间体加氢", ["5c"]), ("CO 插入中间体 C–O 断裂", ["5d"]), ("烯醇缩合", ["5e"]), ("烯基机理", ["5f"]),
          ("产物脱附 / 再吸附", ["5h"])]
stage = lambda n: max(1, int(n))
rx = Counter()
for r in net.RX_ALL:
    if stage(r["n"]) <= 3:
        g = next(name for name, keys in GROUPS if r["category"][:2] in keys)
        if r["equation"] in ("CH4* ⇌ CH4(g) + *",):
            g = "产物脱附 / 再吸附"
        rx[(g, stage(r["n"]))] += 1
SMALLS = {"H*", "O*", "OH*", "H2O*", "CO*", "CO2*", "COOH*", "HCOO*"}


def sp_type(s):
    if s.sites == 0:
        return "气相分子"
    if s.name in SMALLS:
        return "小分子吸附物种"
    if s.kind == "η3-烯丙基":
        return "η3-烯丙基"
    return "含氧中间体" if s.O else "烃类中间体（CHx、CnHm）"


SP_ROWS = ["烃类中间体（CHx、CnHm）", "含氧中间体", "η3-烯丙基", "小分子吸附物种", "气相分子"]
spc = Counter((sp_type(s), stage(s.n)) for s in net.SPEC_ALL.values() if stage(s.n) <= 3)
tot_sp = [sum(spc[(t, k)] for t in SP_ROWS) for k in (1, 2, 3)]
tot_rx = [sum(rx[(g, k)] for g, _ in GROUPS) for k in (1, 2, 3)]
assert tot_sp == [24, 36, 33] and tot_rx == [35, 81, 77], (tot_sp, tot_rx)

# ------------------------------------------------------------------ 面板标题
titles = ["① C1：只有单体化学", "② 加入 C2：出现 C–C 成键", "③ 加入 C3：进入重复单元"]
for k, xo in enumerate(XOFF):
    cs, cr = sum(tot_sp[:k + 1]), sum(tot_rx[:k + 1])
    ax.text(xo + 0.6, 0.95, titles[k], fontsize=30, fontweight="bold", color=HEAD, va="center")
    sub = f"{tot_sp[k]} 个物种 · {tot_rx[k]} 个基元步骤" if k == 0 else \
        f"新增 {tot_sp[k]} 个物种 · {tot_rx[k]} 个基元步骤   →   累计 {cs} · {cr}"
    ax.text(xo + 0.6, 1.75, sub, fontsize=18, color=MUTE, va="center")
for xo in XOFF[1:]:
    ax.plot([xo - 0.25, xo - 0.25], [0.5, 17.6], color="#DDE3E6", lw=1.5)
ax.text(0.6, 18.15, "彩色 = 这一步新增的节点和箭头；灰色 = 之前已有的部分。箭头旁的小数字是反应编号。",
        fontsize=15, color=MUTE, va="center")


# ------------------------------------------------------------------ 表格
def table(x0, y0, title, rows, data, totals, c0w):
    cols = ["C1", "+ C2", "+ C3", "累计"]
    cw, rh = 1.55, 0.52
    ax.text(x0, y0, title, fontsize=22, fontweight="bold", color=HEAD, va="center")
    y = y0 + 0.85
    for j, c in enumerate(cols):
        ax.text(x0 + c0w + (j + 0.5) * cw, y, c, fontsize=16, fontweight="bold", color=INK, ha="center", va="center")
    ax.plot([x0, x0 + c0w + 4 * cw], [y + 0.33, y + 0.33], color=INK, lw=1.2)
    for i, name in enumerate(rows):
        yy = y + 0.72 + i * rh
        if i % 2 == 0:
            ax.add_patch(FancyBboxPatch((x0, yy - rh / 2), c0w + 4 * cw, rh, boxstyle="square,pad=0", fc="#F3F6F7", ec="none", zorder=0))
        ax.text(x0 + 0.15, yy, name, fontsize=15, color=INK, va="center")
        vals = [data[(name, k)] for k in (1, 2, 3)]
        first = next((j for j, v in enumerate(vals) if v), None)
        for j, v in enumerate(vals + [sum(vals)]):
            is_first = j == first and j > 0
            ax.text(x0 + c0w + (j + 0.5) * cw, yy, str(v) if v else "–", fontsize=15, ha="center", va="center",
                    color=MUTE if not v else INK, fontweight="bold" if (is_first or j == 3) else "normal")
            if is_first:                                   # 该类别首次出现
                ax.add_patch(FancyBboxPatch((x0 + c0w + j * cw + 0.3, yy - 0.2), cw - 0.6, 0.4, boxstyle="round,pad=0,rounding_size=0.08",
                                            fc="none", ec="#D9580D", lw=1.6))
    yy = y + 0.72 + len(rows) * rh
    ax.plot([x0, x0 + c0w + 4 * cw], [yy - rh / 2, yy - rh / 2], color=INK, lw=1.2)
    ax.text(x0 + 0.15, yy + 0.05, "合计", fontsize=16, fontweight="bold", color=INK, va="center")
    for j, v in enumerate(totals + [sum(totals)]):
        ax.text(x0 + c0w + (j + 0.5) * cw, yy + 0.05, str(v), fontsize=16, fontweight="bold", color=INK, ha="center", va="center")
    return yy


TY = 19.4
table(0.6, TY, "基元步骤的组成", [g for g, _ in GROUPS], rx, tot_rx, 5.6)
yb = table(13.6, TY, "物种的组成", SP_ROWS, spc, tot_sp, 4.9)
ax.add_patch(FancyBboxPatch((13.75, yb + 0.75), 0.75, 0.4, boxstyle="round,pad=0,rounding_size=0.08", fc="none", ec="#D9580D", lw=1.6))
ax.text(14.7, yb + 0.95, "= 该类别首次出现", fontsize=15, color=MUTE, va="center")

NX = 26.2
ax.text(NX, TY, "组成怎么变", fontsize=22, fontweight="bold", color=HEAD, va="center")
notes = [("C1", "没有 C–C 键。35 步全是 CO 活化、除氧和 CHx 逐步加氢，产物只有 CH4、H2O、CO2。"),
         ("C2", "首次出现三类链增长：CHx 两两偶联 10 步、CO 插入 4 步、烯醇缩合 1 步。"),
         ("", "CO 插入带来 24 个含氧中间体，它们的加氢和 C–O 断裂占了 81 步中的 52 步。"),
         ("", "首次出现烯烃产物（乙烯），链终止分成烷烃和烯烃两条路。"),
         ("C3", "新增烯基机理 3 步和 η3-烯丙基；偶联从“CHx 两两组合”变成“3 个片段 × 4 种 CHx”共 12 步。"),
         ("", "含氧中间体降到 18 个（片段从 4 个变 3 个），烃类中间体从 10 个增到 12 个。"),
         ("之后", "C4 起每个碳数都重复 C3 的这一套：33 个物种、77 个基元步骤，组成不再变化。")]
for i, (tag, txt) in enumerate(notes):
    yy = TY + 0.95 + i * 0.78
    if tag:
        ax.add_patch(FancyBboxPatch((NX, yy - 0.24), 1.05, 0.48, boxstyle="round,pad=0,rounding_size=0.1", fc=HEAD, ec="none"))
        ax.text(NX + 0.525, yy, tag, fontsize=15, fontweight="bold", color="white", ha="center", va="center")
    ax.text(NX + 1.35, yy, txt, fontsize=15.5, color=INK, va="center")

out = net.ROOT / "figures" / "ft_network_growth_C1-C3.png"
fig.savefig(out, dpi=180, facecolor="white")
fig.savefig(out.with_suffix(".pdf"), facecolor="white")
print(out, dict(rx), dict(spc))
