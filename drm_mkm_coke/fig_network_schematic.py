"""Schematic of the extended reaction network for the slide deck (matplotlib only)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge, FancyBboxPatch, FancyArrowPatch

BLACK, RED, GREEN, GREY = "#1B2A41", "#B5311C", "#1E7A5F", "#6A7179"
METAL, SUPPORT = "#CBD7E4", "#DCEAD6"
FS_SP, FS_L = 15, 12.5

fig, ax = plt.subplots(figsize=(18, 10), dpi=200)
ax.set_xlim(0, 18); ax.set_ylim(0, 10); ax.axis("off")

ax.add_patch(Rectangle((0.8, 0.4), 14.2, 2.2, facecolor=SUPPORT, edgecolor="#7FA27A", lw=2))
cx, cy, R = 8.5, 2.6, 5.0
ax.add_patch(Wedge((cx, cy), R, 0, 180, facecolor=METAL, edgecolor="#6D8299", lw=2))
ax.plot([cx - R - 0.4, cx + R + 0.4], [cy, cy], color=RED, lw=4, ls=(0, (6, 3)), alpha=0.5, zorder=3)
ax.text(1.0, 0.65, "CeO$_2$(111) support", fontsize=16, color="#2F5D3A", weight="bold")
ax.text(cx - R + 0.35, cy + 0.3, "Ni particle", fontsize=16, color="#3B4C63", weight="bold")
ax.text(0.95, 2.3, "interface ring (NI sites)", fontsize=12.5, color=RED, va="top")

def box(x, y, s, color=BLACK, fc="white", fs=FS_SP, h=0.44):
    n = len(s.replace("$", "").replace("_", "").replace("^", "").replace("{", "").replace("}", ""))
    w = 0.19 * n + 0.4
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                facecolor=fc, edgecolor=color, lw=1.6, zorder=5))
    ax.text(x, y, s, ha="center", va="center", fontsize=fs, color=color, zorder=6, weight="bold")

def arrow(p, q, label, color=BLACK, rad=0.0, ls="-", lw=2.0, lpos=0.5, dx=0.0, dy=0.2, fs=FS_L, shrink=14):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=18, color=color, lw=lw, ls=ls,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=shrink, shrinkB=shrink, zorder=4))
    t = lpos
    mx, my = (1 - t) * p[0] + t * q[0], (1 - t) * p[1] + t * q[1]
    vx, vy = q[0] - p[0], q[1] - p[1]; L = np.hypot(vx, vy)
    mx += -vy / L * rad * L * 0.5; my += vx / L * rad * L * 0.5
    ax.text(mx + dx, my + dy, label, fontsize=fs, color=color, ha="center", va="center", zorder=7,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))

# metal species
box(5.3, 5.6, "CH$_4$*"); box(6.6, 5.6, "CH$_3$*"); box(7.9, 5.6, "CH$_2$*"); box(9.2, 5.6, "CH*"); box(10.5, 5.6, "C*", fc="#FFF1EC")
box(5.9, 4.0, "H*"); box(9.4, 3.8, "O*", fc="#EAF2FF"); box(12.3, 4.6, "CO*"); box(10.2, 6.8, "CHO*")
box(11.6, 3.4, "C$_2$*", color=RED, fc="#FFF1EC"); box(7.4, 3.1, "C$_b$ (bulk)", color=GREEN, fc="#EDF7F0")
# island on the top-left surface
ax.add_patch(Rectangle((5.8, 7.15), 1.6, 0.22, facecolor="#333333", edgecolor="black", zorder=5))
ax.text(9.2, 8.3, "Cg* island (encapsulating carbon):", fontsize=13, color=RED, ha="center", weight="bold")
ax.text(9.2, 7.88, "formed by R20b (C$_2$* + * → 2 Cg*), grows by R24, gasified by R21 (Cg* + O* → CO* + *)", fontsize=11.5, color=RED, ha="center")
# filament at the rear, deposited carbon off-lattice
for i in range(5):
    ax.add_patch(Rectangle((15.2, 1.0 + i * 0.24), 1.0, 0.2, facecolor="#8A8F96", edgecolor="#444", zorder=5))
ax.text(15.7, 0.55, "carbon filament (whisker)\ngrows at the rear, blocks no site", fontsize=12, color=GREEN, ha="center", va="center", weight="bold")
box(16.5, 3.4, "C(dep)", color=RED, fc="#FFF1EC")
# support species and gases
box(2.0, 1.5, "H$_2$O(g)", fc="#F3F3F3"); box(4.3, 1.5, "OH$^v$"); box(6.4, 1.5, "v"); box(8.9, 1.5, "O$^v$", fc="#EAF2FF")
box(11.0, 1.5, "CO$_2^v$"); box(12.3, 0.8, "CO$_3^v$"); box(14.2, 1.5, "CO$_2$(g)", fc="#F3F3F3")
box(5.3, 9.0, "CH$_4$(g)", fc="#F3F3F3"); box(2.3, 5.2, "H$_2$(g)", fc="#F3F3F3"); box(15.6, 5.8, "CO(g)", fc="#F3F3F3")

# published steps
arrow((5.3, 9.0), (5.3, 5.6), "R1", dx=-0.35, dy=0)
arrow((5.3, 5.6), (6.6, 5.6), "R2", dy=0.3); arrow((6.6, 5.6), (7.9, 5.6), "R3", dy=0.3)
arrow((7.9, 5.6), (9.2, 5.6), "R4", dy=0.3); arrow((9.2, 5.6), (10.5, 5.6), "R5", dy=0.3)
ax.text(6.6, 6.3, "−H* at each step", fontsize=11.5, color=GREY, ha="center")
arrow((5.9, 4.0), (2.3, 5.2), "R6  2H* → H$_2$", lpos=0.5, dy=0.3)
arrow((9.2, 5.6), (10.2, 6.8), "R7 +O*", lpos=0.5, dx=-0.5, dy=0.1)
arrow((10.2, 6.8), (12.3, 4.6), "R8", lpos=0.4, dx=0.3, dy=0.15)
arrow((10.5, 5.6), (12.3, 4.6), "R9 +O*", lpos=0.45, dx=-0.15, dy=-0.3)
arrow((12.3, 4.6), (15.6, 5.8), "R10", lpos=0.5, dy=0.3)
arrow((8.9, 1.5), (9.4, 3.8), "R11  O transport", lpos=0.3, dx=1.15, dy=0.0, lw=3)
arrow((14.2, 1.5), (11.0, 1.5), "R12 (+v)", dy=-0.3)
arrow((11.0, 1.5), (8.9, 1.5), "R13  −CO(g)", dy=-0.3)
arrow((14.2, 1.5), (12.3, 0.8), "R14 +O$^v$", lpos=0.55, dx=0.6, dy=-0.05)
arrow((12.3, 0.8), (9.4, 3.8), "R15", lpos=0.35, dx=0.4, dy=0.0)
arrow((5.9, 4.0), (4.3, 1.5), "R16  H spillover (+O$^v$)", lpos=0.4, dx=-1.25, dy=0.0)
arrow((4.3, 1.5), (2.0, 1.5), "R17 +H*", dy=0.3)

# coke extension
arrow((10.5, 5.6), (11.6, 3.4), "R18  C* + C* → C$_2$*", color=RED, lpos=0.55, dx=1.1, dy=0.0)
arrow((11.6, 3.4), (12.3, 4.6), "R19 +O*", color=RED, rad=-0.35, lpos=0.5, dx=0.55, dy=-0.05)
arrow((11.6, 3.4), (16.5, 3.4), "R20  C$_2$* → C(dep) + *   [scenario A, flux]", color=RED, ls="--", lpos=0.5, dx=0.3, dy=-0.3)

# scenario C
arrow((10.5, 5.6), (7.4, 3.1), "R22  C* <-> C$_b$ (dissolution)", color=GREEN, lpos=0.6, dx=-1.05, dy=0.15)
arrow((7.4, 3.1), (15.7, 2.3), "R23  precipitation on the filament, diffusion-limited (k = D$_C$/r$_m^2$), needs S = x/x$_{eq}$ > 1", color=GREEN, lpos=0.3, dx=0.0, dy=-0.27, fs=11.5)
arrow((10.5, 5.6), (7.4, 7.25), "R24  edge growth (C* + edge → island)", color=GREEN, rad=-0.12, lpos=0.55, dx=0.6, dy=0.45, fs=11.5)

# legend and title
lx, ly = 12.3, 9.55
for i, (col, txt) in enumerate([(BLACK, "published steps R1–R17 (Naveen et al., ACS Catal. 2026)"), (RED, "coke extension R18–R21"), (GREEN, "growth-type carbon, scenario C (R22–R24)")]):
    ax.plot([lx, lx + 0.5], [ly - i * 0.38, ly - i * 0.38], color=col, lw=3)
    ax.text(lx + 0.65, ly - i * 0.38, txt, fontsize=12.5, color=col, va="center")
ax.text(0.3, 9.6, "Extended reaction network of the Ni/CeO$_2$ DRM microkinetic model", fontsize=20, weight="bold", color=BLACK)

fig.savefig("/home/user/Usedcardata/drm_mkm_coke/figures/ppt_network_schematic.png", bbox_inches="tight", facecolor="white")
print("saved")
