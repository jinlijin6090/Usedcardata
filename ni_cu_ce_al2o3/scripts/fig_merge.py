"""Figure 5: two small NiCu particles -> pooled Cu/Ni ratio -> large particle at that ratio."""
import os

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from merge_ratio import main as build
from model_lib import COLORS
from plot_style import setup

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REN = os.path.join(ROOT, "renders")


def panel(fig, rect, name, title, sub):
    ax = fig.add_axes(rect)
    ax.imshow(plt.imread(os.path.join(REN, name + ".png")))
    ax.axis("off")
    ax.set_title(title, loc="left")
    ax.text(0.5, -0.02, sub, transform=ax.transAxes, ha="center", va="top", fontsize=10,
            color="#52514e")


def main():
    r = build()
    setup()
    fig = plt.figure(figsize=(15, 6.2))
    s, m = r["small"], r["medium"]
    panel(fig, [0.00, 0.30, 0.22, 0.58], "04_size_small", "(a) 1.3 nm particle",
          f"Cu$_{{{s[0]}}}$Ni$_{{{s[1]}}}$   Cu/Ni = {s[0]/s[1]:.2f}")
    panel(fig, [0.25, 0.30, 0.22, 0.58], "04_size_medium", "(b) 2.0 nm particle",
          f"Cu$_{{{m[0]}}}$Ni$_{{{m[1]}}}$   Cu/Ni = {m[0]/m[1]:.2f}")
    panel(fig, [0.60, 0.16, 0.40, 0.80], "05_large_merged_ratio",
          f"(c) 3.1 nm particle at the pooled ratio Cu/Ni = {r['ratio']:.2f}",
          f"Cu$_{{{r['n_cu']}}}$Ni$_{{{r['n'] - r['n_cu']}}}$: "
          f"corners+edges full ({r['on']['corner'] + r['on']['edge']} Cu), "
          f"{r['on']['face100'] + r['on']['face111']} Cu overflow onto (100)/(111) facets")
    fig.text(0.235, 0.58, "+", fontsize=34, ha="center", va="center", color="#0b0b0b")
    fig.text(0.535, 0.62, "pool atoms", fontsize=10.5, ha="center", color="#52514e")
    fig.text(0.535, 0.535,
             f"Cu$_{{{r['cu']}}}$Ni$_{{{r['ni']}}}$\nCu/Ni = {r['ratio']:.2f}",
             fontsize=11, ha="center", va="top", color="#0b0b0b")
    fig.add_artist(plt.Line2D([0.485, 0.59], [0.58, 0.58], color="#0b0b0b", lw=2))
    fig.add_artist(plt.Line2D([0.578, 0.59, 0.578], [0.60, 0.58, 0.56], color="#0b0b0b", lw=2))
    names = {"Ni": "Ni", "Cu": "Cu", "Al": "Al", "O": "O", "Ce": "Ce", "X": "O vacancy"}
    h = [Line2D([], [], marker="o", ls="", ms=11, mfc=COLORS[k], mec="#8a8984", label=v)
         for k, v in names.items()]
    fig.legend(handles=h, loc="lower left", ncol=6, frameon=False, fontsize=11,
               bbox_to_anchor=(0.03, 0.02))
    fig.savefig(os.path.join(ROOT, "figures", "fig5_merge_to_large.png"), dpi=200)


if __name__ == "__main__":
    main()
