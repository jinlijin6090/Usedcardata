"""Figure 1: composite of the POV-Ray renders with a colour key."""
import json
import os

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from model_lib import COLORS
from plot_style import setup

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REN = os.path.join(ROOT, "renders")


def img(ax, name, title):
    ax.imshow(plt.imread(os.path.join(REN, name + ".png")))
    ax.set_title(title, loc="left")
    ax.axis("off")


def main():
    setup()
    fig = plt.figure(figsize=(14, 9.6))
    top = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.0, 0.95], left=0.01, right=0.99,
                           top=0.96, bottom=0.47, wspace=0.04)
    bot = fig.add_gridspec(1, 3, left=0.01, right=0.99, top=0.42, bottom=0.08, wspace=0.03)
    img(fig.add_subplot(top[0]), "01_overview_oblique",
        "(a) NiCu particle on Ce-doped γ-Al$_2$O$_3$(100) with O vacancies")
    img(fig.add_subplot(top[1]), "02_top_view", "(b) Top view")
    img(fig.add_subplot(top[2]), "03_particle_closeup", "(c) Cu on corners + edges")
    sizes = [("small", "1.3 nm, Cu/Ni = 0.68"), ("medium", "2.0 nm, Cu/Ni = 0.21"),
             ("large", "3.1 nm, Cu/Ni = 0.10")]
    for j, (tag, lab) in enumerate(sizes):
        img(fig.add_subplot(bot[j]), f"04_size_{tag}", f"({'def'[j]}) {lab}")
    names = {"Ni": "Ni", "Cu": "Cu (low-CN sites)", "Al": "Al", "O": "O", "Ce": "Ce (on Al site)",
             "X": "O vacancy"}
    h = [Line2D([], [], marker="o", ls="", ms=11, mfc=COLORS[k], mec="#8a8984", label=v)
         for k, v in names.items()]
    fig.legend(handles=h, loc="lower center", ncol=6, frameon=False, fontsize=11,
               bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 0.055, "(d-f) same support area; Cu fills all corner+edge sites, so the Cu/Ni ratio "
             "needed falls with particle size", ha="center", fontsize=10, color="#52514e")
    fig.savefig(os.path.join(ROOT, "figures", "fig1_model_renders.png"), dpi=200)


if __name__ == "__main__":
    main()
