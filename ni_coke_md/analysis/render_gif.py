#!/usr/bin/env python
"""Render a case trajectory as an animated GIF with matplotlib only (no OpenGL / GPU needed).

    python analysis/render_gif.py CASE_DIR [--frames 150] [--out CASE_DIR/movie.gif]

Left panel : all atoms (Ni green, C dark grey, O red), shaded spheres, perspective.
Right panel: Ni see-through, C-C bonds drawn -> dissolved C, chains and rings become visible.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "md"))
from common import TYPES, iter_dump, unwrap_about_ni  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

STYLE = {TYPES["Ni"]: ("#6fb89c", 1.25), TYPES["C"]: ("#2b2b2b", 0.70), TYPES["O"]: ("#d8432f", 0.66)}
BOND = "#e07b39"


def rot(elev=25, azim=35):
    e, a = np.radians(elev), np.radians(azim)
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(e), -np.sin(e)], [0, np.sin(e), np.cos(e)]])
    return Rx @ Rz


def lighten(hexcol, f):
    c = np.array(matplotlib.colors.to_rgb(hexcol))
    return tuple(c + (1 - c) * f)


def draw(ax, p, types, half, ni_alpha, bonds, cam=60.0, small=1.0):
    order = np.argsort(p[:, 2])                       # far (low z) first
    scale = cam / (cam - p[:, 2])                     # simple perspective
    xy = p[:, :2] * scale[:, None]
    px_per_A = ax.bbox.width / (2 * half)
    if bonds is not None and len(bonds):
        segs = [[xy[i], xy[j]] for i, j in bonds]
        ax.add_collection(LineCollection(segs, colors=BOND, linewidths=3.2, zorder=3.0))
    for i in order:
        col, r = STYLE[types[i]]
        if types[i] != TYPES["Ni"]:
            r *= small
        a = ni_alpha if types[i] == TYPES["Ni"] else 1.0
        s = (2 * r * scale[i] * px_per_A * 72 / ax.figure.dpi) ** 2
        z = 3 + (p[i, 2] + 100) / 1000
        ax.scatter(*xy[i], s=s, c=[col], alpha=a, edgecolors="none", zorder=z)
        if a > 0.5:
            off = 0.28 * r * scale[i]
            ax.scatter(xy[i, 0] - off, xy[i, 1] + off, s=s * 0.18, c=[lighten(col, 0.55)],
                       alpha=0.9, edgecolors="none", zorder=z + 1e-4)
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--frames", type=int, default=150)
    ap.add_argument("--out", default=None)
    ap.add_argument("--px", type=int, default=380, help="panel size in pixels")
    ap.add_argument("--ms", type=int, default=110, help="ms per GIF frame")
    a = ap.parse_args()
    traj = os.path.join(a.case_dir, "traj.lammpstrj")
    steps = [s for s, *_ in iter_dump(traj)]
    keep = set(np.unique(np.linspace(0, len(steps) - 1, min(a.frames, len(steps))).astype(int)))
    R = rot()
    frames_data, half = [], 10.0
    for k, (step, box, ids, types, pos) in enumerate(iter_dump(traj)):
        if k not in keep:
            continue
        p, c = unwrap_about_ni(pos, types, box)
        p = (p - c) @ R.T
        bound = cKDTree(p[types == TYPES["Ni"]]).query(p)[0] < 4.0
        half = max(half, np.abs(p[bound, :2]).max() + 2.5)
        frames_data.append((step, types, p))
    t0 = frames_data[0][0]
    dpi = 100
    imgs = []
    for step, types, p in frames_data:
        fig, axs = plt.subplots(1, 2, figsize=(2 * a.px / dpi, a.px / dpi + 0.45), dpi=dpi,
                                gridspec_kw=dict(wspace=0.02))
        fig.patch.set_facecolor("white")
        fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.86)
        iC = np.where(types == TYPES["C"])[0]
        bonds = []
        if len(iC) > 1:
            pr = cKDTree(p[iC]).query_pairs(1.85, output_type="ndarray")
            bonds = [(iC[i], iC[j]) for i, j in pr]
        draw(axs[0], p, types, half, 1.0, None)
        draw(axs[1], p, types, half, 0.13, bonds, small=0.5)
        nC, nO = int((types == TYPES["C"]).sum()), int((types == TYPES["O"]).sum())
        title = f"t = {(step - t0) * 0.25e-3:7.1f} ps    C: {nC}" + (f"   O: {nO}" if nO else "")
        fig.text(0.01, 0.95, title, fontsize=11, color="#111", family="monospace", va="top")
        fig.text(0.99, 0.95, "Ni green · C dark · O red  |  right: Ni see-through, C–C bonds",
                 fontsize=8, color="#555", ha="right", va="top")
        fig.canvas.draw()
        imgs.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3]))
        plt.close(fig)
    out = a.out or os.path.join(a.case_dir, "movie.gif")
    pal = [im.quantize(colors=128, method=Image.Quantize.MEDIANCUT) for im in imgs]
    pal[0].save(out, save_all=True, append_images=pal[1:], duration=a.ms, loop=0, optimize=True)
    imgs[-1].save(out.replace(".gif", "_last.png"))
    print(f"wrote {out} ({len(imgs)} frames, {os.path.getsize(out)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
