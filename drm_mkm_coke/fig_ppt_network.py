"""Clean schematic of the EXTENDED reaction network for a 16:9 slide.

Output: figures/ppt_network.png  (3000 x 1700 px, dpi 200, white background)
matplotlib only.  Layout coordinates are in inches (1 unit = 200 px).
Run with --check to print label/box/arrow collisions.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge, Ellipse, FancyBboxPatch, FancyArrowPatch

OUT = "/home/user/Usedcardata/drm_mkm_coke/figures/ppt_network.png"
BLACK, RED, GREEN, GREY, ORANGE = "#1B2A41", "#B5311C", "#1E7A5F", "#6A7179", "#C97A1F"
METAL_FC, METAL_EC = "#D6E0EC", "#6D8299"
SUPP_FC, SUPP_EC = "#E3EFDC", "#7FA27A"
GAS_FC = "#F2F2F2"
FS_SP, FS_L, FS_NOTE = 14, 11.5, 10.5

W, H = 15.0, 8.5
fig = plt.figure(figsize=(W, H), dpi=200)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

# ------------------------------------------------------------------ catalyst
CX, CY, R = 7.2, 2.3, 5.6
ax.add_patch(Rectangle((0.3, 0.3), 12.6, 2.0, facecolor=SUPP_FC, edgecolor=SUPP_EC, lw=2, zorder=1))
# interface ring: a flattened ellipse at the particle foot, half hidden by the particle
ax.add_patch(Ellipse((CX, CY), 2 * R + 0.5, 0.6, facecolor="#F8DDB6", edgecolor=ORANGE, lw=2,
                     ls=(0, (5, 3)), zorder=2))
ax.add_patch(Wedge((CX, CY), R, 0, 180, facecolor=METAL_FC, edgecolor=METAL_EC, lw=2, zorder=3))
ax.text(0.5, 0.42, "CeO$_2$ (support)", fontsize=15, color="#2F5D3A", weight="bold", zorder=20)
ax.text(4.6, 6.35, "Ni (metal)", fontsize=15, color="#3B4C63", weight="bold", zorder=20)
ax.text(2.25, 1.97, "interface ring (NI sites)", fontsize=FS_NOTE, color=ORANGE, ha="center", va="center", zorder=20,
        style="italic")

BOXES, TEXTS, ARROWS = {}, [], []


def box(name, x, y, s, color=BLACK, fc="white", fs=FS_SP, h=0.44, w=None):
    n = len(s.replace("$", "").replace("_", "").replace("^", "").replace("{", "").replace("}", ""))
    w = w or 0.17 * n + 0.36
    p = FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                       facecolor=fc, edgecolor=color, lw=1.6, zorder=8)
    ax.add_patch(p)
    t = ax.text(x, y, s, ha="center", va="center", fontsize=fs, color=color, zorder=9, weight="bold")
    BOXES[name] = dict(patch=p, c=(x, y), w=w, h=h, text=t)
    return (x, y)


def edge_point(name, toward, pad=0.05):
    """Point on the (padded) box boundary along the ray centre -> toward."""
    b = BOXES[name]; cx, cy = b["c"]; hw, hh = b["w"] / 2 + pad, b["h"] / 2 + pad
    dx, dy = toward[0] - cx, toward[1] - cy
    if dx == 0 and dy == 0:
        return (cx, cy)
    tx = hw / abs(dx) if dx else np.inf
    ty = hh / abs(dy) if dy else np.inf
    t = min(tx, ty)
    return (cx + dx * t, cy + dy * t)


def label(x, y, s, color=BLACK, fs=FS_L, ha="center", va="center", bbox=True, **kw):
    t = ax.text(x, y, s, fontsize=fs, color=color, ha=ha, va=va, zorder=12,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.9) if bbox else None, **kw)
    TEXTS.append(t)
    return t


def arrow(a, b, lab, color=BLACK, rad=0.0, ls="-", lw=2.0, lpos=0.5, dx=0.0, dy=0.0, fs=FS_L, mscale=18):
    ca = BOXES[a]["c"] if isinstance(a, str) else a
    cb = BOXES[b]["c"] if isinstance(b, str) else b
    # arc3 control point: midpoint + rad * (dy, -dx)
    mx, my = (ca[0] + cb[0]) / 2, (ca[1] + cb[1]) / 2
    vx, vy = cb[0] - ca[0], cb[1] - ca[1]
    ctrl = (mx + rad * vy, my - rad * vx)
    pa = edge_point(a, ctrl if rad else cb) if isinstance(a, str) else a
    pb = edge_point(b, ctrl if rad else ca) if isinstance(b, str) else b
    p = FancyArrowPatch(pa, pb, arrowstyle="-|>", mutation_scale=mscale, color=color, lw=lw, ls=ls,
                        connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=2, zorder=6)
    ax.add_patch(p)
    ARROWS.append((p, a, b, lab))
    if lab:
        t = lpos
        # quadratic Bezier point at parameter t (control recomputed from the clipped endpoints)
        mx, my = (pa[0] + pb[0]) / 2, (pa[1] + pb[1]) / 2
        vx, vy = pb[0] - pa[0], pb[1] - pa[1]
        c = (mx + rad * vy, my - rad * vx)
        lx = (1 - t) ** 2 * pa[0] + 2 * (1 - t) * t * c[0] + t ** 2 * pb[0]
        ly = (1 - t) ** 2 * pa[1] + 2 * (1 - t) * t * c[1] + t ** 2 * pb[1]
        label(lx + dx, ly + dy, lab, color=color, fs=fs)
    return p


# ------------------------------------------------------------------ species
# metal
box("CH4*", 4.0, 5.3, "CH$_4$*"); box("CH3*", 5.4, 5.3, "CH$_3$*"); box("CH2*", 6.8, 5.3, "CH$_2$*")
box("CH*", 8.2, 5.3, "CH*"); box("C*", 9.6, 5.3, "C*", fc="#FFF1EC")
box("H*", 4.9, 3.7, "H*"); box("O*", 7.2, 3.5, "O*", fc="#EAF2FF")
box("CHO*", 8.8, 6.8, "CHO*"); box("CO*", 10.8, 5.9, "CO*")
box("C2*", 11.2, 4.4, "C$_2$*", color=RED, fc="#FFF1EC")
box("Cb", 8.6, 3.05, "C$_b$ (bulk)", color=GREEN, fc="#EDF7F0")
# support
box("OHv", 3.7, 1.5, "OH$^v$"); box("Ov", 7.2, 1.5, "O$^v$", fc="#EAF2FF"); box("CO2v", 9.2, 1.6, "CO$_2^v$")
box("v", 5.4, 0.75, "v"); box("CO3v", 10.0, 0.65, "CO$_3^v$")
# gas phase
box("CH4g", 2.3, 6.9, "CH$_4$(g)", fc=GAS_FC); box("H2g", 1.2, 4.3, "H$_2$(g)", fc=GAS_FC)
box("COg", 12.4, 6.25, "CO(g)", fc=GAS_FC); box("CO2g", 12.0, 1.15, "CO$_2$(g)", fc=GAS_FC)
box("H2Og", 1.6, 1.0, "H$_2$O(g)", fc=GAS_FC)
box("Cdep", 14.1, 4.45, "C(dep)", color=RED, fc="#FFF1EC")

# coke-extension shaded region (translucent, dashed red edge)
ax.add_patch(FancyBboxPatch((10.6, 2.85), 4.1, 2.15, boxstyle="round,pad=0.02,rounding_size=0.15",
                            facecolor="#FBE9E4", edgecolor=RED, lw=1.5, ls=(0, (4, 3)), alpha=0.75, zorder=4))
label(14.6, 5.18, "coke extension (R18–R21)", color=RED, fs=FS_L, weight="bold", bbox=False, ha="right")

# graphene island (Cg*) on the right flank of the particle, 7-20 degrees
ax.add_patch(Wedge((CX, CY), R + 0.24, 7, 20, width=0.42, facecolor="#333333", edgecolor="black", lw=1, zorder=7))
ISL_IN = lambda deg, r=R - 0.2: (CX + r * np.cos(np.radians(deg)), CY + r * np.sin(np.radians(deg)))
label(13.9, 3.45, "Cg* graphene island\n(formed by R20b, grows\nby R24, gasified by R21;\nblocks Ni sites)",
      color=RED, fs=FS_NOTE, bbox=False)

# whisker / filament at the rear (right) of the particle, growing away from the surface
ang = np.radians(-14); d = np.array([np.cos(ang), np.sin(ang)]); nrm = np.array([-d[1], d[0]])
root = np.array([12.75, 2.5])
for i in range(6):
    off = (i - 2.5) * 0.105
    p0 = root + nrm * off; p1 = p0 + d * 1.55
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#5A6068", lw=3.6, solid_capstyle="butt", zorder=7)
label(13.8, 1.3, "carbon filament (whisker):\ngrows at the rear,\nblocks no Ni site", color=GREEN, fs=10, bbox=False)

# ------------------------------------------------------------------ published steps (black)
arrow("CH4g", "CH4*", "R1", lpos=0.5, dx=0.25, dy=0.15)
arrow("CH4*", "CH3*", "R2", dy=0.28); arrow("CH3*", "CH2*", "R3", dy=0.28)
arrow("CH2*", "CH*", "R4", dy=0.28); arrow("CH*", "C*", "R5", dy=0.28)
label(6.1, 4.62, "R2–R5: one H* released per step", color=GREY, fs=FS_NOTE, bbox=False)
arrow("H*", "H2g", "R6  2H* → H$_2$", lpos=0.5, dy=0.3)
arrow("CH*", "CHO*", "R7 +O*", lpos=0.5, dx=-0.55, dy=0.0)
arrow("CHO*", "CO*", "R8", lpos=0.5, dx=0.0, dy=0.25)
arrow("C*", "CO*", "R9 +O*", lpos=0.45, dx=-0.35, dy=0.18)
arrow("CO*", "COg", "R10", lpos=0.5, dy=0.25)
arrow("Ov", "O*", "R11  O transport", lpos=0.5, dx=-0.95, dy=0.05, lw=3)
arrow("CO2g", "CO2v", "R12 (+v)", lpos=0.5, dy=0.27)
arrow("CO2v", "Ov", "R13  −CO(g)", lpos=0.5, dy=-0.3)
arrow("CO2g", "CO3v", "R14 +O$^v$", lpos=0.5, dx=0.0, dy=-0.28)
arrow("CO3v", "CO2v", "R15 (+* → O*)", lpos=0.5, dx=-0.75, dy=-0.25)
arrow("H*", "OHv", "R16  H spillover (+O$^v$)", lpos=0.55, dx=0.95, dy=0.25)
arrow("OHv", "H2Og", "R17 +H*", lpos=0.5, dy=0.3)

# ------------------------------------------------------------------ coke extension (red)
arrow("C*", "C2*", "R18", color=RED, lpos=0.5, dx=0.05, dy=0.17)
arrow("C2*", "CO*", "R19 +O*", color=RED, lpos=0.45, dx=0.15, dy=0.0)
arrow("C2*", ISL_IN(13), "R20b", color=RED, lpos=0.5, dx=-0.15, dy=-0.25)
arrow("C2*", "Cdep", "R20 (flux)", color=RED, ls="--", lpos=0.62, dy=0.27)
arrow(ISL_IN(19), "CO*", "R21  Cg* + O* → CO*", color=RED, rad=0.45, lpos=0.5, dx=0.95, dy=0.2)

# ------------------------------------------------------------------ scenario C (green)
arrow("C*", "Cb", "R22  dissolution", color=GREEN, lpos=0.5, dx=-0.8, dy=0.1)
arrow("Cb", tuple(root), None, color=GREEN)
label(9.4, 2.5, "R23  diffusion-limited precipitation\non the rear filament, k = D$_C$/r$_m^2$", color=GREEN, fs=FS_NOTE)
arrow("C*", ISL_IN(9), "R24  C* + island edge", color=GREEN, rad=0.4, lpos=0.5, dx=0.15, dy=-0.3)

# ------------------------------------------------------------------ legend
lx, ly, dyl = 9.95, 8.25, 0.27
entries = [(BLACK, ["published steps R1–R17 (Naveen et al., ACS Catal. 2026)"]),
           (RED, ["coke extension R18–R21:  C* + C* → C$_2$* + * (R18),",
                  "C$_2$* + O* → CO* + C* (R19),  C$_2$* → C(dep) + * (R20, flux) or",
                  "C$_2$* + * → 2 Cg* (R20b, blocking),  Cg* + O* → CO* + * (R21)"]),
           (GREEN, ["scenario C, growth-type carbon:  C* $\\rightleftharpoons$ C$_b$ (R22, dissolution),",
                    "C$_b$ → rear filament (R23),  C* + island edge → island (R24)"])]
y = ly
for col, lines in entries:
    ax.plot([lx, lx + 0.45], [y, y], color=col, lw=3)
    for ln in lines:
        label(lx + 0.6, y, ln, color=col, fs=10, ha="left", bbox=False)
        y -= dyl

fig.savefig(OUT, dpi=200, facecolor="white")
print("saved", OUT)

# ------------------------------------------------------------------ collision check
if "--check" in sys.argv:
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    items = []
    for name, b in BOXES.items():
        items.append((f"box:{name}", b["patch"].get_window_extent(rend)))
    for t in TEXTS:
        items.append((f"txt:{t.get_text()[:28]}", t.get_window_extent(rend)))
    def ov(a, b, pad=2):
        return not (a.x1 < b.x0 + pad or b.x1 < a.x0 + pad or a.y1 < b.y0 + pad or b.y1 < a.y0 + pad)
    print("--- text/box overlaps")
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if ov(items[i][1], items[j][1]):
                print("  OVERLAP", items[i][0], "<->", items[j][0])
    print("--- arrows passing through boxes/labels")
    for p, a, b, lab in ARROWS:
        polys = p.get_path().to_polygons(closed_only=False)
        verts = np.vstack(polys)
        disp = ax.transData.transform(verts)
        for name, bb in items:
            if name in (f"box:{a}", f"box:{b}", f"txt:{(lab or '')[:28]}"):
                continue
            inside = ((disp[:, 0] > bb.x0) & (disp[:, 0] < bb.x1) & (disp[:, 1] > bb.y0) & (disp[:, 1] < bb.y1))
            if inside.any():
                print(f"  ARROW {lab or (a, b)} crosses {name}")
    print("--- out-of-canvas")
    for name, bb in items:
        if bb.x1 > W * 200 or bb.x0 < 0 or bb.y1 > H * 200 or bb.y0 < 0:
            print("  OUTSIDE", name)
