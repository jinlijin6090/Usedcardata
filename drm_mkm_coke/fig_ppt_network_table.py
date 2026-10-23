"""Table image of all reactions of the extended network (R1-R17 published, R18-R21 coke
extension, R22-R24 scenario C) for a 16:9 slide: figures/ppt_network_table.png (3000 x 1700 px).
Barriers for R1-R17 are read from MKM_parameters.xlsx (sheet 'Refined params'); the coke /
scenario-C values are the defaults in mkmc/coke.py, filament.py and s7_filament.py.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, "/home/user/Usedcardata/drm_mkm_coke")
from mkmc.parameters import arrhenius_from_excel
from mkmc.coke import CokeParams
from filament import FilamentParams

OUT = "/home/user/Usedcardata/drm_mkm_coke/figures/ppt_network_table.png"
BLACK, RED, GREEN, GREY = "#1B2A41", "#B5311C", "#1E7A5F", "#6A7179"

rc = arrhenius_from_excel()          # published 34 interleaved barriers
Ef, Er = rc.Ea[0::2], rc.Ea[1::2]
cp, fp = CokeParams(), FilamentParams()
from mkmc.config import R_EV, DEFAULT_T
kT = R_EV * DEFAULT_T
x_eq = np.exp(-fp.dG_gr / kT)
k_diff = fp.D0 * np.exp(-fp.E_D / kT) / (4.0e-9) ** 2


def g3(x):
    return f"{x:.3g}"


pub_eq = [
    "CH$_4$(g) + * → CH$_4$*", "CH$_4$* + * → CH$_3$* + H*", "CH$_3$* + * → CH$_2$* + H*",
    "CH$_2$* + * → CH* + H*", "CH* + * → C* + H*", "2 H* → H$_2$(g) + 2 *",
    "CH* + O* → CHO*", "CHO* + * → CO* + H*", "C* + O* → CO* + *", "CO* → CO(g) + *",
    "O$^v$ + * → O* + v", "CO$_2$(g) + v → CO$_2^v$", "CO$_2^v$ → CO(g) + O$^v$",
    "CO$_2$(g) + O$^v$ → CO$_3^v$", "CO$_3^v$ + * → CO$_2^v$ + O*",
    "H* + O$^v$ → OH$^v$ + *", "OH$^v$ + H* → H$_2$O(g) + v + *",
]
pub_site = ["1 · NSm"] + ["Zmm · NSm"] * 8 + ["1 · NSm", "Zms · NIs", "1 · NSs", "1 · NSs", "1 · NSs",
                                              "Zms · NIs", "Zsm · NIm", "Zsm · NIm"]
pub_role = ["CH$_4$ adsorption", "1st C–H scission", "2nd C–H scission", "3rd C–H scission", "4th C–H scission",
            "H$_2$ desorption", "CH oxidation", "CHO decomposition", "C* oxidation", "CO desorption",
            "interfacial O transport", "CO$_2$ adsorption (support)", "CO$_2$ dissociation (support)",
            "carbonate formation", "carbonate decomposition", "H spillover (OH formation)", "H$_2$O formation"]

rows = []
for i in range(17):
    rows.append((f"R{i + 1}", pub_eq[i], pub_site[i], "published", pub_role[i],
                 f"ΔH‡f = {g3(Ef[i])},  ΔH‡r = {g3(Er[i])}"))
rows += [
    ("R18", "C* + C* $\\rightleftharpoons$ C$_2$* + *", "Zmm · NSm", "coke A + B", "C–C coupling, coke nucleation",
     f"E‡f = {cp.Ea_cc:.2f} (scanned 1.0–2.8; cases 1.5, 2.2),  ΔH = {str(f'{cp.dH_cc:+.2f}').replace('-', '−')},  A = k$_B$T/h"),
    ("R19", "C$_2$* + O* $\\rightleftharpoons$ CO* + C*", "Zmm · NSm", "coke A + B", "oxidative scission of the dimer",
     f"E‡f = {cp.Ea_ox:.2f};  K$_{{19}}$ = K$_9$ / K$_{{18}}$ (thermodynamic closure with R9)"),
    ("R20", "C$_2$* → C(dep) + *", "1 · NSm", "coke A (flux)", "irreversible growth, C off-lattice",
     f"E‡f = {cp.Ea_dep:.2f},  k$_r$ = 0 (off-lattice sink = coke flux)"),
    ("R20b", "C$_2$* + * → 2 Cg*", "Zmm · NSm", "coke B (blocking)", "graphitic C stays on Ni sites",
     f"E‡f = {cp.Ea_dep:.2f} (as R20);  ΔH(C* → Cg*) = {str(f'{cp.dH_g:+.2f}').replace('-', '−')} per C (−0.40 in C)"),
    ("R21", "Cg* + O* $\\rightleftharpoons$ CO* + *", "Zmm · NSm", "coke B (blocking)", "slow gasification of graphitic C",
     f"E‡f = {cp.Ea_gas:.2f} (1.6–2.0 scanned);  K$_{{21}}$ = K$_9$ / K$_g$"),
    ("R22", "C* + bulk site $\\rightleftharpoons$ C$_b$ + *", "k$_f$ · NSm · θ$_{C*}$(1 − x)", "scenario C", "dissolution into the Ni bulk",
     f"E‡f = {fp.Ea_sol:.2f},  ΔH = {fp.dH_sol:+.2f};  bulk sites N$_b$ = Ni atoms\nx$_{{eq}}$ = exp(−ΔG$_{{gr}}$/kT), ΔG$_{{gr}}$ = {fp.dG_gr:.2f}  (x$_{{eq}}$ = {100*x_eq:.2f} % at 973 K)"),
    ("R23", "C$_b$ → C$_{fil}$ (rear filament)", "f · (D$_C$/r$_m^2$) · N$_b$ · (x − x$_{eq}$)", "scenario C",
     "diffusion-limited precipitation",
     f"D$_C$ = 1.2×10$^{{-5}}$ exp(−{fp.E_D:.2f}/kT) m²/s;  k = D$_C$/r$_m^2$ = {k_diff/1e4:.1f}×10$^{{4}}$ s$^{{-1}}$ (r$_m$ = 4 nm)\nJ$_{{nuc}}$ = A exp(−E$_{{nuc}}$/kT)(S−1)²,  E$_{{nuc}}$ = {fp.E_nuc:.2f} (ODE) / 3.4 (population model)"),
    ("R24", "C* + Cg(edge)* → 2 Cg*", "Zmm · NSm · θ$_{edge}$", "scenario C",
     "perimeter-limited edge growth",
     f"E‡f = {fp.Ea_edge:.2f};  θ$_{{edge}}$ = 2√π √(n$_{{isl}}$ θ$_{{Cg}}$);  island nucleation E$_{{isl}}$ = 3.6 (population)"),
]
assert len(rows) == 25

W, H = 15.0, 8.5
fig = plt.figure(figsize=(W, H), dpi=200)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

ax.text(0.35, 8.24, "Elementary steps of the extended Ni/CeO$_2$ DRM microkinetic model",
        fontsize=16, weight="bold", color=BLACK, va="center")
ax.text(14.65, 8.24, "energies in eV, T = 973 K", fontsize=11, color=GREY, va="center", ha="right")

# column layout: (left x, header, align)
cols = [(0.35, "ID", "left"), (0.95, "Elementary step", "left"), (4.2, "Site · coordination factor", "left"),
        (6.45, "Origin", "left"), (7.85, "Role", "left"), (10.15, "Base parameters (eV unless noted)", "left")]
XR = 14.65
y_top, rh, hh = 7.95, 0.245, 0.34
# header
ax.add_patch(Rectangle((0.3, y_top - hh), XR - 0.3 + 0.05, hh, facecolor=BLACK, edgecolor="none", zorder=1))
for x, h, al in cols:
    ax.text(x, y_top - hh / 2, h, fontsize=11.5, weight="bold", color="white", va="center", ha=al, zorder=3)

COL = {"published": BLACK, "coke": RED, "scenario": GREEN}
FS = 10.4
y = y_top - hh
prev_group = None
for k, (rid, eq, site, origin, role, par) in enumerate(rows):
    group = "published" if origin == "published" else ("coke" if origin.startswith("coke") else "scenario")
    if group != prev_group and prev_group is not None:
        ax.plot([0.3, XR + 0.05], [y, y], color=COL[group], lw=1.4, zorder=2)
    prev_group = group
    c = COL[group]
    rhk = rh * (1.75 if "\n" in par else 1.0)
    if k % 2 == 0:
        ax.add_patch(Rectangle((0.3, y - rhk), XR - 0.3 + 0.05, rhk, facecolor="#F3F5F8", edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((0.3, y - rhk), 0.07, rhk, facecolor=c, edgecolor="none", zorder=1))
    yc = y - rhk / 2
    ax.text(cols[0][0] + 0.1, yc, rid, fontsize=FS, weight="bold", color=c, va="center", zorder=3)
    ax.text(cols[1][0], yc, eq, fontsize=FS, color=BLACK, va="center", zorder=3)
    ax.text(cols[2][0], yc, site, fontsize=FS - 0.8, color=BLACK, va="center", zorder=3)
    ax.text(cols[3][0], yc, origin, fontsize=FS - 0.8, color=c, va="center", zorder=3)
    ax.text(cols[4][0], yc, role, fontsize=FS - 0.8, color=GREY, va="center", zorder=3)
    ax.text(cols[5][0], yc, par, fontsize=FS - 1.0, color=BLACK, va="center", zorder=3, linespacing=1.25)
    y -= rhk
ax.plot([0.3, XR + 0.05], [y, y], color=BLACK, lw=1.2, zorder=2)

foot = ("Rate of step j:  r$_j$ = k$_j$ · (coordination × site count / N$_A$) · Π θ · Π p.   NSm / NSs = surface metal / support (O$^v$) sites, "
        "NIm / NIs = interfacial metal / support sites, Zmm, Zms, Zsm = coordination numbers (mkmc/network.py).\n"
        "R1–R17: forward / reverse activation enthalpies from MKM_parameters.xlsx, sheet 'Refined params' (Naveen et al., ACS Catal. 2026, 16, 2535).\n"
        "R18–R21: CokeParams defaults (mkmc/coke.py); R18 + R19 = R9 and R21 close the thermodynamic cycle with R9.   "
        "R22–R24: FilamentParams defaults (filament.py);\n"
        "S = x / x$_{eq}$ = bulk-carbon supersaturation, f = fraction of particles carrying a filament; population-model barriers E$_{isl}$ = 3.6, E$_{nuc}$ = 3.4 from s7_filament.py.")
ax.text(0.35, y - 0.1, foot, fontsize=8.6, color=GREY, va="top", linespacing=1.35)

fig.savefig(OUT, dpi=200, facecolor="white")
print("saved", OUT, "bottom y =", round(y, 2))

if "--check" in sys.argv:
    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    # check that no cell text runs into the next column or off the canvas
    texts = [t for t in ax.texts]
    bad = 0
    for t in texts:
        bb = t.get_window_extent(rend)
        if bb.x1 > W * 200 - 10 or bb.y0 < 0:
            print("  OUTSIDE:", t.get_text()[:50]); bad += 1
    # per-row column overflow
    for t in texts:
        x0 = t.get_position()[0]
        nxt = [c[0] for c in cols if c[0] > x0 + 0.01]
        if nxt and t.get_fontsize() <= FS and t.get_fontsize() > 9:
            bb = t.get_window_extent(rend)
            if bb.x1 / 200 > nxt[0] - 0.08:
                print(f"  OVERFLOW into next column: {t.get_text()[:50]} ends at {bb.x1/200:.2f} > {nxt[0]:.2f}"); bad += 1
    print("checked", len(texts), "texts; problems:", bad)
