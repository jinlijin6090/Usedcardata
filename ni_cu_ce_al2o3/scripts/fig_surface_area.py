"""Figure 6 (+ individual panels): specific surface area and exposed-Ni statistics."""
import csv
import os

import matplotlib.pyplot as plt
import numpy as np

from plot_style import SERIES, setup

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
ROWS = list(csv.DictReader(open(os.path.join(ROOT, "data", "surface_area.csv"))))
SCEN = ["pure Ni", "Cu on corners+edges", "Cu/Ni = 0.10", "Cu/Ni = 0.27"]
STYLE = {s: (SERIES[i], m) for i, (s, m) in enumerate(zip(SCEN, "osD^"))}
RHO_NI = 8.908  # g/cm3


def series(scen, key):
    r = [x for x in ROWS if x["series"] == "size scan" and x["scenario"] == scen]
    return np.array([float(x["d_nm"]) for x in r]), np.array([float(x[key]) for x in r])


def lines(a, key, scale=1.0, scen=SCEN):
    for s in scen:
        d, y = series(s, key)
        c, m = STYLE[s]
        a.plot(d, y * scale, "-", color=c, marker=m, ms=5, label=s)
    a.set_xlabel("Particle diameter (nm)")


def p_ssa(a):
    lines(a, "SSA_metal_m2g", scen=["pure Ni"])
    d = np.linspace(1.4, 10, 100)
    a.plot(d, 6000 / (RHO_NI * d), "--", color="#52514e", lw=1.5, label="sphere, 6/(ρd)")
    a.set_ylabel("Metal specific surface area (m²/g metal)")
    a.set_title("Metal specific surface area", loc="left")
    a.legend(frameon=False)


def p_disp(a):
    lines(a, "dispersion", 100, scen=["pure Ni"])
    a.set_ylabel("Dispersion: exposed / all atoms (%)")
    a.set_title("Dispersion", loc="left")


def p_ssa_ni(a):
    lines(a, "SSA_Ni_m2g")
    a.set_ylabel("Exposed Ni surface area (m²/g metal)")
    a.set_title("Exposed Ni surface area", loc="left")
    a.legend(frameon=False, fontsize=9)


def p_surf_ni(a):
    lines(a, "surface_Ni_frac", 100)
    a.set_ylabel("Ni share of the exposed sites (%)")
    a.set_ylim(0, 105)
    a.set_title("Surface Ni fraction", loc="left")
    a.legend(frameon=False, fontsize=9, loc="lower right")


def p_ni_exp(a):
    lines(a, "Ni_exposed_frac", 100)
    a.set_ylabel("Ni atoms exposed / all Ni atoms (%)")
    a.set_title("Ni exposure (Ni dispersion)", loc="left")
    a.legend(frameon=False, fontsize=9)


def p_bars(a):
    picks = [("1.3 nm render", "Cu on corners+edges", "1.3 nm\nCu/Ni 0.68"),
             ("2.0 nm render", "Cu on corners+edges", "2.0 nm\nCu/Ni 0.21"),
             ("3.1 nm render", "Cu on corners+edges", "3.1 nm\nCu/Ni 0.10"),
             ("3.1 nm merged-ratio render", "Cu/Ni = 0.27", "3.1 nm\nCu/Ni 0.27")]
    labs, ni_e, cu_e, bur = [], [], [], []
    for ser, sc, lab in picks:
        r = next(x for x in ROWS if x["series"] == ser and x["scenario"] == sc)
        n = int(r["N"])
        labs.append(lab)
        ni_e.append(100 * int(r["N_Ni_exposed"]) / n)
        cu_e.append(100 * int(r["N_Cu_exposed"]) / n)
        bur.append(100 - ni_e[-1] - cu_e[-1])
    x = np.arange(len(labs))
    a.bar(x, ni_e, 0.6, color=SERIES[0], edgecolor="white", lw=2, label="exposed Ni")
    a.bar(x, cu_e, 0.6, bottom=ni_e, color=SERIES[1], edgecolor="white", lw=2, label="exposed Cu")
    a.bar(x, bur, 0.6, bottom=np.add(ni_e, cu_e), color="#d6d5cf", edgecolor="white", lw=2,
          label="interior / interface")
    for i in range(len(x)):
        a.text(x[i], ni_e[i] / 2, f"{ni_e[i]:.0f}%", ha="center", va="center", color="white",
               fontsize=9, fontweight="bold")
        a.text(x[i], ni_e[i] + cu_e[i] / 2, f"{cu_e[i]:.0f}%", ha="center", va="center",
               color="#0b0b0b", fontsize=9)
    a.set_xticks(x, labs, fontsize=9)
    a.set_ylabel("Share of all metal atoms (%)")
    a.set_ylim(0, 100)
    a.set_title("Rendered particles: where the atoms are", loc="left")
    a.legend(frameon=False, fontsize=9, loc="upper right")
    a.grid(axis="x", visible=False)


PANELS = [("6a_metal_SSA", p_ssa), ("6b_dispersion", p_disp), ("6c_exposed_Ni_area", p_ssa_ni),
          ("6d_surface_Ni_fraction", p_surf_ni), ("6e_Ni_exposure", p_ni_exp),
          ("6f_rendered_particles", p_bars)]


def main():
    setup()
    for name, f in PANELS:
        fig, a = plt.subplots(figsize=(5.6, 4.3))
        f(a)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, f"fig{name}.png"), dpi=220)
        plt.close(fig)
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.2))
    for a, (name, f) in zip(ax.ravel(), PANELS):
        f(a)
        a.set_title(f"({name[1]}) " + a.get_title(loc="left"), loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig6_surface_area_Ni_exposure.png"), dpi=200)


if __name__ == "__main__":
    main()
