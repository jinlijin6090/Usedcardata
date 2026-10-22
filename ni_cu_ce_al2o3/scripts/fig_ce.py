"""Figure 3: Ce doping of gamma-Al2O3 and surface O vacancies (MACE-MP-0 energies)."""
import json
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom

from plot_style import BLUE_RAMP, SERIES, setup

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = json.load(open(os.path.join(ROOT, "data", "ce_vacancy_mace.json")))
KB = 8.617e-5
Z_CAT = 2  # surface cation neighbours of a top-layer O on the (100) termination


def evf_by_k():
    """Mean E_vf for O with k Ce neighbours (all Ce contents pooled) and the k=0 reference."""
    runs = R["runs"]
    ks = sorted({r["k_ce"] for r in runs})
    return {k: np.mean([r["E_vf"] for r in runs if r["k_ce"] == k]) for k in ks}


def vac_enhancement(x_ce, T, E):
    """[V_O](x)/[V_O](0) for a random surface Ce fraction x (dilute vacancies, Boltzmann)."""
    e0 = E[0]
    kmax = max(E)
    w = sum(binom.pmf(k, Z_CAT, x_ce) * np.exp(-(E[min(k, kmax)] - e0) / (KB * T))
            for k in range(Z_CAT + 1))
    return w


def main():
    setup()
    runs = R["runs"]
    E = evf_by_k()
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))

    # (a) every computed vacancy
    a = ax[0]
    n_ces = sorted({r["n_ce"] for r in runs})
    for j, n in enumerate(n_ces):
        rr = [r for r in runs if r["n_ce"] == n]
        k = np.array([r["k_ce"] for r in rr]) + (j - 1.5) * 0.07
        a.plot(k, [r["E_vf"] for r in rr], "o", ms=8, color=SERIES[j], mec="white", mew=1.2,
               label=f"{100*rr[0]['ce_frac_top']:.1f} % Ce (surface)")
    ks = sorted(E)
    a.plot(ks, [E[k] for k in ks], "-", color="#52514e", lw=1.5, label="mean")
    a.set_xticks(ks)
    a.set_xlabel("Ce neighbours of the removed O")
    a.set_ylabel(r"$E_{vf}$ (eV, vs. ½ O$_2$)")
    a.set_title("(a) O-vacancy formation energy", loc="left")
    a.legend(frameon=False, fontsize=8.5)

    # (b) lowest / mean E_vf vs Ce content of the slab
    a = ax[1]
    frac = [100 * [r for r in runs if r["n_ce"] == n][0]["ce_frac_top"] for n in n_ces]
    mins = [min(r["E_vf"] for r in runs if r["n_ce"] == n) for n in n_ces]
    means = [np.mean([r["E_vf"] for r in runs if r["n_ce"] == n]) for n in n_ces]
    a.plot(frac, mins, "-o", color=SERIES[0], ms=6, label="easiest O site")
    a.plot(frac, means, "-s", color=SERIES[1], ms=6, label="mean of sampled sites")
    a.set_xlabel("Ce / (Al+Ce) in the surface layer (%)")
    a.set_ylabel(r"$E_{vf}$ (eV)")
    a.set_title("(b) Vacancy energy vs. Ce loading", loc="left")
    a.legend(frameon=False, fontsize=9)

    # (c) equilibrium vacancy enhancement vs Ce loading
    a = ax[2]
    x = np.linspace(0, 0.30, 61)
    for j, T in enumerate([773, 873, 973, 1073]):
        a.plot(100 * x, [vac_enhancement(xi, T, E) for xi in x], color=BLUE_RAMP[2 + j - (j > 2)],
               lw=2, label=f"{T} K")
    a.set_yscale("log")
    a.set_xlabel("Ce / (Al+Ce) in the surface layer (%)")
    a.set_ylabel(r"[V$_O$] relative to undoped γ-Al$_2$O$_3$")
    a.set_title("(c) Equilibrium surface V$_O$ enhancement", loc="left")
    a.legend(title="T", frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "fig3_Ce_doping_O_vacancy.png"), dpi=220)
    summ = dict(E_vf_by_k={int(k): float(v) for k, v in E.items()},
                enh_973K={f"{p}%": float(vac_enhancement(p / 100, 973, E)) for p in (0, 5, 10, 20, 30)})
    json.dump(summ, open(os.path.join(ROOT, "data", "ce_summary.json"), "w"), indent=1)
    print(summ)


if __name__ == "__main__":
    main()
