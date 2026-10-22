"""Figure 4: particle-size effect on CH4-activation activity and coking (site-descriptor model).

Ingredients (all per particle, from the site counts of fig 2):
* CH4 dissociation barrier per site class (DFT literature range for Ni; Cu-covered sites
  are effectively inactive): rate_k = exp(-Ea_k / kT).
* Cu occupancy of each site class from the EMT segregation model (fig 2).
* Carbon (graphene) nucleation needs free Ni step/corner sites (Bengaard et al., J. Catal.
  209 (2002) 365) and is suppressed on small particles (Gibbs-Thomson destabilisation of
  carbon / no room for a stable graphene island): logistic factor S(d), d_c = 4 nm.
* Carbon gasification by lattice O supplied at the metal/oxide perimeter, boosted by the
  O-vacancy enhancement that Ce gives (fig 3).
Absolute numbers are arbitrary units; the trends are the point.
"""
import json
import os

import matplotlib.pyplot as plt
import numpy as np

from fig_cu import E_EMT, SIZES, occ, required_ratio
from plot_style import SERIES, setup

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KB = 8.617e-5
T = 973.0
EA = dict(corner=0.75, edge=0.80, face100=0.90, face111=1.00)  # eV, CH4 -> CH3 + H on Ni
EA_CU = 1.45  # eV, CH4 activation on Cu-like sites
EA_PERIM = 0.85
D_C, W_C = 4.0, 0.8  # nm, size threshold for stable carbon nucleation/whisker growth
K_GAS = 0.6  # gasification strength per perimeter site (relative to nucleation term)


def ce_enhancement(x_ce):
    p = os.path.join(ROOT, "data", "ce_summary.json")
    if x_ce == 0 or not os.path.exists(p):
        return 1.0
    e = json.load(open(p))["enh_973K"]
    xs = np.array([float(k[:-1]) for k in e]) / 100
    return float(np.interp(x_ce, xs, list(e.values())))


def descriptors(s, cu_ni, x_ce=0.0):
    f = occ(s, E_EMT, cu_ni) if cu_ni > 0 else {k: 0.0 for k in s}
    kT = KB * T
    rate = 0.0
    for k, ea in EA.items():
        n = s[k]
        rate += n * ((1 - f[k]) * np.exp(-ea / kT) + f[k] * np.exp(-EA_CU / kT))
    rate += s["perimeter"] * np.exp(-EA_PERIM / kT)
    n_ni = s["N"] / (1 + cu_ni)  # Ni atoms in a particle that holds N metal atoms total
    activity = rate / s["N"]  # per metal atom (~ per gram)
    free_steps = sum(s[k] * (1 - f[k]) for k in ("corner", "edge"))
    size_f = 1 / (1 + np.exp(-(s["d_nm"] - D_C) / W_C))
    nucleation = rate * free_steps / max(s["corner"] + s["edge"], 1) * size_f
    gasif = K_GAS * s["perimeter"] * np.exp(-EA_PERIM / kT) * ce_enhancement(x_ce) ** 0.5
    coke = nucleation / (nucleation + gasif)
    return activity, coke, n_ni


def main():
    setup()
    d = np.array([s["d_nm"] for s in SIZES])
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.6))
    cases = [("Ni / γ-Al$_2$O$_3$", 0.0, 0.0, "-"),
             ("NiCu (Cu/Ni = 0.10) / γ-Al$_2$O$_3$", 0.10, 0.0, "-"),
             ("NiCu (Cu/Ni = 0.10) / 10 % Ce-γ-Al$_2$O$_3$", 0.10, 0.10, "--"),
             ("NiCu (size-matched Cu, 90 % edge cover) / 10 % Ce", "match", 0.10, ":")]
    act_ref = descriptors(SIZES[4], 0.0)[0]
    for j, (lab, cu, ce, ls) in enumerate(cases):
        A, C = [], []
        for s in SIZES:
            r = required_ratio(s, E_EMT, 0.90) if cu == "match" else cu
            a_, c_, _ = descriptors(s, r, ce)
            A.append(a_ / act_ref)
            C.append(c_)
        ax[0, 0].plot(d, A, ls, color=SERIES[j], marker="o", ms=4, label=lab)
        ax[0, 1].plot(d, C, ls, color=SERIES[j], marker="o", ms=4, label=lab)
    ax[0, 0].set_yscale("log")
    ax[0, 0].set_xlabel("Ni particle diameter (nm)")
    ax[0, 0].set_ylabel("CH$_4$ activation rate per metal atom\n(rel. to Ni, 4.3 nm)")
    ax[0, 0].set_title(f"(a) Mass-specific activity, T = {T:.0f} K", loc="left")
    ax[0, 1].set_xlabel("Ni particle diameter (nm)")
    ax[0, 1].set_ylabel("Coke propensity index (0-1)")
    ax[0, 1].set_title("(b) Coking tendency", loc="left")
    ax[0, 1].set_ylim(0, 1)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.515))

    # (c,d) maps over size x Cu/Ni at 10 % Ce
    ratios = np.logspace(-2.5, 0, 40)
    A = np.zeros((len(ratios), len(SIZES)))
    C = np.zeros_like(A)
    for i, r in enumerate(ratios):
        for j, s in enumerate(SIZES):
            a_, c_, _ = descriptors(s, r, 0.10)
            A[i, j], C[i, j] = a_ / act_ref, c_
    ext = dict(shading="gouraud")
    m = ax[1, 0].pcolormesh(d, ratios, C, cmap="Oranges", vmin=0, vmax=1, **ext)
    fig.colorbar(m, ax=ax[1, 0], label="Coke propensity index")
    cs = ax[1, 0].contour(d, ratios, C, levels=[0.2, 0.5], colors="#0b0b0b", linewidths=1)
    ax[1, 0].clabel(cs, fmt="%.1f", fontsize=8.5)
    ax[1, 0].set_title("(c) Coking map, 10 % Ce-γ-Al$_2$O$_3$", loc="left")
    score = A * (1 - C)
    score_n = score / score.max(axis=0, keepdims=True)
    m = ax[1, 1].pcolormesh(d, ratios, score_n, cmap="Blues", vmin=0, vmax=1, **ext)
    fig.colorbar(m, ax=ax[1, 1], label="Activity × (1 − coke index),\nnormalised at each size")
    best = ratios[np.argmax(score, axis=0)]
    ax[1, 1].plot(d, best, "-", color=SERIES[1], lw=2, label="best Cu/Ni at each size")
    ax[1, 1].legend(frameon=False, fontsize=9, loc="lower left")
    ax[1, 1].set_title("(d) Stable-activity map, 10 % Ce", loc="left")
    r90 = [required_ratio(s, E_EMT, 0.90) for s in SIZES]
    for a in ax[1]:
        a.set_yscale("log")
        a.set_xlabel("Ni particle diameter (nm)")
        a.set_ylabel("Cu/Ni molar ratio")
        a.plot(d, r90, "--", color="#52514e", lw=1.3)
        a.grid(False)
    ax[1, 0].text(d[-1], r90[-1] * 1.15, "90 % edge cover", ha="right", fontsize=8.5, color="#0b0b0b")

    fig.tight_layout(h_pad=5.0)
    fig.savefig(os.path.join(ROOT, "figures", "fig4_size_activity_coking.png"), dpi=220)
    print("best Cu/Ni per size:", dict(zip(np.round(d, 1), np.round(best, 3))))


if __name__ == "__main__":
    main()
