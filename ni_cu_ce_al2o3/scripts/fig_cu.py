"""Figure 2: Cu/Ni ratio needed to occupy the corner+edge sites at different particle sizes."""
import json
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from cu_segregation import KB, TYPES
from plot_style import BLUE_RAMP, SERIES, setup

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = json.load(open(os.path.join(ROOT, "data", "cu_segregation.json")))
SIZES = D["sizes"]
T = 973.0  # K, typical reduction / reforming temperature

E_EMT = dict(D["e_seg_emt"], interface=0.0, bulk=0.0)
# bond-cutting alternative: linear in missing bonds, slope from the EMT (111) value
eps = -E_EMT["face111"] / 3
E_BC = dict(corner=-6 * eps, edge=-5 * eps, face100=-4 * eps, face111=-3 * eps,
            subsurface=0.0, interface=0.0, bulk=0.0)
LOW = ("corner", "edge")


def occ(s, E, cu_ni, T=T):
    n_cu = s["N"] * cu_ni / (1 + cu_ni)
    kT = KB * T

    def f(mu):
        return {k: 1 / (1 + np.exp(np.clip((E[k] - mu) / kT, -200, 200))) for k in TYPES}

    mu = brentq(lambda m: sum(s[k] * f(m)[k] for k in TYPES) - n_cu, -10, 10)
    return f(mu)


def theta_low(s, E, cu_ni, T=T):
    f = occ(s, E, cu_ni, T)
    return sum(s[k] * f[k] for k in LOW) / sum(s[k] for k in LOW)


def required_ratio(s, E, target, T=T):
    return brentq(lambda x: theta_low(s, E, x, T) - target, 1e-5, 50)


def main():
    setup()
    d = np.array([s["d_nm"] for s in SIZES])
    N = np.array([s["N"] for s in SIZES])
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.6))

    # (a) site statistics
    a = ax[0, 0]
    frac = {k: np.array([s[k] for s in SIZES]) / N for k in TYPES}
    lab = {"corner": "corner (CN≤6)", "edge": "edge (CN=7)", "face100": "(100) terrace",
           "face111": "(111) terrace"}
    for i, k in enumerate(["corner", "edge", "face100", "face111"]):
        a.plot(d, 100 * frac[k], "-o", color=SERIES[i], ms=4, label=lab[k])
    a.plot(d, 100 * (frac["corner"] + frac["edge"]), "--", color="#52514e", lw=1.5,
           label="corner+edge")
    a.set_xlabel("Ni particle diameter (nm)")
    a.set_ylabel("Fraction of all Ni atoms (%)")
    a.set_yscale("log")
    a.set_title("(a) Site statistics, supported truncated octahedra", loc="left")
    a.legend(frameon=False, fontsize=9)

    # (b) coverage of corner+edge by Cu vs Cu/Ni
    a = ax[0, 1]
    x = np.logspace(-3, 0.3, 160)
    picks = [0, 1, 2, 4, 8, 12]
    for j, i in enumerate(picks):
        s = SIZES[i]
        th = [theta_low(s, E_EMT, xi) for xi in x]
        col = BLUE_RAMP[j]
        a.plot(x, 100 * np.array(th), color=col, lw=2, label=f"{s['d_nm']:.1f} nm")
        geo = (s["corner"] + s["edge"]) / (s["N"] - s["corner"] - s["edge"])
        a.axvline(geo, color=col, lw=0.8, ls=":")
    a.set_xscale("log")
    a.set_xlabel("Overall Cu/Ni molar ratio")
    a.set_ylabel("Corner+edge sites occupied by Cu (%)")
    a.set_title(f"(b) Cu occupancy of low-CN sites, T = {T:.0f} K", loc="left")
    a.legend(title="Ni diameter", frameon=False, fontsize=9, loc="upper left")
    a.text(0.98, 0.04, "dotted: geometric minimum\n(all Cu exactly on corners+edges)",
           transform=a.transAxes, ha="right", va="bottom", fontsize=8.5, color="#52514e")

    # (c) required Cu/Ni vs size
    a = ax[1, 0]
    geo = np.array([(s["corner"] + s["edge"]) / (s["N"] - s["corner"] - s["edge"]) for s in SIZES])
    r90 = np.array([required_ratio(s, E_EMT, 0.90) for s in SIZES])
    r95 = np.array([required_ratio(s, E_EMT, 0.95) for s in SIZES])
    r90bc = np.array([required_ratio(s, E_BC, 0.90) for s in SIZES])
    a.plot(d, geo, "-o", color=SERIES[0], ms=4, label="geometric minimum (T = 0, ideal)")
    a.plot(d, r90, "-s", color=SERIES[1], ms=4, label="90 % coverage, EMT energies")
    a.plot(d, r95, "-^", color=SERIES[2], ms=4, label="95 % coverage, EMT energies")
    a.plot(d, r90bc, "--", color=SERIES[1], lw=1.5, label="90 % coverage, bond-cutting model")
    a.set_yscale("log")
    a.set_xlabel("Ni particle diameter (nm)")
    a.set_ylabel("Required Cu/Ni molar ratio")
    a.set_title("(c) Cu/Ni needed to decorate all corners+edges", loc="left")
    a.legend(frameon=False, fontsize=8.5, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    txt = "Cu/Ni for 90 % coverage (EMT):\n" + "\n".join(
        f"  {d[i]:.1f} nm: {r90[i]:.2f}  ({100*r90[i]/(1+r90[i]):.0f} at.% Cu)" for i in (0, 4, 12))
    a.text(0.03, 0.03, txt, transform=a.transAxes, fontsize=8.5, color="#52514e", va="bottom",
           family="DejaVu Sans Mono")

    # (d) where the Cu goes at the 90 % point, for each size
    a = ax[1, 1]
    order = ["corner", "edge", "face100", "face111", "subsurface", "bulk"]
    names = ["corner", "edge", "(100)", "(111)", "sub-surface", "core/interface"]
    shares = []
    for s, r in zip(SIZES, r90):
        f = occ(s, E_EMT, r)
        n = {k: s[k] * f[k] for k in TYPES}
        n["bulk"] += n["interface"]
        tot = sum(n[k] for k in order)
        shares.append([n[k] / tot for k in order])
    shares = np.array(shares) * 100
    bottom = np.zeros(len(d))
    w = np.gradient(d) * 0.8
    for k in range(len(order)):
        a.bar(d, shares[:, k], width=w, bottom=bottom, color=SERIES[k], edgecolor="white",
              lw=1.0, label=names[k])
        bottom += shares[:, k]
    a.set_xlabel("Ni particle diameter (nm)")
    a.set_ylabel("Share of the Cu atoms (%)")
    a.set_ylim(0, 100)
    a.set_title("(d) Where the Cu sits at 90 % corner+edge coverage", loc="left")
    a.legend(frameon=False, fontsize=8.5, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.13))

    fig.tight_layout()
    out = os.path.join(ROOT, "figures", "fig2_CuNi_ratio_vs_size.png")
    fig.savefig(out, dpi=220)
    rows = [dict(d_nm=round(s["d_nm"], 2), N=s["N"], geo=round(g, 4), r90=round(r, 4),
                 r95=round(q, 4), r90_bondcut=round(b, 4))
            for s, g, r, q, b in zip(SIZES, geo, r90, r95, r90bc)]
    json.dump(rows, open(os.path.join(ROOT, "data", "required_CuNi.json"), "w"), indent=1)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
