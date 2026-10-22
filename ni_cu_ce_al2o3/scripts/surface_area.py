"""Specific surface area and exposed-Ni statistics of supported NiCu particles (atom counting).

Exposed atoms = corner, edge, (100), (111) sites (the bottom layer touching the oxide is not
counted as exposed). Area per exposed atom: Ni 6.49 A^2 (1.54e19 at/m^2), Cu 6.85 A^2
(1.47e19 at/m^2) - the usual chemisorption site densities.
Cu is placed with the same rule as the renders: corners -> edges -> (100) -> (111), facet rows
nearest the edges first.

Scenarios: pure Ni; Cu exactly on corners+edges (size-matched); fixed Cu/Ni = 0.10; fixed
Cu/Ni = 0.27 (the pooled ratio of the 1.3 + 2.0 nm particles).
"""
import csv
import os

import numpy as np

from merge_ratio import fill_cu
from model_lib import classify_sites, ni_particle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NA = 6.02214e23
M = {"Ni": 58.693, "Cu": 63.546}
AREA = {"Ni": 6.49e-20, "Cu": 6.85e-20}  # m^2 per surface atom
EXPOSED = ["corner", "edge", "face100", "face111"]
A_ATOM_VOL = 3.524 ** 3 / 4  # A^3 per Ni atom


def stats(p, t):
    sym = np.array(p.get_chemical_symbols())
    exp = np.isin(t, EXPOSED)
    mass_g = sum(M[s] for s in sym) / NA
    n_ni = (sym == "Ni").sum()
    n_ni_exp = (exp & (sym == "Ni")).sum()
    n_cu_exp = (exp & (sym == "Cu")).sum()
    area = n_ni_exp * AREA["Ni"] + n_cu_exp * AREA["Cu"]
    d_eq = 2 * (3 * len(p) * A_ATOM_VOL / (4 * np.pi)) ** (1 / 3) / 10
    return dict(
        d_nm=round(d_eq, 3), N=len(p), N_Ni=int(n_ni), N_Cu=int(len(p) - n_ni),
        CuNi=round((len(p) - n_ni) / n_ni, 4),
        N_exposed=int(exp.sum()), N_Ni_exposed=int(n_ni_exp), N_Cu_exposed=int(n_cu_exp),
        dispersion=round(exp.sum() / len(p), 4),                 # exposed / all metal atoms
        surface_Ni_frac=round(n_ni_exp / exp.sum(), 4),          # Ni share of the exposed sites
        Ni_exposed_frac=round(n_ni_exp / n_ni, 4),               # fraction of Ni atoms exposed
        SSA_metal_m2g=round(area / mass_g, 2),                   # total metal area per g metal
        SSA_Ni_m2g=round(n_ni_exp * AREA["Ni"] / mass_g, 2),    # exposed Ni area per g metal
    )


def particle(L, c, cut):
    p = ni_particle(L, c, cut)
    _, t = classify_sites(p)
    return p, t


def scenarios(p, t):
    n_low = int(np.isin(t, ["corner", "edge"]).sum())
    out = {"pure Ni": p.copy(), "Cu on corners+edges": fill_cu(p, t, n_low)}
    for r in (0.10, 0.27):
        out[f"Cu/Ni = {r:.2f}"] = fill_cu(p, t, int(round(len(p) * r / (1 + r))))
    return out


def main():
    rows = []
    for c in range(2, 15):
        p, t = particle(3 * c + 1, c, max(1, c // 2))
        for name, q in scenarios(p, t).items():
            rows.append(dict(series="size scan", scenario=name, **stats(q, t)))
    # the four particles that were rendered
    for label, (L, c, cut) in {"1.3 nm render": (5, 1, 0), "2.0 nm render": (9, 3, 1),
                               "3.1 nm render": (13, 4, 2)}.items():
        p, t = particle(L, c, cut)
        n_low = int(np.isin(t, ["corner", "edge"]).sum())
        rows.append(dict(series=label, scenario="Cu on corners+edges", **stats(fill_cu(p, t, n_low), t)))
        rows.append(dict(series=label, scenario="pure Ni", **stats(p, t)))
    p, t = particle(13, 4, 2)
    q = fill_cu(p, t, int(round(len(p) * (100 / 368) / (1 + 100 / 368))))
    rows.append(dict(series="3.1 nm merged-ratio render", scenario="Cu/Ni = 0.27", **stats(q, t)))

    out = os.path.join(ROOT, "data", "surface_area.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        if r["series"] != "size scan" or r["scenario"] in ("pure Ni", "Cu on corners+edges"):
            print({k: r[k] for k in ("series", "scenario", "d_nm", "N", "CuNi", "dispersion",
                                     "surface_Ni_frac", "Ni_exposed_frac", "SSA_metal_m2g",
                                     "SSA_Ni_m2g")})
    return rows


if __name__ == "__main__":
    main()
