"""Merge the two small NiCu particles' compositions and build a large particle at that Cu/Ni ratio.

Cu fills the lowest-coordinated exposed sites first (corner -> edge -> (100) -> (111));
within a facet class, the rows nearest the edges are filled first on every facet, so the
overflow spreads inward from all edges evenly. No energetics, just the site-filling rule used for the renders.
"""
import os

import numpy as np
from ase.io import write
from scipy.spatial import cKDTree

from build_and_render import RENDER, STRUCT, render
from model_lib import (A_NI, classify_sites, decorate_cu, dope_ce_and_vacancies,
                       gamma_alumina_slab, ni_particle, place_on)

ORDER = ["corner", "edge", "face100", "face111"]


def composition(L, c, cut):
    p = ni_particle(L, c, cut)
    _, t = classify_sites(p)
    p = decorate_cu(p, t)
    return int((p.symbols == "Cu").sum()), int((p.symbols == "Ni").sum())


def fill_cu(p, types, n_cu, seed=0):
    """Fill Cu by class; inside a facet class, rows closest to the edges go first, on every facet."""
    p = p.copy()
    rng = np.random.default_rng(seed)
    low = np.where(np.isin(types, ORDER[:2]))[0]
    dist_edge = cKDTree(p.positions[low]).query(p.positions)[0]
    picked = []
    for typ in ORDER:
        idx = np.where(types == typ)[0]
        key = np.round(dist_edge[idx], 1) + rng.random(len(idx)) * 1e-3  # ring by ring, random ties
        picked += list(idx[np.argsort(key)])
    # surface full -> remaining Cu goes in shell by shell, nearest to the surface first
    rest = np.setdiff1d(np.arange(len(p)), picked)
    depth = cKDTree(p.positions[picked]).query(p.positions[rest])[0]
    picked += list(rest[np.argsort(depth, kind="stable")])
    p.symbols[picked[:n_cu]] = "Cu"
    return p


def main():
    small = composition(5, 1, 0)
    medium = composition(9, 3, 1)
    cu, ni = small[0] + medium[0], small[1] + medium[1]
    ratio = cu / ni
    print(f"small Cu{small[0]}Ni{small[1]} + medium Cu{medium[0]}Ni{medium[1]} "
          f"-> Cu{cu}Ni{ni}, Cu/Ni = {ratio:.3f}")

    big = ni_particle(13, 4, cut_bottom_layers=2)
    _, t = classify_sites(big)
    n_cu = int(round(len(big) * ratio / (1 + ratio)))
    big = fill_cu(big, t, n_cu)
    on = {k: int(((big.symbols == "Cu") & (t == k)).sum()) for k in ORDER}
    print(f"large: {len(big)} atoms -> Cu{n_cu}Ni{len(big) - n_cu}; Cu by site: {on}; "
          f"corner+edge sites available: {int(np.isin(t, ORDER[:2]).sum())}")

    slab = gamma_alumina_slab(6, 6, n_mixed=3, seed=0)
    centre = slab.cell[0, :2] / 2 + slab.cell[1, :2] / 2
    slab, vac, _ = dope_ce_and_vacancies(slab, 16, 0.6, centre, exclude_r=16, seed=4)
    m = place_on(slab, big)
    write(os.path.join(STRUCT, "size_large_mergedratio.xyz"), m)
    render(m, "05_large_merged_ratio", "35z,-58x", width=1200, vac_pos=vac)
    render(big, "05_large_merged_ratio_particle", "35z,-60x", width=1200)
    return dict(small=small, medium=medium, cu=cu, ni=ni, ratio=ratio, n_cu=n_cu,
                n=len(big), on=on)


if __name__ == "__main__":
    main()
