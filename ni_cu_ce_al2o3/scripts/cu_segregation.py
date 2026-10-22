"""Cu occupation of low-coordinated Ni sites vs. particle size and Cu/Ni ratio.

1. EMT (Jacobsen) segregation energies of a single Cu in a Ni truncated octahedron,
   E_seg(type) = E(Cu at site) - E(Cu in core), with local relaxation.
2. Site statistics of supported truncated-octahedral Ni particles, d = 1-10 nm.
3. Site-occupancy (Fermi-Dirac / McLean) model: for a given total Cu/Ni ratio the Cu
   chemical potential mu is solved so that sum_i N_i f(E_i - mu) = N_Cu.
   Coverage of corner+edge sites is then read off. Cu-Cu interactions are neglected.
"""
import json
import os

import numpy as np
from ase.calculators.emt import EMT
from ase.optimize import FIRE
from scipy.optimize import brentq
from scipy.spatial import cKDTree

from model_lib import A_NI, classify_sites, ni_particle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
os.makedirs(DATA, exist_ok=True)
KB = 8.617e-5
TYPES = ["corner", "edge", "face100", "face111", "subsurface", "interface", "bulk"]


def emt_segregation(length=10, cutoff=3, n_per_type=4, seed=0):
    """Mean EMT segregation energy per site class (several symmetry-inequivalent picks)."""
    from ase.constraints import FixAtoms
    rng = np.random.default_rng(seed)
    p = ni_particle(length, cutoff, cut_bottom_layers=0)
    cn, t = classify_sites(p, support_contact=False)
    com = p.positions.mean(axis=0)
    r = np.linalg.norm(p.positions - com, axis=1)
    surf = cn < 12
    near_surf = np.array([np.any(np.linalg.norm(p.positions[surf] - x, axis=1) < 2.7) for x in p.positions])
    t = t.astype(object)
    t[(cn == 12) & near_surf] = "subsurface"

    def energy(i):
        a = p.copy()
        a.symbols[i] = "Cu"
        a.calc = EMT()
        d = np.linalg.norm(a.positions - a.positions[i], axis=1)
        a.set_constraint(FixAtoms(indices=np.where(d > 6.0)[0]))
        FIRE(a, logfile=None).run(fmax=0.02, steps=500)
        return a.get_potential_energy()

    e_core = energy(int(np.argmin(r)))
    seg, spread = {}, {}
    for typ in ["corner", "edge", "face100", "face111", "subsurface"]:
        idx = np.where(t == typ)[0]
        idx = rng.choice(idx, min(n_per_type, len(idx)), replace=False)
        vals = [energy(int(i)) - e_core for i in idx]
        seg[typ], spread[typ] = float(np.mean(vals)), float(np.std(vals))
    return seg, spread, len(p)


def site_counts(c, cut=None):
    """Supported regular truncated octahedron with edge parameter c (length = 3c+1)."""
    L = 3 * c + 1
    if cut is None:
        cut = max(1, c // 2)  # flatten particles a bit: ~ wetting of the oxide
    p = ni_particle(L, c, cut_bottom_layers=cut)
    cn, t = classify_sites(p)
    t = t.astype(object)
    tree = cKDTree(p.positions)
    nn = A_NI / np.sqrt(2)
    exposed = np.isin(t, ["corner", "edge", "face100", "face111"])
    # first sub-surface shell: fully coordinated atoms bonded to an exposed atom
    nb = tree.query_ball_point(p.positions, 1.2 * nn)
    sub = np.array([t[i] == "bulk" and exposed[nb[i]].any() for i in range(len(p))])
    t[sub] = "subsurface"
    span = np.ptp(p.positions[:, 0]) + A_NI / np.sqrt(2)
    v_at = A_NI ** 3 / 4
    d_eq = 2 * (3 * len(p) * v_at / (4 * np.pi)) ** (1 / 3)
    counts = {k: int((t == k).sum()) for k in TYPES}
    # metal/oxide perimeter = bottom-layer atoms with < 4 in-plane neighbours
    z = p.positions[:, 2]
    bot = np.where(t == "interface")[0]
    inplane = np.array([sum(abs(z[j] - z[i]) < 0.3 for j in nb[i]) - 1 for i in bot])
    perim = int((inplane < 4).sum())
    return dict(c=c, N=len(p), d_nm=d_eq / 10, span_nm=span / 10, perimeter=perim, **counts)


def occupancy(counts, e_seg, cu_ni, T):
    """Fermi-Dirac occupation of each site class for a Cu/Ni ratio."""
    N = counts["N"]
    n_cu = N * cu_ni / (1 + cu_ni)
    E = {k: e_seg.get(k, 0.0) for k in TYPES}
    E["bulk"] = 0.0
    E["interface"] = 0.0  # oxide-bonded bottom layer: no surface driving force assumed
    kT = KB * T

    def f(mu):
        return {k: 1 / (1 + np.exp(np.clip((E[k] - mu) / kT, -200, 200))) for k in TYPES}

    def resid(mu):
        fk = f(mu)
        return sum(counts[k] * fk[k] for k in TYPES) - n_cu

    mu = brentq(resid, -5, 5)
    return f(mu)


def main():
    seg, spread, n_ref = emt_segregation()
    print("EMT segregation energies (eV), Ni%d host:" % n_ref,
          {k: f"{seg[k]:.3f}+-{spread[k]:.3f}" for k in seg})
    sizes = [site_counts(c) for c in range(2, 15)]
    for s in sizes:
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in s.items()})
    json.dump(dict(e_seg_emt=seg, e_seg_emt_std=spread, sizes=sizes), open(os.path.join(DATA, "cu_segregation.json"), "w"),
              indent=1)


if __name__ == "__main__":
    main()
