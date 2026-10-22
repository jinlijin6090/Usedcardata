"""O-vacancy formation energies on Ce-doped gamma-Al2O3(100) with the MACE-MP-0 potential.

E_vf = E(slab - O) + 1/2 E(O2) - E(slab), surface O sites grouped by the number of Ce
first neighbours. Ce substitutes surface octahedral Al (isovalent Ce_Al in a neutral cell).
Bottom ~3 A of the slab is frozen. Results are appended to data/ce_vacancy_mace.json.
"""
import json
import os
import sys

import numpy as np
import torch
from ase import Atoms
from ase.constraints import FixAtoms
from ase.optimize import FIRE

from model_lib import gamma_alumina_slab

torch.set_num_threads(4)
from mace.calculators import mace_mp  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "ce_vacancy_mace.json")
CALC = mace_mp(model="medium", default_dtype="float32", device="cpu")
FMAX, STEPS = 0.10, 150


def relax(a):
    a = a.copy()
    a.calc = CALC
    a.set_constraint(FixAtoms(indices=np.where(a.positions[:, 2] < 3.0)[0]))
    FIRE(a, logfile=None).run(fmax=FMAX, steps=STEPS)
    return a, a.get_potential_energy()


def o2_energy():
    m = Atoms("O2", positions=[(0, 0, 0), (0, 0, 1.21)], cell=[12, 12, 12])
    m.calc = CALC
    FIRE(m, logfile=None).run(fmax=0.01)
    return m.get_potential_energy()


def pick_ce(slab, n, rng):
    top = slab.positions[:, 2].max()
    cand = np.where((slab.symbols == "Al") & (slab.positions[:, 2] > top - 0.8))[0]
    best = None
    for _ in range(400):  # spread Ce as far apart as possible
        s = rng.choice(cand, n, replace=False)
        dmin = min([slab.get_distance(i, j, mic=True) for i in s for j in s if i < j] or [99])
        if best is None or dmin > best[0]:
            best = (dmin, s)
    return list(best[1])


def main():
    rng = np.random.default_rng(7)
    res = json.load(open(OUT)) if os.path.exists(OUT) else {"runs": []}
    done = {(r["n_ce"], r["o_index"]) for r in res["runs"]}
    if "E_O2" not in res:
        res["E_O2"] = o2_energy()
    half_o2 = 0.5 * res["E_O2"]

    slab0 = gamma_alumina_slab(2, 2, n_mixed=4, vacuum=10.0, seed=0)
    n_top_cat = int(((slab0.symbols == "Al") & (slab0.positions[:, 2] > slab0.positions[:, 2].max() - 0.8)).sum())
    res["n_top_cations"] = n_top_cat
    host, _ = relax(slab0)
    for n_ce in [0, 1, 2, 4]:
        s = host.copy()
        if n_ce:
            s.symbols[pick_ce(s, n_ce, rng)] = "Ce"
        s, e_host = relax(s)
        top = s.positions[s.symbols == "O", 2].max()
        o_top = np.where((s.symbols == "O") & (s.positions[:, 2] > top - 1.2))[0]
        ce = np.where(s.symbols == "Ce")[0]
        k = np.array([sum(s.get_distance(o, c, mic=True) < 2.8 for c in ce) for o in o_top])
        picks = []
        for kk in sorted(set(k)):
            idx = o_top[k == kk]
            picks += [(int(i), int(kk)) for i in rng.choice(idx, min(2 if n_ce else 4, len(idx)), replace=False)]
        for o, kk in picks:
            if (n_ce, o) in done:
                continue
            d = s.copy()
            del d[o]
            d, e_def = relax(d)
            ev = e_def + half_o2 - e_host
            res["runs"].append(dict(n_ce=n_ce, ce_frac_top=n_ce / n_top_cat, o_index=o, k_ce=kk,
                                    E_vf=float(ev)))
            json.dump(res, open(OUT, "w"), indent=1)
            print(f"n_ce={n_ce} O{o} k={kk}  E_vf={ev:.2f} eV", flush=True)


if __name__ == "__main__":
    main()
