"""Phase-1 benchmark of universal ML potentials on dry-reforming energetics.

  python benchmark_umlip.py --calcs mace:small mace:medium mace:medium-omat-0 uma:uma-s-1p1/oc20 \
      [--dft-ref ../data/dft_ref.json] [--parity ../dft/labelled.extxyz] --out ../data/bench.json

1) Adsorption energies on Ni(111) p(3x3), 4 layers (bottom 2 frozen), best of a few sites,
   relative to CH4, H2, H2O, CO, CO2 gas references:
     H* : E - E_slab - 1/2 H2          C* : E - E_slab - (CH4 - 2 H2)
     O* : E - E_slab - (H2O - H2)      CH*, CH2*, CH3* : ... - (CH4 - (4-n)/2 H2)
     CO*: E - E_slab - CO              OH*: E - E_slab - (H2O - 1/2 H2)
2) Reaction energies (separated adsorbates):
     CH4(g) -> CH3* + H*,  CO2(g) -> CO* + O*,  CO* -> C* + O*,  C* + O* -> CO(g)
3) Optional parity on your own DFT-labelled frames (extxyz with REF_energy / REF_forces):
   energy MAE per atom and force MAE / RMSE (overall and per element).
Supply DFT numbers (same functional you will fine-tune on) in --dft-ref to get errors.
Barriers (NEB on CH4 -> CH3 + H, C + O -> CO, C + C -> C2) are the next step; see PLAN.md.
"""
import argparse
import json
import os
import sys

import numpy as np
from ase import Atoms
from ase.build import add_adsorbate, fcc111
from ase.constraints import FixAtoms
from ase.io import read
from ase.optimize import FIRE

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calculators import get_calculator  # noqa: E402

FMAX, STEPS = 0.05, 400

ADS = {  # name: (atoms builder, list of (site, height of first atom))
    "H": (lambda: Atoms("H"), [("fcc", 0.9), ("hcp", 0.9)]),
    "C": (lambda: Atoms("C"), [("hcp", 1.0), ("fcc", 1.0)]),
    "O": (lambda: Atoms("O"), [("fcc", 1.1), ("hcp", 1.1)]),
    "OH": (lambda: Atoms("OH", positions=[(0, 0, 0), (0, 0, 0.97)]), [("fcc", 1.3)]),
    "CH": (lambda: Atoms("CH", positions=[(0, 0, 0), (0, 0, 1.10)]), [("hcp", 1.1), ("fcc", 1.1)]),
    "CH2": (lambda: Atoms("CH2", positions=[(0, 0, 0), (0.9, 0, 0.6), (-0.9, 0, 0.6)]),
            [("fcc", 1.2), ("hcp", 1.2)]),
    "CH3": (lambda: Atoms("CH3", positions=[(0, 0, 0), (1.03, 0, 0.36), (-0.51, 0.89, 0.36),
                                            (-0.51, -0.89, 0.36)]), [("fcc", 1.5), ("hcp", 1.5)]),
    "CO": (lambda: Atoms("CO", positions=[(0, 0, 0), (0, 0, 1.17)]), [("hcp", 1.3), ("fcc", 1.3)]),
}


def relax(at, calc, fix=None):
    at = at.copy()
    at.calc = calc
    if fix is not None:
        at.set_constraint(FixAtoms(indices=fix))
    FIRE(at, logfile=None).run(fmax=FMAX, steps=STEPS)
    return at.get_potential_energy()


def gas(calc):
    box = [14.0, 14.5, 15.0]
    mols = {
        "H2": Atoms("H2", positions=[(0, 0, 0), (0, 0, 0.74)]),
        "H2O": Atoms("H2O", positions=[(0, 0, 0), (0.76, 0, 0.59), (-0.76, 0, 0.59)]),
        "CH4": Atoms("CH4", positions=np.vstack([[0, 0, 0], 1.09 / np.sqrt(3) * np.array(
            [[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]])])),
        "CO": Atoms("CO", positions=[(0, 0, 0), (0, 0, 1.13)]),
        "CO2": Atoms("CO2", positions=[(0, 0, 0), (0, 0, 1.16), (0, 0, -1.16)]),
    }
    out = {}
    for k, m in mols.items():
        m.cell = box
        m.center()
        m.pbc = True
        out[k] = relax(m, calc)
    return out


def adsorption(calc):
    slab = fcc111("Ni", (3, 3, 4), a=3.524, vacuum=10.0)
    fix = [i for i, t in enumerate(slab.get_tags()) if t >= 3]
    e_slab = relax(slab, calc, fix)
    E = {}
    for name, (mk, sites) in ADS.items():
        best = np.inf
        for site, h in sites:
            s = slab.copy()
            add_adsorbate(s, mk(), h, site)
            best = min(best, relax(s, calc, fix))
        E[name] = best - e_slab
    return E, e_slab


def summarise(Eads_raw, g):
    """Convert raw E(ads+slab)-E(slab) to adsorption / reaction energies."""
    r = Eads_raw
    ads = {
        "H*": r["H"] - 0.5 * g["H2"],
        "C*": r["C"] - (g["CH4"] - 2 * g["H2"]),
        "O*": r["O"] - (g["H2O"] - g["H2"]),
        "OH*": r["OH"] - (g["H2O"] - 0.5 * g["H2"]),
        "CH*": r["CH"] - (g["CH4"] - 1.5 * g["H2"]),
        "CH2*": r["CH2"] - (g["CH4"] - g["H2"]),
        "CH3*": r["CH3"] - (g["CH4"] - 0.5 * g["H2"]),
        "CO*": r["CO"] - g["CO"],
    }
    rxn = {
        "CH4(g)->CH3*+H*": r["CH3"] + r["H"] - g["CH4"],
        "CH3*->CH2*+H*": r["CH2"] + r["H"] - r["CH3"],
        "CH2*->CH*+H*": r["CH"] + r["H"] - r["CH2"],
        "CH*->C*+H*": r["C"] + r["H"] - r["CH"],
        "CO2(g)->CO*+O*": r["CO"] + r["O"] - g["CO2"],
        "CO*->C*+O*": r["C"] + r["O"] - r["CO"],
        "C*+O*->CO(g)": g["CO"] - r["C"] - r["O"],
        "CO2(g)+H2(g)->CO(g)+H2O(g)": g["CO"] + g["H2O"] - g["CO2"] - g["H2"],
        "CH4(g)+CO2(g)->2CO(g)+2H2(g)": 2 * g["CO"] + 2 * g["H2"] - g["CH4"] - g["CO2"],
    }
    return ads, rxn


def parity(calc, path):
    frames = read(path, ":")
    de, df, els = [], [], []
    for f in frames:
        x = f.copy()
        x.calc = calc
        de.append((x.get_potential_energy() - f.info["REF_energy"]) / len(f))
        df.append(x.get_forces() - f.arrays["REF_forces"])
        els += list(f.get_chemical_symbols())
    df = np.vstack(df)
    els = np.array(els)
    per_el = {e: float(np.abs(df[els == e]).mean()) for e in np.unique(els)}
    return dict(E_MAE_meV_atom=1000 * float(np.abs(de).mean()),
                F_MAE=float(np.abs(df).mean()), F_RMSE=float(np.sqrt((df ** 2).mean())),
                F_MAE_per_element=per_el, n_frames=len(frames))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calcs", nargs="+", default=["mace:small", "mace:medium"])
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--dtype", default="float64")
    ap.add_argument("--dft-ref", default=None, help="json {ads:{..}, rxn:{..}} from your DFT")
    ap.add_argument("--parity", default=None)
    ap.add_argument("--out", default="../data/benchmark_umlip.json")
    a = ap.parse_args()
    ref = json.load(open(a.dft_ref)) if a.dft_ref else None
    res = json.load(open(a.out)) if os.path.exists(a.out) else {}
    for spec in a.calcs:
        calc = get_calculator(spec, a.device, a.dtype)
        g = gas(calc)
        raw, _ = adsorption(calc)
        ads, rxn = summarise(raw, g)
        entry = dict(ads=ads, rxn=rxn)
        if ref:
            entry["err_ads"] = {k: ads[k] - v for k, v in ref.get("ads", {}).items() if k in ads}
            entry["err_rxn"] = {k: rxn[k] - v for k, v in ref.get("rxn", {}).items() if k in rxn}
        if a.parity:
            entry["parity"] = parity(calc, a.parity)
        res[spec] = entry
        json.dump(res, open(a.out, "w"), indent=1)
        print(spec)
        for k, v in list(ads.items()) + list(rxn.items()):
            print("  %-30s %8.3f eV" % (k, v))


if __name__ == "__main__":
    main()
