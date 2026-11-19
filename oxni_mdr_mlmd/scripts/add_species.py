"""Add CH4/CO2 gas or C*/O*/H* adatoms to a (typically annealed) structure.

  python add_species.py ../runs/S1_anneal_final.extxyz --gas 12 12 --height 24 --out ../structures/S1_annealed_gas.extxyz
  python add_species.py ../runs/S1_anneal_final.extxyz --adatoms C=4,O=4 --radius 5 --out ../structures/S2_annealed_CO.extxyz
Pore centre/radius are read from atoms.info (written by build_models.py and kept by run_md.py).
"""
import argparse
import os
import sys

import numpy as np
from ase.io import read, write

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import meo_lib as ml  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("structure")
ap.add_argument("--gas", nargs=2, type=int, metavar=("N_CH4", "N_CO2"))
ap.add_argument("--height", type=float, default=24.0)
ap.add_argument("--adatoms", default=None, help="e.g. C=4,O=4")
ap.add_argument("--radius", type=float, default=None, help="disc radius for adatoms (default: pore radius - 1)")
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--out", required=True)
a = ap.parse_args()

at = read(a.structure)
for k in ("momenta",):
    if k in at.arrays:
        del at.arrays[k]
if a.gas:
    at.info.pop("wall_z", None)
    at = ml.add_gas(at, a.gas[0], a.gas[1], height=a.height, seed=a.seed)
if a.adatoms:
    counts = {k: int(v) for k, v in (x.split("=") for x in a.adatoms.split(","))}
    c = at.info.get("pore_center_xy")
    r = a.radius or (float(at.info["pore_radius"]) - 1.0 if "pore_radius" in at.info else None)
    at = ml.add_adatoms(at, counts, None if c is None else np.asarray(c), r, seed=a.seed)
    at.info["wall_z"] = float(at.positions[:, 2].max() + 10.0)
    at.cell[2, 2] = max(at.cell[2, 2], at.info["wall_z"] + 4.0)
write(a.out, at)
print(ml.summary(at))
