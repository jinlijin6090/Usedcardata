"""Collect VASP single points / AIMD into an extxyz training set (REF_energy, REF_forces).

  python vasp_to_extxyz.py runs/*/vasprun.xml --out labelled.extxyz [--every 10]
REF_energy is the force-consistent free energy (TOTEN), which matches the forces.
Original tags are re-attached if a sibling 'source.extxyz' exists in the run folder.
"""
import argparse
import os

from ase.io import read, write

ap = argparse.ArgumentParser()
ap.add_argument("files", nargs="+")
ap.add_argument("--every", type=int, default=1)
ap.add_argument("--out", required=True)
a = ap.parse_args()
out = []
for f in a.files:
    frames = read(f, index="::%d" % a.every)
    src = os.path.join(os.path.dirname(f), "source.extxyz")
    tags = read(src).get_tags() if os.path.exists(src) else None
    for fr in frames:
        e = fr.calc.results.get("free_energy", fr.get_potential_energy())
        F = fr.get_forces()
        fr.calc = None
        fr.info["REF_energy"] = e
        fr.arrays["REF_forces"] = F
        fr.info["config_type"] = os.path.basename(os.path.dirname(os.path.abspath(f)))
        if tags is not None:
            fr.set_tags(tags)
        out.append(fr)
write(a.out, out)
print("wrote", len(out), "frames to", a.out)
