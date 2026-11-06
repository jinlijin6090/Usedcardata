#!/usr/bin/env python
"""Sanity checks of the ReaxFF parameter file against experiment / typical GGA-DFT values.

Run BEFORE production (especially before the C+O co-deposition task):
    python validate/ff_checks.py [--ffield ff/ffield.reax.FC] [--out validate/ff_report.md]

Reference values are approximate (experiment where available, otherwise the usual range of
GGA-DFT literature values). Replace the 'ref' entries with your own DFT numbers once you have them.
"""
import argparse
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "md"))
from common import write_data, TYPES

KCAL2EV = 1 / 23.060548


class Calc:
    def __init__(self, ffield, threads):
        self.ff = ffield
        self.threads = threads
        self.tmp = tempfile.mkdtemp(prefix="ffcheck_")

    def energy(self, pos, types, box, fixed=None, relax=True, box_relax=None):
        """Return (E_eV, final positions, final box)."""
        from lammps import lammps
        args = ["-log", "none", "-screen", "none", "-nocite"]
        if self.threads > 1:
            args += ["-sf", "omp", "-pk", "omp", str(self.threads)]
        L = lammps(cmdargs=args)
        f = os.path.join(self.tmp, "s.data")
        write_data(f, np.asarray(pos, float), np.asarray(types, int), np.asarray(box, float))
        for c in ["units real", "atom_style charge", "boundary p p p", "atom_modify map yes",
                  f"read_data {f}", "pair_style reaxff NULL checkqeq yes safezone 3.0 mincap 200",
                  f"pair_coeff * * {self.ff} C Ni O",
                  "fix qeq all qeq/reaxff 1 0.0 10.0 1.0e-8 reaxff", "neighbor 2.0 bin",
                  "neigh_modify every 1 delay 0 check yes", "thermo_style custom step pe"]:
            L.command(c)
        if fixed is not None and np.any(fixed):
            ids = " ".join(str(i + 1) for i in np.where(fixed)[0])
            L.command(f"group fx id {ids}")
            L.command("fix fz fx setforce 0.0 0.0 0.0")
        if box_relax:
            L.command(f"fix br all box/relax {box_relax} 0.0 vmax 0.002")
        if relax:
            L.command("min_style cg")
            L.command("minimize 1.0e-12 1.0e-8 20000 200000")
            if box_relax:
                L.command("unfix br")
                L.command("minimize 1.0e-12 1.0e-8 20000 200000")
        else:
            L.command("run 0")
        e = L.get_thermo("pe") * KCAL2EV
        n = L.get_natoms()
        x = np.array(L.numpy.extract_atom("x")[:n], float)
        idx = np.argsort(np.array(L.numpy.extract_atom("id")[:n]))
        lo, hi = L.extract_box()[:2]
        L.close()
        return e, x[idx], np.array(hi) - np.array(lo)


def fcc_bulk(a, n=4):
    base = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return ((g[:, None, :] + base[None]).reshape(-1, 3)) * a, np.full(3, n * a)


def rect_ni111(a, nx=3, nyr=2, nl=4, vac=20.0):
    """Orthogonal Ni(111) slab: nx x (2*nyr) atoms per layer, ABC stacking along z."""
    d = a / np.sqrt(2)
    h = a / np.sqrt(3)
    Lx, Ly = nx * d, nyr * d * np.sqrt(3)
    base = np.array([[0, 0], [d / 2, d * np.sqrt(3) / 2]])
    stack = [np.zeros(2), np.array([d / 2, d * np.sqrt(3) / 6]), np.array([d, d * np.sqrt(3) / 3])]
    pos = []
    for l in range(nl):
        for i in range(nx):
            for j in range(nyr):
                for b in base:
                    xy = b + np.array([i * d, j * d * np.sqrt(3)]) + stack[l % 3]
                    pos.append([xy[0] % Lx, xy[1] % Ly, 5.0 + l * h])
    pos = np.array(pos)
    box = np.array([Lx, Ly, 5.0 + (nl - 1) * h + vac])
    top_z = 5.0 + (nl - 1) * h
    fcc_xy = (np.array([d, d * np.sqrt(3) / 2]) + stack[nl % 3]) % [Lx, Ly]      # above layer nl-3
    hcp_xy = (np.array([d, d * np.sqrt(3) / 2]) + stack[(nl - 2) % 3]) % [Lx, Ly]  # above layer nl-2
    top_xy = (np.array([d, d * np.sqrt(3) / 2]) + stack[(nl - 1) % 3]) % [Lx, Ly]
    fixed = pos[:, 2] < 5.0 + 1.5 * h
    return pos, box, top_z, dict(fcc=fcc_xy, hcp=hcp_xy, top=top_xy), fixed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ffield", default=os.path.join(ROOT, "ff", "ffield.reax.FC"))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("OMP_NUM_THREADS", "1")))
    ap.add_argument("--out", default=os.path.join(HERE, "ff_report.md"))
    a = ap.parse_args()
    C = Calc(a.ffield, a.threads)
    Cc, Ni, O = TYPES["C"], TYPES["Ni"], TYPES["O"]
    rows = []

    def add(name, val, ref, unit, note=""):
        rows.append((name, val, ref, unit, note))
        print(f"{name:45s} {val:10.3f}   ref {ref:>14s} {unit}  {note}", flush=True)

    big = np.full(3, 30.0)
    e_at = {e: C.energy([[15, 15, 15]], [TYPES[e]], big, relax=False)[0] for e in ("C", "Ni", "O")}

    # --- Ni bulk
    pos, box = fcc_bulk(3.524)
    e, _, bx = C.energy(pos, np.full(len(pos), Ni), box, box_relax="iso")
    a0 = bx[0] / 4
    add("Ni fcc lattice constant", a0, "3.524 (exp)", "A")
    add("Ni cohesive energy", -(e / len(pos) - e_at["Ni"]), "4.44 (exp)", "eV/atom")

    # --- graphene
    d0 = 1.42
    gp, gt = [], []
    for i in range(6):
        for j in range(4):
            for fx, fy in [(0, 0), (0.5, 1 / 6), (0.5, 0.5), (0, 2 / 3)]:
                gp.append([(i + fx) * d0 * np.sqrt(3), (j + fy) * 3 * d0, 10.0])
    gbox = np.array([6 * d0 * np.sqrt(3), 4 * 3 * d0, 20.0])
    e, _, bx = C.energy(gp, np.full(len(gp), Cc), gbox, box_relax="x 0.0 y")
    add("graphene C-C bond", bx[0] / 6 / np.sqrt(3), "1.42 (exp)", "A")
    add("graphene cohesive energy", -(e / len(gp) - e_at["C"]), "7.9 (exp, ~7.4-7.9)", "eV/atom")

    # --- gas molecules
    def mol(xyz, ty):
        e, x, _ = C.energy(np.array(xyz) + 15, ty, big)
        return e, x
    eCO, x = mol([[0, 0, 0], [0, 0, 1.13]], [Cc, O])
    add("CO bond length", np.linalg.norm(x[1] - x[0]), "1.128 (exp)", "A")
    add("CO dissociation energy", e_at["C"] + e_at["O"] - eCO, "11.1 (exp)", "eV")
    eO2, x = mol([[0, 0, 0], [0, 0, 1.21]], [O, O])
    add("O2 bond length", np.linalg.norm(x[1] - x[0]), "1.207 (exp)", "A")
    add("O2 dissociation energy", 2 * e_at["O"] - eO2, "5.1-5.2 (exp)", "eV")
    eCO2, x = mol([[0, 0, 0], [0, 0, 1.16], [0, 0, -1.16]], [Cc, O, O])
    add("CO2 C=O bond length", np.linalg.norm(x[1] - x[0]), "1.160 (exp)", "A")
    add("CO2 atomization energy", e_at["C"] + 2 * e_at["O"] - eCO2, "16.6 (exp)", "eV")
    add("CO2 -> CO + O(atom) energy", eCO + e_at["O"] - eCO2, "5.5 (exp)", "eV")

    # --- Ni(111) adsorption
    spos, sbox, ztop, sites, fixed = rect_ni111(a0)
    st = np.full(len(spos), Ni)
    e_slab, _, _ = C.energy(spos, st, sbox, fixed=fixed)

    def ads(species, site, heights):
        p = np.vstack([spos] + [[*sites[site], ztop + h] for h in heights])
        t = np.concatenate([st, species])
        fx = np.concatenate([fixed, np.zeros(len(species), bool)])
        return C.energy(p, t, sbox, fixed=fx)[0]
    for el, ref in (("C", "-6.5 to -7.0 (GGA)"), ("O", "-5.0 to -5.6 (GGA)")):
        for site in ("fcc", "hcp"):
            e = ads([TYPES[el]], site, [1.1])
            add(f"{el} on Ni(111) {site} (vs gas atom)", e - e_slab - e_at[el], ref, "eV")
    for site in ("top", "fcc"):
        e = ads([Cc, O], site, [1.8 if site == "top" else 1.3, 2.95 if site == "top" else 2.45])
        add(f"CO on Ni(111) {site} (vs gas CO)", e - e_slab - eCO,
            "-1.3 (exp) / -1.5 to -1.9 (GGA)", "eV")
    eC = ads([Cc], "hcp", [1.1]) - e_slab
    eOads = ads([O], "fcc", [1.1]) - e_slab
    eCOads = ads([Cc, O], "fcc", [1.3, 2.45]) - e_slab
    add("C* + O* -> CO* (surface reaction energy)", eCOads - eC - eOads, "~ -1 to +0.5 (GGA, coverage dep.)", "eV")
    add("C* + O* -> CO(g)", eCO - eC - eOads, "~ +0.3 to +1.0 (GGA)", "eV")

    with open(a.out, "w") as f:
        f.write(f"# ReaxFF sanity checks\n\nforce field: `{os.path.relpath(a.ffield, ROOT)}`\n\n")
        f.write("| quantity | ReaxFF | reference | unit |\n|---|---|---|---|\n")
        for n, v, r, u, _ in rows:
            f.write(f"| {n} | {v:.3f} | {r} | {u} |\n")
        f.write("\nReference values are approximate; replace with your own DFT numbers (same "
                "functional for all entries) before drawing quantitative conclusions.\n")
    print("\nwritten", a.out)


if __name__ == "__main__":
    main()
