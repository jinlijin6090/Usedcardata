#!/usr/bin/env python
"""Build a free-standing fcc Ni cluster as a LAMMPS data file (C/Ni/O 3-type convention).

  truncated octahedron (default, same construction as ase.cluster.Octahedron):
      python build_cluster.py --shape octahedron --length 6 --cutoff 2 -o structures/ni116.data
  roughly spherical cluster with exactly N atoms:
      python build_cluster.py --shape sphere --natoms 201 -o structures/ni201.data
"""
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import write_data, TYPES

A_NI = 3.524


def fcc_points(nmax):
    r = np.arange(-nmax, nmax + 1)
    g = np.array(np.meshgrid(r, r, r, indexing="ij")).reshape(3, -1).T
    return g[(g.sum(1) % 2) == 0]          # fcc sites in units of a/2


def octahedron(length, cutoff):
    """Same atom set as ase.cluster.Octahedron(length, cutoff): bounded by {111} at
    layer length/2 (even) or (length-1)/2 (odd) and {100} at layer length-1-cutoff."""
    if cutoff < 0 or length < 2 * cutoff + 1:
        raise ValueError("need 0 <= cutoff <= (length-1)/2")
    pts = fcc_points(2 * length + 4)                      # integer coords, units of a/2
    centre = np.array([1, 1, 1]) if length % 2 == 0 else np.zeros(3, int)
    l111 = length / 2 if length % 2 == 0 else (length - 1) / 2
    l100 = length - 1 - cutoff
    rel = pts - centre
    s111 = np.max(np.abs(rel @ np.array([[1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1]]).T), axis=1)
    keep = (s111 <= 2 * l111 + 1e-9) & (np.abs(rel).max(1) <= l100 + 1e-9)
    return rel[keep] * A_NI / 2


def sphere(n):
    pts = fcc_points(12) * A_NI / 2
    d = np.linalg.norm(pts, axis=1)
    return pts[np.argsort(d, kind="stable")[:n]]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shape", choices=["octahedron", "sphere"], default="octahedron")
    ap.add_argument("--length", type=int, default=6)
    ap.add_argument("--cutoff", type=int, default=2)
    ap.add_argument("--natoms", type=int, default=116)
    ap.add_argument("--box", type=float, default=50.0)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    pos = octahedron(a.length, a.cutoff) if a.shape == "octahedron" else sphere(a.natoms)
    pos = pos - pos.mean(0) + a.box / 2
    write_data(a.out, pos, np.full(len(pos), TYPES["Ni"]), np.full(3, a.box),
               comment=f"Ni{len(pos)} {a.shape} cluster")
    r = np.linalg.norm(pos - pos.mean(0), axis=1).max()
    print(f"wrote {a.out}: Ni{len(pos)}, radius {r:.2f} A, box {a.box} A")
