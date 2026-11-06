"""Shared helpers: element table, LAMMPS data/dump I/O, fragment detection."""
import json
import os
import numpy as np
from scipy.spatial import cKDTree

# fixed LAMMPS type mapping used everywhere in this package
TYPES = {"C": 1, "Ni": 2, "O": 3}
ELEM = {v: k for k, v in TYPES.items()}
MASS = {"C": 12.011, "Ni": 58.6934, "O": 15.999}
KB_KCAL = 0.0019872041          # kcal/mol/K
AMU_A2_FS2_TO_KCAL = 2390.0573  # 1 amu*(A/fs)^2 in kcal/mol

# bond cutoffs (A) for fragment / bonding analysis
CUT = {("C", "C"): 1.85, ("C", "O"): 1.60, ("O", "O"): 1.60,
       ("C", "Ni"): 2.30, ("Ni", "O"): 2.30, ("Ni", "Ni"): 3.00}
MAXCUT = max(CUT.values())


def cut(a, b):
    return CUT[tuple(sorted((a, b)))]


def write_json_atomic(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1)
    os.replace(tmp, path)


def write_data(path, pos, types, box, charges=None, vel=None, comment="ni_coke_md"):
    """Write a LAMMPS data file (atom_style charge) with the 3-type C/Ni/O convention."""
    n = len(pos)
    q = np.zeros(n) if charges is None else charges
    with open(path, "w") as f:
        f.write(f"{comment}\n\n{n} atoms\n3 atom types\n\n")
        for k, ax in enumerate("xyz"):
            f.write(f"0.0 {box[k]:.6f} {ax}lo {ax}hi\n")
        f.write("\nMasses\n\n")
        for e, t in TYPES.items():
            f.write(f"{t} {MASS[e]}\n")
        f.write("\nAtoms # charge\n\n")
        for i in range(n):
            f.write(f"{i+1} {types[i]} {q[i]:.6f} {pos[i,0]:.6f} {pos[i,1]:.6f} {pos[i,2]:.6f}\n")
        if vel is not None:
            f.write("\nVelocities\n\n")
            for i in range(n):
                f.write(f"{i+1} {vel[i,0]:.8f} {vel[i,1]:.8f} {vel[i,2]:.8f}\n")


def iter_dump(path):
    """Yield (step, box_lengths, ids, types, pos) for every complete frame of a LAMMPS custom dump
    written as `id type x y z q`."""
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(f.readline())
            f.readline()
            n = int(f.readline())
            f.readline()
            lo_hi = [f.readline().split() for _ in range(3)]
            hdr = f.readline().split()[2:]
            rows = [f.readline() for _ in range(n)]
            if len(rows[-1].split()) != len(hdr):      # truncated last frame
                return
            d = np.loadtxt(rows, ndmin=2)
            box = np.array([float(b[1]) - float(b[0]) for b in lo_hi])
            col = {h: i for i, h in enumerate(hdr)}
            yield (step, box, d[:, col["id"]].astype(int), d[:, col["type"]].astype(int),
                   d[:, [col["x"], col["y"], col["z"]]])


def unwrap_about_ni(pos, types, box):
    """Shift every atom to the periodic image nearest the Ni centre; return (pos, centre)."""
    ni = pos[types == TYPES["Ni"]]
    ref = ni[0]
    ni_u = ref + (ni - ref) - box * np.round((ni - ref) / box)
    c = ni_u.mean(0)
    pos = c + (pos - c) - box * np.round((pos - c) / box)
    return pos, c


def fragments(pos, types):
    """Connected components of the distance-bond graph. pos must already be unwrapped.
    Returns list of index arrays and the neighbour pair list."""
    import networkx as nx
    el = np.array([ELEM[t] for t in types])
    tree = cKDTree(pos)
    pairs = tree.query_pairs(MAXCUT, output_type="ndarray")
    if len(pairs):
        d = np.linalg.norm(pos[pairs[:, 0]] - pos[pairs[:, 1]], axis=1)
        keep = np.array([d[k] <= cut(el[i], el[j]) for k, (i, j) in enumerate(pairs)], dtype=bool)
        pairs = pairs[keep]
    g = nx.Graph()
    g.add_nodes_from(range(len(pos)))
    g.add_edges_from(map(tuple, pairs))
    return [np.array(sorted(c)) for c in nx.connected_components(g)], pairs


def formula(types_subset):
    from collections import Counter
    c = Counter(ELEM[t] for t in types_subset)
    s = ""
    for e in ("C", "O", "Ni"):
        if c[e]:
            s += e + (str(c[e]) if c[e] > 1 else "")
    return s
