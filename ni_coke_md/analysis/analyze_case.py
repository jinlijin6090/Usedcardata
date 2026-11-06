#!/usr/bin/env python
"""Per-case analysis: trajectory + event log -> metrics.csv (one row per analysed frame).

    python analysis/analyze_case.py CASE_DIR [--stride 1] [--force]

Columns
  t_ps, step                 simulation time
  C_ins, O_ins               cumulative atoms inserted (events.csv)
  C_rem, O_rem, C_rem_CO, C_rem_CO2, C_rem_other   cumulative atoms removed as desorbed fragments
  nC, nO                     atoms present
  C_interior / C_surface / C_detached   C inside the Ni convex hull / outside but in the Ni-bound
                             fragment / in a fragment not touching Ni
  O_interior / O_surface / O_detached   same for O
  C1, C2, C3_5, C6p          number of C atoms whose C-C cluster has size 1, 2, 3-5, >=6
  largest                    size of the largest C-C cluster
  r5, r6, r7, rings          rings in the minimum cycle basis of the C-C graph
  sp, sp2, sp3               C atoms with 2, 3, >=4 C neighbours
  CO_bonds                   C-O bonds (any location)
Also writes rings_final.csv (ring-size histogram of the last frame).
"""
import argparse
import csv
import os
import sys
from collections import Counter

import networkx as nx
import numpy as np
from scipy.spatial import Delaunay, cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "md"))
from common import TYPES, iter_dump, unwrap_about_ni, fragments  # noqa: E402

C, NI, O = TYPES["C"], TYPES["Ni"], TYPES["O"]


def read_events(path):
    ev = []
    if not os.path.exists(path):
        return ev
    with open(path) as f:
        for r in csv.DictReader(f):
            ev.append((int(r["step"]), r["event"], r["species"], int(r["nC"]), int(r["nO"])))
    return ev


def cumulative(ev, step):
    out = Counter()
    for s, kind, sp, nc, no in ev:
        if s > step:
            break
        if kind == "insert":
            out["C_ins"] += nc
            out["O_ins"] += no
        else:
            out["C_rem"] += nc
            out["O_rem"] += no
            key = "C_rem_CO" if sp == "CO" else "C_rem_CO2" if sp == "CO2" else "C_rem_other"
            out[key] += nc
    return out


def frame_metrics(types, pos, box):
    pos, _ = unwrap_about_ni(pos, types, box)
    m = Counter()
    iC, iO = np.where(types == C)[0], np.where(types == O)[0]
    m["nC"], m["nO"] = len(iC), len(iO)
    ni = pos[types == NI]
    hull = Delaunay(ni)
    comps, pairs = fragments(pos, types)
    bound = np.zeros(len(types), bool)
    for comp in comps:
        if (types[comp] == NI).any():
            bound[comp] = True
    inside = hull.find_simplex(pos) >= 0
    for el, idx in (("C", iC), ("O", iO)):
        for i in idx:
            m[f"{el}_" + ("interior" if inside[i] else "surface" if bound[i] else "detached")] += 1
    # C-C graph
    gC = nx.Graph()
    gC.add_nodes_from(int(i) for i in iC)
    for i, j in pairs.tolist():
        ti, tj = types[i], types[j]
        if ti == C and tj == C:
            gC.add_edge(i, j)
        elif {ti, tj} == {C, O}:
            m["CO_bonds"] += 1
    largest = 0
    for comp in nx.connected_components(gC):
        k = len(comp)
        largest = max(largest, k)
        m["C1" if k == 1 else "C2" if k == 2 else "C3_5" if k < 6 else "C6p"] += k
    m["largest"] = largest
    deg = np.array([gC.degree(int(i)) for i in iC]) if len(iC) else np.array([])
    m["sp"], m["sp2"], m["sp3"] = int((deg == 2).sum()), int((deg == 3).sum()), int((deg >= 4).sum())
    sizes = Counter(len(c) for c in nx.minimum_cycle_basis(gC)) if gC.number_of_edges() else Counter()
    m["r5"], m["r6"], m["r7"] = sizes[5], sizes[6], sizes[7]
    m["rings"] = sum(sizes.values())
    return m, sizes


COLS = ["t_ps", "step", "C_ins", "O_ins", "C_rem", "O_rem", "C_rem_CO", "C_rem_CO2", "C_rem_other",
        "nC", "nO", "C_interior", "C_surface", "C_detached", "O_interior", "O_surface", "O_detached",
        "C1", "C2", "C3_5", "C6p", "largest", "r5", "r6", "r7", "rings", "sp", "sp2", "sp3", "CO_bonds"]


def analyze(case, stride=1, force=False, dt_fs=0.25):
    traj = os.path.join(case, "traj.lammpstrj")
    out = os.path.join(case, "metrics.csv")
    if not os.path.exists(traj):
        return None
    if (not force and os.path.exists(out) and os.path.getmtime(out) > os.path.getmtime(traj)):
        return out
    ev = read_events(os.path.join(case, "events.csv"))
    rows, last_sizes = [], Counter()
    for k, (step, box, ids, types, pos) in enumerate(iter_dump(traj)):
        if k % stride:
            continue
        m, last_sizes = frame_metrics(types, pos, box)
        m.update(cumulative(ev, step))
        m["t_ps"], m["step"] = step * dt_fs / 1000, step
        rows.append([m.get(c, 0) for c in COLS])
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLS)
        w.writerows(rows)
    with open(os.path.join(case, "rings_final.csv"), "w") as f:
        f.write("ring_size,count\n")
        for s in sorted(last_sizes):
            f.write(f"{s},{last_sizes[s]}\n")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    print(analyze(a.case_dir, a.stride, a.force))
