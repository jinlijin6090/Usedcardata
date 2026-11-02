"""Exhaustively enumerate C/H/O closed-shell molecules level by level (heavy atoms)
and aggregate per-(C,O) statistics of radicals and elementary-reaction channels.

usage: python enumerate_space.py --hmax 7 --cmax 12 --omax 2 [--tag full] [--no-analyze-above H]
Writes data/species/<tag>_h<h>.smi.gz and data/stats_<tag>.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import time
from collections import defaultdict
from multiprocessing import Pool

from chem_space import analyze, children

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_OPT = {}


def _init(omax, cmax, sym):
    _OPT.update(omax=omax, cmax=cmax, sym=sym)


def _children(smi):
    return children(smi, _OPT["omax"], _OPT["cmax"], _OPT["sym"])


def _analyze(smi):
    r = analyze(smi, want_products=False)
    r["sub_hist"] = {f"{c},{o}": v for (c, o), v in r["sub_hist"].items()}
    return r


KEYS = ["rad", "recomb", "beta_heavy", "beta_H", "ring_open", "co_add", "hshift", "hshift_near",
        "elim4", "taut13", "elim11", "rda", "birad"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmax", type=int, default=7)
    ap.add_argument("--cmax", type=int, default=12)
    ap.add_argument("--omax", type=int, default=2)
    ap.add_argument("--tag", default="full")
    ap.add_argument("--analyze-max", type=int, default=99, help="skip channel analysis above this h")
    ap.add_argument("--save-max", type=int, default=99, help="do not write SMILES above this h")
    ap.add_argument("--nosym", action="store_true")
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    a = ap.parse_args()
    os.makedirs(os.path.join(DATA, "species"), exist_ok=True)

    seeds = ["C"] + (["O"] if a.omax >= 1 else [])
    level = {1: set(seeds)}
    stats = defaultdict(lambda: defaultdict(int))
    timing = {}
    with Pool(a.procs, initializer=_init, initargs=(a.omax, a.cmax, not a.nosym)) as pool:
        for h in range(1, a.hmax + 1):
            t0 = time.time()
            if h > 1:
                cur = set()
                for ch in pool.imap_unordered(_children, sorted(level[h - 1]), chunksize=64):
                    cur |= ch
                level[h] = cur
                del level[h - 1]  # keep memory low; parents no longer needed
            t1 = time.time()
            mols = sorted(level[h])
            if h <= a.save_max:
                with gzip.open(os.path.join(DATA, "species", f"{a.tag}_h{h}.smi.gz"), "wt") as f:
                    f.write("\n".join(mols) + "\n")
            nmol = defaultdict(int)
            if h <= a.analyze_max:
                for r in pool.imap_unordered(_analyze, mols, chunksize=64):
                    k = f"{r['C']},{r['O']}"
                    s = stats[k]
                    s["mol"] += 1
                    s["mol_p"] += int(r["plaus"])
                    for key in KEYS:
                        s[key] += r[key]
                        s[key + "_p"] += r[key + "_p"]
                    for hk, v in r["sub_hist"].items():
                        s["subA:" + hk] += v
                        if r["plaus"]:
                            s["subA_p:" + hk] += v
            else:
                from rdkit import Chem
                for smi in mols:
                    m = Chem.MolFromSmiles(smi)
                    c = sum(1 for x in m.GetAtoms() if x.GetSymbol() == "C")
                    nmol[f"{c},{h - c}"] += 1
                for k, v in nmol.items():
                    stats[k]["mol"] += v
            timing[h] = (round(t1 - t0, 1), round(time.time() - t1, 1), len(mols))
            print(f"h={h}: {len(mols)} molecules  gen {t1 - t0:.1f}s  analyze {time.time() - t1:.1f}s",
                  flush=True)
            out = {"args": vars(a), "timing": timing,
                   "stats": {k: dict(v) for k, v in sorted(stats.items())}}
            with open(os.path.join(DATA, f"stats_{a.tag}.json"), "w") as f:
                json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
