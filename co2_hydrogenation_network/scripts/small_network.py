"""Explicit full networks for small carbon caps (validation of the counting formulas
and an inspectable reaction list).  usage: python small_network.py 2 3"""
import gzip
import json
import os
import sys

from rdkit import Chem

from chem_space import counts
from reactions import build, write_csv

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def species_upto(cmax, omax=2, tag="CHO_O2"):
    out = []
    for h in range(1, cmax + omax + 1):
        fn = os.path.join(DATA, "species", f"{tag}_h{h}.smi.gz")
        for s in gzip.open(fn, "rt").read().split():
            c, o = counts(Chem.MolFromSmiles(s))
            if c <= cmax and o <= omax:
                out.append(s)
    return out


if __name__ == "__main__":
    summary = {}
    for cmax in map(int, sys.argv[1:]):
        sp = species_upto(cmax)
        net, info = build(sp, cmax=cmax)
        write_csv(net, os.path.join(DATA, f"network_full_C{cmax}.csv"))
        with open(os.path.join(DATA, f"species_full_C{cmax}.txt"), "w") as f:
            f.write("\n".join(sorted(net.species())) + "\n")
        summary[cmax] = {"info": info, "species_in_network": len(net.species()),
                         "reactions": len(net.rxn), "families": dict(net.fam)}
        print(cmax, json.dumps(summary[cmax]))
    with open(os.path.join(DATA, "small_networks_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
