"""Network reduction when the product distribution is known.

Reduced species set for a list of observed products P:
  * the C1 entry chemistry (all C<=1, O<=2 species: CO2, CO, HCOOH, CH2O, CH3OH, ...)
  * "skeleton closure" of P: every connected sub-skeleton of every product, at every
    hydrogenation level between fully saturated and the product's own bond orders
    (aromatic bonds may be single or double, so cyclohexene/cyclohexadiene-type
    hydro-aromatics are kept). Computed as the downward closure under
    "cut one bond (keep each fragment)" and "lower one bond order by 1".
  * optional "+O": each skeleton fragment with one extra O on any carbon (C-OH or C=O),
    i.e. the oxygenates of a CO-insertion / alcohol-route growth mechanism.
  * all monoradicals of the above.
Reactions: all unimolecular and association families inside this set; for the
bimolecular H-/atom-transfer families the partner Q is restricted to the C0-C1
radical pool (H, OH, CH3, CH2OH, CH3O, HCO, HOCO, HCOO, ...).

usage: python prune.py
"""
from __future__ import annotations

import json
import os
from math import comb

from rdkit import Chem

from chem_space import BT, analyze, bond_order, check, counts, kekule_forms
from reactions import build

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _canon(rw):
    m = Chem.Mol(rw)
    for a in m.GetAtoms():
        a.SetNoImplicit(False)
        a.SetNumExplicitHs(0)
        a.SetNumRadicalElectrons(0)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return []
    return [Chem.MolToSmiles(f) for f in Chem.GetMolFrags(m, asMols=True)]


def skeleton_closure(products):
    seen = set()
    stack = [Chem.MolToSmiles(Chem.MolFromSmiles(p)) for p in products]
    while stack:
        s = stack.pop()
        if s in seen:
            continue
        seen.add(s)
        mol = Chem.MolFromSmiles(s)
        for km in kekule_forms(mol):
            for b in km.GetBonds():
                i, j, o = b.GetBeginAtomIdx(), b.GetEndAtomIdx(), bond_order(b)
                rw = Chem.RWMol(km)
                rw.RemoveBond(i, j)
                for f in _canon(rw):
                    if f not in seen:
                        stack.append(f)
                if o > 1:
                    rw = Chem.RWMol(km)
                    rw.GetBondBetweenAtoms(i, j).SetBondType(BT[o - 1])
                    for f in _canon(rw):
                        if f not in seen:
                            stack.append(f)
    return seen


def add_one_oxygen(species, omax=2):
    out = set()
    for s in species:
        mol = Chem.MolFromSmiles(s)
        if counts(mol)[1] >= omax:
            continue
        for km in kekule_forms(mol):
            for a in km.GetAtoms():
                if a.GetSymbol() != "C":
                    continue
                h = a.GetTotalNumHs()
                for o in (1, 2):
                    if h < o:
                        continue
                    rw = Chem.RWMol(km)
                    j = rw.AddAtom(Chem.Atom("O"))
                    rw.AddBond(a.GetIdx(), j, BT[o])
                    if not check(rw)[0]:
                        continue
                    for f in _canon(rw):
                        out.add(f)
    return out


def c1_entry():
    from small_network import species_upto
    return set(species_upto(1))


def alk(n):
    return "C" * n


SCENARIOS = {
    "A_methanation (Ni)": ["C", "O", "CC"],
    "B_methanol (Cu/ZnO, In2O3)": ["CO", "COC", "COC=O", "O"],
    "C_CO2-FTS (Fe-K/Na)": (
        [alk(n) for n in range(1, 13)]
        + ["C=C" + alk(n - 2) for n in range(2, 13)]
        + ["CC=C" + alk(n - 3) for n in range(4, 13)]
        + ["CC(C)" + alk(n - 3) for n in range(4, 13)]
        + [alk(n) + "O" for n in range(1, 7)] + ["O"]),
    "D_CO2-to-aromatics (oxide/zeolite)": [
        "CO", "COC", "O", "C=C", "C=CC", "C=CCC", "CC=CC", "C=C(C)C", "CC", "CCC", "CCCC", "CC(C)C",
        "c1ccccc1", "Cc1ccccc1", "Cc1ccccc1C", "Cc1cccc(C)c1", "Cc1ccc(C)cc1", "CCc1ccccc1",
        "Cc1cccc(C)c1C", "Cc1ccc(C)c(C)c1", "Cc1cc(C)cc(C)c1", "CCc1ccccc1C", "CCc1cccc(C)c1",
        "CCc1ccc(C)cc1", "CCCc1ccccc1", "CC(C)c1ccccc1", "Cc1ccc(C)c(C)c1C", "Cc1cc(C)c(C)c(C)c1",
        "Cc1cc(C)c(C)cc1C", "Cc1cc(C)c(C)c(C)c1C", "Cc1c(C)c(C)c(C)c(C)c1C",
        "c1ccc2ccccc2c1", "Cc1cccc2ccccc12", "Cc1ccc2ccccc2c1",
        "Cc1cccc2c(C)cccc12", "Cc1ccc2cc(C)ccc2c1", "Cc1ccc2ccc(C)cc2c1"],
}


def run(name, products, with_ox):
    sk = skeleton_closure(products)
    sp = set(sk) | c1_entry()
    if with_ox:
        sp |= add_one_oxygen(sk)
    sp = {s for s in sp if check(Chem.RWMol(Chem.MolFromSmiles(s)))[0]}
    an = [analyze(x, want_products=True, want_birad=False) for x in sorted(sp)
          if counts(Chem.MolFromSmiles(x))[0] <= 1]
    c1rads = sorted({"[H]"} | {r for a in an if a["C"] <= 1 for r in a["products"]["rad"]})
    net, info = build(sorted(sp), pool=c1rads, quadratic=True, store_quadratic=len(sp) < 20000)
    nrad = info["radicals"]
    q = {"h_abstr": comb(nrad, 2), "dispro": info["beta_H_total"] * nrad}
    lin = sum(v for k, v in net.fam.items() if k not in ("h_abstr", "dispro", "sh2"))
    out = {"products": len(products), "closed_shell": info["closed"], "radicals": info["radicals"],
           "species": info["closed"] + info["radicals"], "pool": c1rads,
           "reactions_pool": len(net), "linear_reactions": lin, "families": dict(net.fam),
           "quadratic_unrestricted": q,
           "reactions_unrestricted_approx": lin + q["h_abstr"] + q["dispro"] + net.fam.get("sh2", 0),
           "maxC": max(counts(Chem.MolFromSmiles(s))[0] for s in sp)}
    tag = name.split(" ")[0] + ("_O" if with_ox else "")
    from reactions import write_csv
    write_csv(net, os.path.join(DATA, f"network_reduced_{tag}.csv" + (".gz" if len(net.rxn) > 50000 else "")))
    with open(os.path.join(DATA, f"species_reduced_{tag}.txt"), "w") as f:
        f.write("\n".join(sorted(net.species())) + "\n")
    return tag, out


if __name__ == "__main__":
    res = {}
    for name, prods in SCENARIOS.items():
        for ox in (False, True):
            tag, out = run(name, prods, ox)
            out["scenario"] = name
            res[tag] = out
            print(tag, {k: v for k, v in out.items() if k not in ("pool", "families")}, flush=True)
    with open(os.path.join(DATA, "reduced_networks.json"), "w") as f:
        json.dump(res, f, indent=1)
