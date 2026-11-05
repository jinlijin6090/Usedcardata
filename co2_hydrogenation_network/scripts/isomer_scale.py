"""Scale of an isomer-resolved *mechanistic* CO2-hydrogenation network up to C_max.

Unlike chem_space.py (all valence-allowed C/H/O graphs), this counts only the
species classes that appear in the one-page network (isomer_network.js) and
enumerates the reaction families used there, explicitly, for C1..C_max:

  aliphatic   alkanes, alkenes (+E/Z), dienes, naphthenes (5/6 rings),
              cycloalkenes, cyclodienes
  carbenium   every protonated olefin/diene/cycloolefin (= alkyl, alkenyl,
              cycloalkyl ions); radicals and adsorbed alkyls have the same
              topology classes and are counted from the same sets
  oxygenates  alcohols, aldehydes, ketones, acids (and alkoxy*, acyl*)
  aromatics   alkylbenzenes, indanes, tetralins, naphthalenes (sat. side chains)

Reaction families enumerated explicitly (species-level, unordered, deduplicated):
  zeolite     protonation, 1,2-H shift, alkyl shift (type A), PCP (type B),
              beta-scission / oligomerisation, SMS methylation, ring closure
              of alkenyl ions; arene methylation, ring isomerisation,
              Friedel-Crafts alkylation / dealkylation; HT counted by formula
  metal       beta-H elimination, CH2 insertion, CO insertion, hydrogenation
  radicals    beta-scission, H-abstraction by a small pool (H, CH3, OH)

Usage: python isomer_scale.py [CMAX]   -> data/isomer_scale_C{CMAX}.json
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

_BRACKET = re.compile(r"\[[^\]]*\]")


def has_ring(s) -> bool:
    """ring-closure digits outside bracket atoms ([CH2+] / [CH3+] carry H-count digits, not rings)"""
    return any(ch.isdigit() for ch in _BRACKET.sub("", s))

RDLogger.DisableLog("rdApp.*")
CMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def canon(m) -> str | None:
    try:
        m2 = Chem.Mol(m)
        Chem.SanitizeMol(m2)
        return Chem.MolToSmiles(m2)
    except Exception:
        return None


def mol(s):
    m = Chem.MolFromSmiles(s)
    for a in m.GetAtoms():  # let RDKit recompute H after every edit (bracket atoms default to fixed H)
        a.SetNoImplicit(False)
        a.SetNumExplicitHs(0)
    m.UpdatePropertyCache(False)
    Chem.FastFindRings(m)
    return m


def nC(m) -> int:
    return sum(a.GetSymbol() == "C" for a in m.GetAtoms())


def classes(m) -> dict[int, int]:
    """symmetry class -> representative atom index"""
    rk = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    rep = {}
    for i, r in enumerate(rk):
        rep.setdefault(r, i)
    return rep


def grow(seeds: list[str], cmax: int, ok=lambda a: True) -> dict[int, set[str]]:
    """add one sp3 carbon at a time to any atom with H that passes ok()"""
    out = defaultdict(set)
    frontier = {s for s in seeds if nC(mol(s)) <= cmax}
    for s in frontier:
        out[nC(mol(s))].add(s)
    while frontier:
        nxt = set()
        for s in frontier:
            m = mol(s)
            if nC(m) >= cmax:
                continue
            for i in classes(m).values():
                a = m.GetAtomWithIdx(i)
                if a.GetTotalNumHs() == 0 or not ok(a):
                    continue
                rw = Chem.RWMol(m)
                j = rw.AddAtom(Chem.Atom(6))
                rw.AddBond(i, j, Chem.BondType.SINGLE)
                c = canon(rw)
                if c and c not in out[nC(m) + 1]:
                    out[nC(m) + 1].add(c)
                    nxt.add(c)
        frontier = nxt
    return out


def add_double(smis, allow_ring_db=True) -> set[str]:
    """one more non-cumulated, non-aromatic C=C"""
    res = set()
    for s in smis:
        m = mol(s)
        for b in m.GetBonds():
            if b.GetBondType() != Chem.BondType.SINGLE or b.GetIsAromatic():
                continue
            a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
            if a1.GetIsAromatic() or a2.GetIsAromatic():
                continue
            if a1.GetTotalNumHs() == 0 or a2.GetTotalNumHs() == 0:
                continue
            if any(x.GetBondType() == Chem.BondType.DOUBLE for a in (a1, a2) for x in a.GetBonds()):
                continue
            if not allow_ring_db and b.IsInRing():
                continue
            rw = Chem.RWMol(m)
            rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
            c = canon(rw)
            if c:
                res.add(c)
    return res


def n_ez(s: str) -> int:
    m = mol(s)
    opts = StereoEnumerationOptions(unique=True, onlyUnassigned=True)
    seen = set()
    for iso in EnumerateStereoisomers(m, options=opts):
        for a in iso.GetAtoms():
            a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        seen.add(Chem.MolToSmiles(iso))
    return max(1, len(seen))


def in_domain(m) -> bool:
    """aliphatic domain: C<=CMAX, <=1 ring of size 5/6, <=2 C=C, no cumulenes, +1 at most"""
    if nC(m) > CMAX or nC(m) == 0:
        return False
    ri = m.GetRingInfo()
    if ri.NumRings() > 1 or any(len(r) not in (5, 6) for r in ri.AtomRings()):
        return False
    nd = 0
    for a in m.GetAtoms():
        d = sum(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds())
        if d > 1 or a.GetIsAromatic():
            return False
    nd = sum(b.GetBondType() == Chem.BondType.DOUBLE for b in m.GetBonds())
    return nd <= 2


def cation_idx(m) -> int:
    return next(a.GetIdx() for a in m.GetAtoms() if a.GetFormalCharge() == 1)


def is_primary(m) -> bool:
    a = m.GetAtomWithIdx(cation_idx(m))
    allyl = any(b.GetBondType() == Chem.BondType.DOUBLE for n in a.GetNeighbors() for b in n.GetBonds())
    return a.GetDegree() <= 1 and not allyl


# ---------------------------------------------------------------- species
log("growing alkanes / naphthenes / aromatics")
ALK = grow(["C"], CMAX)
NAPH = grow(["C1CCCC1", "C1CCCCC1"], CMAX, ok=lambda a: True)
# keep monocyclic only (grow never closes rings, so this holds)
ARO = grow(["c1ccccc1"], CMAX)
ARO2 = grow(["c1ccc2c(c1)CCC2", "c1ccc2c(c1)CCCC2", "c1ccc2ccccc2c1"], CMAX)
alkanes = set().union(*ALK.values())
naph = set().union(*NAPH.values())
arom = set().union(*ARO.values()) | set().union(*ARO2.values())

log("olefins")
alkenes = add_double(alkanes)
dienes = add_double(alkenes)
cyc1 = add_double(naph)
cyc2 = {s for s in add_double(cyc1) if all(  # both C=C endocyclic -> aromatic precursor
    b.IsInRing() for b in mol(s).GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)}
olefinic = alkenes | dienes | cyc1 | cyc2
SPEC = {"alkane": alkanes, "alkene": alkenes, "diene": dienes, "naphthene": naph,
        "cycloalkene": cyc1, "cyclodiene": cyc2, "aromatic": arom}


def by_c(ss):
    c = Counter(nC(mol(s)) for s in ss)
    return [c.get(n, 0) for n in range(1, CMAX + 1)]


# ---------------------------------------------------------------- reactions on cations
def protonate(s):
    """olefin -> set of carbenium ions (H to one end, + on the other)"""
    m = mol(s)
    out = set()
    for b in m.GetBonds():
        if b.GetBondType() != Chem.BondType.DOUBLE:
            continue
        for i in (b.GetBeginAtomIdx(), b.GetEndAtomIdx()):
            rw = Chem.RWMol(m)
            rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(i).SetFormalCharge(1)
            c = canon(rw)
            if c:
                out.add(c)
    return out


log("protonation")
prot_edges = set()
ions = set()
for s in olefinic:
    for c in protonate(s):
        prot_edges.add((s, c))
        ions.add(c)
# methyl cation and ions of alkanes not reachable from an olefin (only CH3+)
ions.add("[CH3+]")
log("ions", len(ions))


def ion_moves(s):
    """yield (family, products tuple) for unimolecular ion steps"""
    m = mol(s)
    a = cation_idx(m)
    A = m.GetAtomWithIdx(a)
    for B in A.GetNeighbors():
        b = B.GetIdx()
        bond_ab = m.GetBondBetweenAtoms(a, b)
        if bond_ab.GetBondType() != Chem.BondType.SINGLE:
            continue
        # 1,2-H shift
        if B.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m)
            rw.GetAtomWithIdx(a).SetFormalCharge(0)
            rw.GetAtomWithIdx(b).SetFormalCharge(1)
            yield "hshift", (canon(rw),)
        for C in B.GetNeighbors():
            c = C.GetIdx()
            if c == a or m.GetBondBetweenAtoms(b, c).GetBondType() != Chem.BondType.SINGLE:
                continue
            # alkyl shift (type A): c moves from b to a, + to b
            rw = Chem.RWMol(m)
            rw.RemoveBond(b, c)
            rw.AddBond(a, c, Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(a).SetFormalCharge(0)
            rw.GetAtomWithIdx(b).SetFormalCharge(1)
            yield "alkylshift", (canon(rw),)
            # PCP (type B): a moves from b to c, H c->b, + stays on a
            if C.GetTotalNumHs() > 0:
                rw = Chem.RWMol(m)
                rw.RemoveBond(a, b)
                rw.AddBond(a, c, Chem.BondType.SINGLE)
                yield "pcp", (canon(rw),)
            # beta-scission: break b-c, a=b, + on c
            rw = Chem.RWMol(m)
            rw.RemoveBond(b, c)
            rw.GetBondBetweenAtoms(a, b).SetBondType(Chem.BondType.DOUBLE)
            rw.GetAtomWithIdx(a).SetFormalCharge(0)
            rw.GetAtomWithIdx(c).SetFormalCharge(1)
            m2 = rw.GetMol()
            try:
                Chem.SanitizeMol(m2)
            except Exception:
                continue
            frags = [Chem.MolToSmiles(f) for f in Chem.GetMolFrags(m2, asMols=True)]
            yield ("beta" if len(frags) == 2 else "ringopen"), tuple(frags)
    # ring closure of alkenyl ions: a attacks C=C (d=e) -> 5/6 ring, + on far end
    if A.IsInRing():
        return
    for bnd in m.GetBonds():
        if bnd.GetBondType() != Chem.BondType.DOUBLE:
            continue
        for d, e in ((bnd.GetBeginAtomIdx(), bnd.GetEndAtomIdx()), (bnd.GetEndAtomIdx(), bnd.GetBeginAtomIdx())):
            path = Chem.GetShortestPath(m, a, d)
            if len(path) not in (5, 6):
                continue
            rw = Chem.RWMol(m)
            rw.AddBond(a, d, Chem.BondType.SINGLE)
            rw.GetBondBetweenAtoms(d, e).SetBondType(Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(a).SetFormalCharge(0)
            rw.GetAtomWithIdx(e).SetFormalCharge(1)
            yield "cyclise", (canon(rw),)


log("ion rearrangements / scission")
E = defaultdict(set)
ions_l = sorted(ions)
for k, s in enumerate(ions_l):
    if k % 20000 == 0:
        log(f"  {k}/{len(ions_l)}")
    if s == "[CH3+]":
        continue
    for fam, prods in ion_moves(s):
        if any(p is None for p in prods):
            continue
        pm = [mol(p) for p in prods]
        if not all(in_domain(x) for x in pm):
            continue
        if fam in ("hshift", "alkylshift", "pcp", "cyclise", "ringopen"):
            p = prods[0]
            if p == s or p not in ions:
                continue
            E[fam].add(frozenset((s, p)))
        elif fam == "beta":
            ion_p = next(p for p in prods if "+" in p)
            ole = next(p for p in prods if "+" not in p)
            if ole not in olefinic or ion_p not in ions:
                continue
            fam2 = "methylation (SMS)" if ion_p == "[CH3+]" else "beta-scission/oligo"
            E[fam2].add((s, ion_p, ole))

# primary-free variant (classical zeolite chemistry: no primary carbenium ions)
prim = {s for s in ions if s != "[CH3+]" and is_primary(mol(s))}


def no_prim(fam, x):
    items = list(x) if isinstance(x, frozenset) else [x[0], x[1]] if fam != "prot" else [x[1]]
    return not any(i in prim for i in items)


Z = {"protonation": len(prot_edges)}
Z.update({k: len(v) for k, v in E.items()})
Z_np = {"protonation": sum(c not in prim for _, c in prot_edges)}
for k, v in E.items():
    Z_np[k] = sum(no_prim(k, x) for x in v)
n_ion = len(ions) - 1
n_ion_np = n_ion - len(prim)

# ---------------------------------------------------------------- aromatics
log("aromatic reactions")
aro_meth, aro_iso, aro_fc = set(), set(), set()
alkyl_ions = {s for s in ions if not has_ring(s) and "=" not in s and s not in prim}
ions_by_c = defaultdict(list)
for s in alkyl_ions:
    ions_by_c[nC(mol(s))].append(s)
arene_ring_sites = 0
for s in arom:
    m = mol(s)
    n = nC(m)
    for i in classes(m).values():
        A = m.GetAtomWithIdx(i)
        if not A.GetIsAromatic():
            continue
        if A.GetTotalNumHs() > 0:
            arene_ring_sites += 1
            # methylation / alkylation at a ring CH
            for k in range(1, CMAX - n + 1):
                for ion in ([None] if k == 1 else ions_by_c[k]):
                    rw = Chem.RWMol(m)
                    if ion is None:
                        j = rw.AddAtom(Chem.Atom(6))
                        rw.AddBond(i, j, Chem.BondType.SINGLE)
                    else:
                        im = mol(ion)
                        off = rw.GetNumAtoms()
                        for at in im.GetAtoms():
                            na = Chem.Atom(6)
                            rw.AddAtom(na)
                        for bd in im.GetBonds():
                            rw.AddBond(off + bd.GetBeginAtomIdx(), off + bd.GetEndAtomIdx(), bd.GetBondType())
                        rw.AddBond(i, off + cation_idx(im), Chem.BondType.SINGLE)
                    p = canon(rw)
                    if p in arom:
                        (aro_meth if ion is None else aro_fc).add((s, ion, p))
        else:
            # 1,2-shift of the substituent to a neighbouring ring CH
            subs = [x for x in A.GetNeighbors() if not x.GetIsAromatic()]
            for X in subs:
                for R in A.GetNeighbors():
                    if not R.GetIsAromatic() or R.GetTotalNumHs() == 0:
                        continue
                    rw = Chem.RWMol(m)
                    rw.RemoveBond(i, X.GetIdx())
                    rw.AddBond(R.GetIdx(), X.GetIdx(), Chem.BondType.SINGLE)
                    p = canon(rw)
                    if p in arom and p != s:
                        aro_iso.add(frozenset((s, p)))
aro_ions = arene_ring_sites + sum(  # arenium ions: protonation at every ring C class
    sum(1 for i in classes(mol(s)).values() if mol(s).GetAtomWithIdx(i).GetIsAromatic()
        and mol(s).GetAtomWithIdx(i).GetTotalNumHs() == 0) for s in arom)

# transalkylation (methyl transfer arene -> arene), counted by carbon number
donors, accept = Counter(), Counter()
for s in arom:
    m = mol(s)
    n = nC(m)
    for i in classes(m).values():
        A = m.GetAtomWithIdx(i)
        if A.GetIsAromatic() and A.GetTotalNumHs() > 0:
            accept[n] += 1
        if (not A.GetIsAromatic()) and A.GetDegree() == 1 and any(x.GetIsAromatic() for x in A.GetNeighbors()):
            donors[n] += 1
transalk = sum(donors[n] * accept[k] for n in donors for k in accept if n - 1 >= 6 and k + 1 <= CMAX) // 2

# ---------------------------------------------------------------- metal / oxygenates / radicals
log("metal + oxygenates + radicals")
site_all = site_ch3 = site_ch2 = 0
for s in alkanes:
    m = mol(s)
    for i in classes(m).values():
        h = m.GetAtomWithIdx(i).GetTotalNumHs()
        if h:
            site_all += 1
        site_ch3 += h == 3 or s == "C"
        site_ch2 += h == 2
n_alkyl = site_all                     # R*, R*, R+ topology classes (incl. primary)
oxy = {"alcohol (ROH) / alkoxy*": site_all,
       "aldehyde (RCHO) / acyl*": site_ch3,
       "ketone": site_ch2, "carboxylic acid": site_ch3}
# aldehyde with C carbons = CH3 sites of alkane with C carbons (CHO replaces CH3); same for acids
beta_h = sum(len(protonate(s)) for s in alkenes)   # alkyl* <-> olefin* pairs (incl. primary)
metal = {
    "alkyl* + H* -> alkane": n_alkyl,
    "beta-H elimination alkyl* <-> olefin*": beta_h,
    "CH2* insertion (chain growth, all branches)": sum(1 for s in alkanes for i in classes(mol(s)).values()
                                                       if mol(s).GetAtomWithIdx(i).GetTotalNumHs() and nC(mol(s)) < CMAX),
    "olefin adsorption (pi / di-sigma)": len(alkenes),
    "CO insertion -> acyl*": sum(1 for s in alkanes for i in classes(mol(s)).values()
                                 if mol(s).GetAtomWithIdx(i).GetTotalNumHs() and nC(mol(s)) < CMAX),
    "acyl* -> aldehyde -> alkoxy* -> alcohol (3 steps)": 3 * oxy["aldehyde (RCHO) / acyl*"],
    "alcohol dehydration / ketone hydrogenation": site_all + site_ch2,
}
rad = {
    "radicals R. (alkyl, all sites)": n_alkyl,
    "beta-scission R. -> olefin + R'.": Z["beta-scission/oligo"] + Z.get("methylation (SMS)", 0),
    "H-abstraction by H./CH3./OH (3 per site)": 3 * n_alkyl,
    "recombination with H./CH3.": 2 * n_alkyl,
    "1,2-H shift (R. <-> R'.)": Z["hshift"],
}

# hydride transfer (bimolecular, R+ + R'H): one reaction per unordered ion pair
ht_pairs = n_ion_np * (n_ion_np - 1) // 2

species = {k: by_c(v) for k, v in SPEC.items()}
species["alkene E/Z resolved"] = None
log("E/Z count")
ez = Counter()
for s in alkenes | dienes:
    ez[("alkene" if s in alkenes else "diene", nC(mol(s)))] += n_ez(s)
species["alkene E/Z resolved"] = [ez[("alkene", n)] for n in range(1, CMAX + 1)]
species["diene E/Z resolved"] = [ez[("diene", n)] for n in range(1, CMAX + 1)]
species["carbenium ions (all)"] = by_c(ions - {"[CH3+]"})
species["carbenium ions (no primary)"] = by_c(ions - prim - {"[CH3+]"})

out = {
    "cmax": CMAX,
    "species_by_carbon_number": species,
    "totals": {
        "closed-shell hydrocarbons": sum(len(v) for v in SPEC.values()),
        "closed-shell hydrocarbons, E/Z resolved": len(alkanes) + len(naph) + len(cyc1) + len(cyc2) + len(arom)
        + sum(ez.values()),
        "oxygenates": oxy,
        "carbenium ions (all)": n_ion,
        "carbenium ions (no primary)": n_ion_np,
        "arenium ions": aro_ions,
        "alkyl radicals = alkyl* = alkoxy* classes": n_alkyl,
    },
    "reactions": {
        "zeolite (all ions)": Z,
        "zeolite (no primary ions)": Z_np,
        "zeolite hydride transfer (ion pairs, formula)": ht_pairs,
        "aromatics": {"methylation": len(aro_meth), "ring 1,2-shift isomerisation": len(aro_iso),
                      "Friedel-Crafts alkylation/dealkylation": len(aro_fc),
                      "transalkylation (methyl, formula)": transalk},
        "metal": metal,
        "radicals": rad,
    },
    "runtime_s": round(time.time() - T0, 1),
}
Path(__file__).resolve().parents[1].joinpath("data", f"isomer_scale_C{CMAX}.json").write_text(
    json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps(out, indent=1, ensure_ascii=False))


# ---------------------------------------------------------------- explicit graph dump (DUMP=path.json.gz)
def _dump_graph(path):
    import gzip

    S, D = Chem.BondType.SINGLE, Chem.BondType.DOUBLE
    nodes, edges = {}, set()

    def skel(s):
        m = Chem.RWMol(mol(s))
        for b in m.GetBonds():
            if not b.GetIsAromatic():
                b.SetBondType(S)
        for a in m.GetAtoms():
            a.SetFormalCharge(0)
        for i in sorted((a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O"), reverse=True):
            m.RemoveAtom(i)
        return canon(m) or s

    def node(cls, s):
        nid = f"{cls}|{s}"
        if nid not in nodes:
            nodes[nid] = ("hub", 1, s) if cls == "hub" else (cls, nC(mol(s)), skel(s))
        return nid

    def edge(u, v, fam):
        if u != v:
            edges.add((min(u, v), max(u, v), fam))

    def edit(s, fn):
        rw = Chem.RWMol(mol(s))
        try:
            fn(rw)
        except Exception:
            return None
        return canon(rw)

    for cls, ss in SPEC.items():
        for s in ss:
            node(cls, s)
    hub = {k: node("hub", k) for k in ("CH2*", "CO*", "CH3OH/SMS")}
    # ---- carbenium ions (no primary)
    def ion_cls(s):
        return "ion_cyc" if has_ring(s) else "ion_alkenyl" if "=" in s else "ion_alkyl"
    for o, c in prot_edges:
        if c not in prim:
            edge(node(ion_cls(c), c), next(f"{k}|{o}" for k in SPEC if o in SPEC[k]), "Z:prot")
    for fam, v in E.items():
        if fam == "ringopen":
            continue
        for x in v:
            if not no_prim(fam, x):
                continue
            if fam in ("beta-scission/oligo", "methylation (SMS)"):
                s, ip, ole = x
                on = next(f"{k}|{ole}" for k in SPEC if ole in SPEC[k])
                sn = node(ion_cls(s), s)
                edge(sn, on, "Z:beta")
                edge(sn, hub["CH3OH/SMS"] if ip == "[CH3+]" else node(ion_cls(ip), ip), "Z:beta" if ip != "[CH3+]" else "Z:meth")
            else:
                a, b = tuple(x)
                edge(node(ion_cls(a), a), node(ion_cls(b), b), "Z:" + fam)
    # ---- alkyl sites: adsorbed alkyl R*, radical R., oxygenates
    sites = set()
    for s in alkanes:
        m = mol(s)
        for i in classes(m).values():
            if m.GetAtomWithIdx(i).GetTotalNumHs():
                rw = Chem.RWMol(m)
                rw.GetAtomWithIdx(i).SetFormalCharge(1)
                sites.add((canon(rw), s))
    for key, alk in sites:
        ads, rad = node("ads", key), node("rad", key)
        an = f"alkane|{alk}"
        edge(ads, an, "M:hyd")
        edge(rad, an, "R:abs")
        m = mol(key)
        a = cation_idx(m)
        A = m.GetAtomWithIdx(a)
        nh = A.GetTotalNumHs()
        # beta-H elimination -> alkene
        for B in A.GetNeighbors():
            if B.GetTotalNumHs():
                def f(rw, b=B.GetIdx()):
                    rw.GetBondBetweenAtoms(a, b).SetBondType(D)
                    rw.GetAtomWithIdx(a).SetFormalCharge(0)
                ol = edit(key, f)
                if ol in alkenes:
                    edge(ads, f"alkene|{ol}", "M:betaH")
        if nC(m) < CMAX:
            def ins(rw):
                j = rw.AddAtom(Chem.Atom(6))
                rw.AddBond(a, j, S)
                rw.GetAtomWithIdx(a).SetFormalCharge(0)
                rw.GetAtomWithIdx(j).SetFormalCharge(1)
            k2 = edit(key, ins)
            if k2:
                edge(ads, node("ads", k2), "M:CH2")
                edge(hub["CH2*"], node("ads", k2), "M:hub")

            def coi(rw):
                j = rw.AddAtom(Chem.Atom(6))
                o = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(a, j, S)
                rw.AddBond(j, o, D)
                rw.GetAtomWithIdx(a).SetFormalCharge(0)
            ald = edit(key, coi)
            if ald:
                edge(ads, node("aldehyde", ald), "M:CO")
                edge(hub["CO*"], node("aldehyde", ald), "M:hub")
        # oxygenates on this site
        def oh(rw):
            o = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(a, o, S)
            rw.GetAtomWithIdx(a).SetFormalCharge(0)
        alc = node("alcohol", edit(key, oh))
        for B in A.GetNeighbors():
            if B.GetTotalNumHs():
                def f(rw, b=B.GetIdx()):
                    rw.GetBondBetweenAtoms(a, b).SetBondType(D)
                    rw.GetAtomWithIdx(a).SetFormalCharge(0)
                ol = edit(key, f)
                if ol in alkenes:
                    edge(alc, f"alkene|{ol}", "O:dehyd")
        def co(rw):
            o = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(a, o, D)
            rw.GetAtomWithIdx(a).SetFormalCharge(0)
        if nh == 2 or (nh == 3 and A.GetDegree() == 0):     # CH3 site -> aldehyde, acid
            ald = node("aldehyde", edit(key, co))
            edge(ald, alc, "O:hyd")
            def acid(rw):
                o = rw.AddAtom(Chem.Atom(8))
                o2 = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(a, o, D)
                rw.AddBond(a, o2, S)
                rw.GetAtomWithIdx(a).SetFormalCharge(0)
            edge(ald, node("acid", edit(key, acid)), "O:ox")
        elif nh == 1 and A.GetDegree() == 2:                # CH2 site -> ketone
            edge(node("ketone", edit(key, co)), alc, "O:hyd")
        # radical beta-scission
        for fam, prods in ion_moves(key):
            if fam != "beta" or any(p is None for p in prods):
                continue
            ip = next(p for p in prods if "+" in p)
            ole = next(p for p in prods if "+" not in p)
            if ole in alkenes and not has_ring(ip) and "=" not in ip:
                edge(rad, f"alkene|{ole}", "R:beta")
                edge(rad, node("rad", ip), "R:beta")
    # ---- hydrogenation / hydrogen transfer between closed-shell classes
    def parents(s):
        m = mol(s)
        out = set()
        for b in m.GetBonds():
            if b.GetBondType() == D:
                rw = Chem.RWMol(m)
                rw.GetBondWithIdx(b.GetIdx()).SetBondType(S)
                out.add(canon(rw))
        return out
    for s in alkenes:
        for p in parents(s):
            edge(f"alkene|{s}", f"alkane|{p}", "H:hyd")
    for s in dienes:
        for p in parents(s):
            if p in alkenes:
                edge(f"diene|{s}", f"alkene|{p}", "Z:HT")
    for s in cyc1:
        for p in parents(s):
            if p in naph:
                edge(f"cycloalkene|{s}", f"naphthene|{p}", "H:hyd")
    for s in cyc2:
        for p in parents(s):
            if p in cyc1:
                edge(f"cyclodiene|{s}", f"cycloalkene|{p}", "Z:HT")
        m = mol(s)
        ring = m.GetRingInfo().AtomRings()[0]
        if len(ring) == 6:
            rw = Chem.RWMol(m)
            rs = set(ring)
            for b in rw.GetBonds():
                if b.GetBeginAtomIdx() in rs and b.GetEndAtomIdx() in rs:
                    b.SetBondType(Chem.BondType.AROMATIC)
                    b.SetIsAromatic(True)
            for i in rs:
                rw.GetAtomWithIdx(i).SetIsAromatic(True)
            ar = canon(rw)
            if ar in arom:
                edge(f"cyclodiene|{s}", f"aromatic|{ar}", "Z:HT")
    # ---- aromatics
    for s, _, p in aro_meth:
        edge(f"aromatic|{s}", f"aromatic|{p}", "A:meth")
        edge(hub["CH3OH/SMS"], f"aromatic|{p}", "A:hub")
    for x in aro_iso:
        a, b = tuple(x)
        edge(f"aromatic|{a}", f"aromatic|{b}", "A:iso")
    for s, ion, p in aro_fc:
        edge(f"aromatic|{s}", f"aromatic|{p}", "A:FC")
        edge(node(ion_cls(ion), ion), f"aromatic|{p}", "A:FC")
    edges_ok = [e for e in edges if e[0] in nodes and e[1] in nodes]
    data = {"cmax": CMAX, "nodes": [[k, *v] for k, v in nodes.items()], "edges": [list(e) for e in edges_ok]}
    with gzip.open(path, "wt") as fh:
        json.dump(data, fh)
    log("graph dumped", len(nodes), "nodes", len(edges_ok), "edges", Counter(e[2] for e in edges_ok).most_common())


if os.environ.get("DUMP"):
    _dump_graph(os.environ["DUMP"])
