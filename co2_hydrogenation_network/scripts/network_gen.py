"""Explicit, isomer-resolved CO2-hydrogenation reaction network up to C_max, with every
elementary reaction recorded with all of its participants.

    python network_gen.py --cmax 16 --full-upto 12 --out ../data/net_C16_shape   # recommended
    python network_gen.py --cmax 16 --full-upto 16 --out ../data/net_C16_full    # everything

Species classes and reaction families are the ones of isomer_scale.py / the one-page decks:

  closed shell : alkane, alkene, diene, naphthene (5/6 ring), cycloalkene, cyclodiene, aromatic
  ions         : every non-primary carbenium ion from protonating any olefinic species
  sites        : adsorbed alkyl R*, gas radical R., alcohol (alkoxy*), aldehyde (acyl*), ketone, acid
  C1 core      : CO2 -> formate / CO* routes -> CH3OH, CH2*, SMS, DME, hydrocarbon pool (HCP)
  pool         : H2, H*, H+, H2O, [H] (hydride-transfer acceptor), H.

Above --full-upto carbons only ZSM-5-compatible skeletons are kept (methyl branches only,
at most two branch points, no quaternary carbon; ring systems with at most two linear
substituents). The filter is closed under removing a terminal carbon, so growing only from
kept parents still finds every kept child.

Outputs (in --out):
  species.tsv.gz    id  class  C  smiles  skeleton
  reactions.tsv.gz  id  family  reactant_ids  product_ids  equation(SMILES)
  stats.json        counts per class x carbon number and per family
Bimolecular hydride transfer (R+ + R'H) and arene-arene transalkylation are combinatorial
pair families; they are counted in stats.json, not listed one by one.
"""
from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
import os
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
S1, D2, AR = Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.AROMATIC
CMAX = 16
FULL = 12
T0 = time.time()
G: dict = {}          # membership sets, shared with forked workers


def log(*a):
    print(f"[{time.time() - T0:8.1f}s]", *a, flush=True)


# ------------------------------------------------------------------ RDKit helpers
def mol(s):
    m = Chem.MolFromSmiles(s)
    for a in m.GetAtoms():
        a.SetNoImplicit(False)
        a.SetNumExplicitHs(0)
    m.UpdatePropertyCache(False)
    Chem.FastFindRings(m)
    return m


def canon(m):
    try:
        m2 = Chem.Mol(m)
        Chem.SanitizeMol(m2)
        return Chem.MolToSmiles(m2)
    except Exception:
        return None


def nC(m):
    return sum(a.GetAtomicNum() == 6 for a in m.GetAtoms())


def reps(m):
    rk = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    out = {}
    for i, r in enumerate(rk):
        out.setdefault(r, i)
    return list(out.values())


def skeleton(m):
    rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        if not b.GetIsAromatic():
            b.SetBondType(S1)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0)
    for i in sorted((a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 8), reverse=True):
        rw.RemoveAtom(i)
    m2 = rw.GetMol()
    try:
        Chem.SanitizeMol(m2)
    except Exception:
        return None
    Chem.FastFindRings(m2)
    return m2


def shape_ok_skel(sk) -> bool:
    """ZSM-5-compatible skeleton (applied above FULL carbons)."""
    ri = sk.GetRingInfo()
    ring_atoms = {i for r in ri.AtomRings() for i in r}
    for a in sk.GetAtoms():
        if a.GetDegree() > 3 and not a.GetIsAromatic():
            return False                                  # quaternary carbon
    if not ring_atoms:
        br = [a for a in sk.GetAtoms() if a.GetDegree() == 3]
        if len(br) > 2:
            return False
        return all(sum(n.GetDegree() > 1 for n in a.GetNeighbors()) <= 2 for a in br)
    # ring systems: <= 2 substituted ring atoms, every side chain linear
    subst = [a for a in sk.GetAtoms() if a.GetIdx() in ring_atoms
             and any(n.GetIdx() not in ring_atoms for n in a.GetNeighbors())]
    if len(subst) > 2:
        return False
    return all(a.GetDegree() <= 2 for a in sk.GetAtoms() if a.GetIdx() not in ring_atoms)


def shape_ok(m) -> bool:
    if nC(m) <= FULL:
        return True
    sk = skeleton(m)
    return sk is not None and shape_ok_skel(sk)


def cation_idx(m):
    return next(a.GetIdx() for a in m.GetAtoms() if a.GetFormalCharge() == 1)


def is_primary(m) -> bool:
    a = m.GetAtomWithIdx(cation_idx(m))
    allyl = any(b.GetBondType() == D2 for n in a.GetNeighbors() for b in n.GetBonds())
    return a.GetDegree() <= 1 and not allyl


_BRACKET = re.compile(r"\[[^\]]*\]")


def has_ring(s) -> bool:
    """ring-closure digits outside bracket atoms ([CH2+] / [CH3+] carry H-count digits, not rings)"""
    return any(ch.isdigit() for ch in _BRACKET.sub("", s))


def ion_cls(s):
    return "ion_cyc" if has_ring(s) else "ion_alkenyl" if "=" in s else "ion_alkyl"


def aliph_domain(m) -> bool:
    if nC(m) > CMAX or nC(m) == 0:
        return False
    ri = m.GetRingInfo()
    if ri.NumRings() > 1 or any(len(r) not in (5, 6) for r in ri.AtomRings()):
        return False
    for a in m.GetAtoms():
        if a.GetIsAromatic() or sum(b.GetBondType() == D2 for b in a.GetBonds()) > 1:
            return False
    return sum(b.GetBondType() == D2 for b in m.GetBonds()) <= 2


# ------------------------------------------------------------------ worker functions (forked)
def w_grow(s):
    m = mol(s)
    out = []
    for i in reps(m):
        if m.GetAtomWithIdx(i).GetTotalNumHs() == 0:
            continue
        rw = Chem.RWMol(m)
        j = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(i, j, S1)
        c = canon(rw)
        if c and shape_ok(mol(c)):
            out.append(c)
    return out


def w_double(s):
    """one more non-cumulated, non-aromatic C=C"""
    m = mol(s)
    out = []
    for b in m.GetBonds():
        if b.GetBondType() != S1 or b.GetIsAromatic():
            continue
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetIsAromatic() or a2.GetIsAromatic() or a1.GetTotalNumHs() == 0 or a2.GetTotalNumHs() == 0:
            continue
        if any(x.GetBondType() == D2 for a in (a1, a2) for x in a.GetBonds()):
            continue
        rw = Chem.RWMol(m)
        rw.GetBondWithIdx(b.GetIdx()).SetBondType(D2)
        c = canon(rw)
        if c:
            out.append(c)
    return out


def w_prot(s):
    """olefinic species -> (reaction, ion) for every non-primary carbenium ion"""
    m = mol(s)
    out = []
    for b in m.GetBonds():
        if b.GetBondType() != D2:
            continue
        for i in (b.GetBeginAtomIdx(), b.GetEndAtomIdx()):
            rw = Chem.RWMol(m)
            rw.GetBondWithIdx(b.GetIdx()).SetBondType(S1)
            rw.GetAtomWithIdx(i).SetFormalCharge(1)
            c = canon(rw)
            if c and not is_primary(mol(c)):
                out.append(c)
    return s, sorted(set(out))


def ion_moves(s):
    m = mol(s)
    a = cation_idx(m)
    A = m.GetAtomWithIdx(a)
    for B in A.GetNeighbors():
        b = B.GetIdx()
        if m.GetBondBetweenAtoms(a, b).GetBondType() != S1:
            continue
        if B.GetTotalNumHs() > 0:
            rw = Chem.RWMol(m)
            rw.GetAtomWithIdx(a).SetFormalCharge(0)
            rw.GetAtomWithIdx(b).SetFormalCharge(1)
            yield "Z:hshift", (canon(rw),)
        for C in B.GetNeighbors():
            c = C.GetIdx()
            if c == a or m.GetBondBetweenAtoms(b, c).GetBondType() != S1:
                continue
            rw = Chem.RWMol(m)
            rw.RemoveBond(b, c)
            rw.AddBond(a, c, S1)
            rw.GetAtomWithIdx(a).SetFormalCharge(0)
            rw.GetAtomWithIdx(b).SetFormalCharge(1)
            yield "Z:alkylshift", (canon(rw),)
            if C.GetTotalNumHs() > 0:
                rw = Chem.RWMol(m)
                rw.RemoveBond(a, b)
                rw.AddBond(a, c, S1)
                yield "Z:pcp", (canon(rw),)
            rw = Chem.RWMol(m)
            rw.RemoveBond(b, c)
            rw.GetBondBetweenAtoms(a, b).SetBondType(D2)
            rw.GetAtomWithIdx(a).SetFormalCharge(0)
            rw.GetAtomWithIdx(c).SetFormalCharge(1)
            m2 = rw.GetMol()
            try:
                Chem.SanitizeMol(m2)
            except Exception:
                continue
            frags = tuple(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(m2, asMols=True))
            if len(frags) == 2:
                yield "beta", frags
    if A.IsInRing():
        return
    for bnd in m.GetBonds():
        if bnd.GetBondType() != D2:
            continue
        for d, e in ((bnd.GetBeginAtomIdx(), bnd.GetEndAtomIdx()), (bnd.GetEndAtomIdx(), bnd.GetBeginAtomIdx())):
            if len(Chem.GetShortestPath(m, a, d)) not in (5, 6):
                continue
            rw = Chem.RWMol(m)
            rw.AddBond(a, d, S1)
            rw.GetBondBetweenAtoms(d, e).SetBondType(S1)
            rw.GetAtomWithIdx(a).SetFormalCharge(0)
            rw.GetAtomWithIdx(e).SetFormalCharge(1)
            yield "Z:cyclise", (canon(rw),)


def okey(s):
    """class|smiles key of a neutral olefinic species, or None"""
    for k in ("alkene", "diene", "cycloalkene", "cyclodiene"):
        if s in G[k]:
            return f"{k}|{s}"
    return None


def w_ion(s):
    ions = G["ions"]
    out = set()
    me = f"{ion_cls(s)}|{s}"
    for fam, prods in ion_moves(s):
        if any(p is None for p in prods):
            continue
        if fam == "beta":
            ip = next((p for p in prods if "+" in p), None)
            ole = next((p for p in prods if "+" not in p), None)
            if ip is None or ole is None:
                continue
            ok = okey(ole)
            if ok is None:
                continue
            if ip == "[CH3+]":
                out.add(("Z:meth", (ok, "core|SMS"), (me,)))
            elif ip in ions:
                out.add(("Z:beta", (me,), (ok, f"{ion_cls(ip)}|{ip}")))
        else:
            p = prods[0]
            if p == s or p not in ions:
                continue
            if fam == "Z:cyclise":
                out.add((fam, (me,), (f"{ion_cls(p)}|{p}",)))
            elif s < p:                                   # symmetric moves: record once
                out.add((fam, (me,), (f"{ion_cls(p)}|{p}",)))
    return list(out)


def edit(s, fn):
    rw = Chem.RWMol(mol(s))
    try:
        fn(rw)
    except Exception:
        return None
    return canon(rw)


def w_site(alk):
    """metal / oxygenate / radical chemistry on every symmetry-distinct C-H site of an alkane"""
    m0 = mol(alk)
    R, spec = [], []
    an = f"alkane|{alk}"
    for i in reps(m0):
        A0 = m0.GetAtomWithIdx(i)
        if A0.GetTotalNumHs() == 0:
            continue
        rw = Chem.RWMol(m0)
        rw.GetAtomWithIdx(i).SetFormalCharge(1)
        key = canon(rw)
        ads, rad = f"ads|{key}", f"rad|{key}"
        spec += [ads, rad]
        R.append(("M:hyd", (ads, "pool|H*"), (an,)))
        R.append(("R:abs", (rad, "pool|H2"), (an, "pool|H.")))
        m = mol(key)
        a = cation_idx(m)
        A = m.GetAtomWithIdx(a)
        nh = A.GetTotalNumHs()
        n = nC(m)

        def to_db(b):
            def f(rw2):
                rw2.GetBondBetweenAtoms(a, b).SetBondType(D2)
                rw2.GetAtomWithIdx(a).SetFormalCharge(0)
            return edit(key, f)
        olefins = []
        for B in A.GetNeighbors():
            if B.GetTotalNumHs():
                ol = to_db(B.GetIdx())
                if ol in G["alkene"]:
                    olefins.append(ol)
        for ol in sorted(set(olefins)):
            R.append(("M:betaH", (ads,), (f"alkene|{ol}", "pool|H*")))
        if n < CMAX:
            def ins(rw2):
                j = rw2.AddAtom(Chem.Atom(6))
                rw2.AddBond(a, j, S1)
                rw2.GetAtomWithIdx(a).SetFormalCharge(0)
                rw2.GetAtomWithIdx(j).SetFormalCharge(1)
            k2 = edit(key, ins)
            if k2 and shape_ok(mol(k2)):
                R.append(("M:CH2", (ads, "core|CH2*"), (f"ads|{k2}",)))

            def coi(rw2):
                j = rw2.AddAtom(Chem.Atom(6))
                o = rw2.AddAtom(Chem.Atom(8))
                rw2.AddBond(a, j, S1)
                rw2.AddBond(j, o, D2)
                rw2.GetAtomWithIdx(a).SetFormalCharge(0)
            ald = edit(key, coi)
            if ald and shape_ok(mol(ald)):
                R.append(("M:CO", (ads, "core|CO*"), (f"aldehyde|{ald}",)))
                spec.append(f"aldehyde|{ald}")

        def oh(rw2):
            o = rw2.AddAtom(Chem.Atom(8))
            rw2.AddBond(a, o, S1)
            rw2.GetAtomWithIdx(a).SetFormalCharge(0)
        alc = f"alcohol|{edit(key, oh)}"
        spec.append(alc)
        R.append(("O:OH", (ads, "pool|H2O"), (alc, "pool|H*")))     # R* + OH* termination
        for ol in sorted(set(olefins)):
            R.append(("O:dehyd", (alc,), (f"alkene|{ol}", "pool|H2O")))

        def co(rw2):
            o = rw2.AddAtom(Chem.Atom(8))
            rw2.AddBond(a, o, D2)
            rw2.GetAtomWithIdx(a).SetFormalCharge(0)
        if nh == 2 or (nh == 3 and A.GetDegree() == 0):
            ald = f"aldehyde|{edit(key, co)}"
            spec.append(ald)
            R.append(("O:hyd", (ald, "pool|H2"), (alc,)))

            def acid(rw2):
                o = rw2.AddAtom(Chem.Atom(8))
                o2 = rw2.AddAtom(Chem.Atom(8))
                rw2.AddBond(a, o, D2)
                rw2.AddBond(a, o2, S1)
                rw2.GetAtomWithIdx(a).SetFormalCharge(0)
            ac = f"acid|{edit(key, acid)}"
            spec.append(ac)
            R.append(("O:ox", (ald, "pool|H2O"), (ac, "pool|H2")))
        elif nh == 1 and A.GetDegree() == 2:
            ket = f"ketone|{edit(key, co)}"
            spec.append(ket)
            R.append(("O:hyd", (ket, "pool|H2"), (alc,)))
        seen = set()
        for fam, prods in ion_moves(key):
            if fam != "beta":
                continue
            ip = next((p for p in prods if "+" in p), None)
            ole = next((p for p in prods if "+" not in p), None)
            if ip is None or ole not in G["alkene"] or has_ring(ip) or "=" in ip:
                continue
            if (ole, ip) in seen:
                continue
            seen.add((ole, ip))
            R.append(("R:beta", (rad,), (f"alkene|{ole}", f"rad|{ip}")))
    return spec, R


def parents(s):
    m = mol(s)
    out = set()
    for b in m.GetBonds():
        if b.GetBondType() == D2:
            rw = Chem.RWMol(m)
            rw.GetBondWithIdx(b.GetIdx()).SetBondType(S1)
            out.add(canon(rw))
    return out


def w_hyd(item):
    cls, s = item
    R = []
    if cls == "alkene":
        for p in parents(s):
            R.append(("H:hyd", (f"alkene|{s}", "pool|H2"), (f"alkane|{p}",)))
    elif cls == "diene":
        for p in parents(s):
            if p in G["alkene"]:
                R.append(("Z:HT", (f"alkene|{p}", "pool|[H]"), (f"diene|{s}",)))
    elif cls == "cycloalkene":
        for p in parents(s):
            if p in G["naphthene"]:
                R.append(("H:hyd", (f"cycloalkene|{s}", "pool|H2"), (f"naphthene|{p}",)))
    elif cls == "cyclodiene":
        for p in parents(s):
            if p in G["cycloalkene"]:
                R.append(("Z:HT", (f"cycloalkene|{p}", "pool|[H]"), (f"cyclodiene|{s}",)))
        m = mol(s)
        ring = m.GetRingInfo().AtomRings()[0]
        if len(ring) == 6:
            rw = Chem.RWMol(m)
            rs = set(ring)
            for b in rw.GetBonds():
                if b.GetBeginAtomIdx() in rs and b.GetEndAtomIdx() in rs:
                    b.SetBondType(AR)
                    b.SetIsAromatic(True)
            for i in rs:
                rw.GetAtomWithIdx(i).SetIsAromatic(True)
            ar = canon(rw)
            if ar in G["aromatic"]:
                R.append(("Z:HT", (f"cyclodiene|{s}", "pool|[H]"), (f"aromatic|{ar}",)))
    return R


def w_arom(s):
    m = mol(s)
    n = nC(m)
    R = set()
    arom, ions_by_c = G["aromatic"], G["alkyl_ions_by_c"]
    me = f"aromatic|{s}"
    for i in reps(m):
        A = m.GetAtomWithIdx(i)
        if not A.GetIsAromatic():
            continue
        if A.GetTotalNumHs() > 0:
            for k in range(1, CMAX - n + 1):
                for ion in ([None] if k == 1 else ions_by_c.get(k, ())):
                    rw = Chem.RWMol(m)
                    if ion is None:
                        j = rw.AddAtom(Chem.Atom(6))
                        rw.AddBond(i, j, S1)
                    else:
                        im = mol(ion)
                        off = rw.GetNumAtoms()
                        for _ in im.GetAtoms():
                            rw.AddAtom(Chem.Atom(6))
                        for bd in im.GetBonds():
                            rw.AddBond(off + bd.GetBeginAtomIdx(), off + bd.GetEndAtomIdx(), bd.GetBondType())
                        rw.AddBond(i, off + cation_idx(im), S1)
                    p = canon(rw)
                    if p in arom:
                        if ion is None:
                            R.add(("A:meth", (me, "core|SMS"), (f"aromatic|{p}", "pool|H+")))
                        else:
                            R.add(("A:FC", (me, f"ion_alkyl|{ion}"), (f"aromatic|{p}", "pool|H+")))
        else:
            for X in A.GetNeighbors():
                if X.GetIsAromatic():
                    continue
                for Rr in A.GetNeighbors():
                    if not Rr.GetIsAromatic() or Rr.GetTotalNumHs() == 0:
                        continue
                    rw = Chem.RWMol(m)
                    rw.RemoveBond(i, X.GetIdx())
                    rw.AddBond(Rr.GetIdx(), X.GetIdx(), S1)
                    p = canon(rw)
                    if p in arom and s < p:
                        R.add(("A:iso", (me,), (f"aromatic|{p}",)))
    # side-chain dehydrocyclisation onto the ortho ring position -> indanes (k = 3) / tetralins (k = 4)
    for i in range(m.GetNumAtoms()):
        A = m.GetAtomWithIdx(i)
        if not A.GetIsAromatic():
            continue
        orthos = [o.GetIdx() for o in A.GetNeighbors() if o.GetIsAromatic() and o.GetTotalNumHs() > 0]
        for X in A.GetNeighbors():
            if X.GetIsAromatic() or X.IsInRing():
                continue
            # every simple path of 3 or 4 side-chain carbons starting at X (branched chains included)
            stack = [[X.GetIdx()]]
            while stack:
                path = stack.pop()
                if len(path) in (3, 4) and m.GetAtomWithIdx(path[-1]).GetTotalNumHs() > 0:
                    for o in orthos:
                        rw = Chem.RWMol(m)
                        rw.AddBond(path[-1], o, S1)
                        p = canon(rw)
                        if p in arom:
                            R.add(("A:cyc", (me,), (f"aromatic|{p}", "pool|H2")))
                if len(path) < 4:
                    for y in m.GetAtomWithIdx(path[-1]).GetNeighbors():
                        if not y.GetIsAromatic() and not y.IsInRing() and y.GetIdx() not in path:
                            stack.append(path + [y.GetIdx()])
    # dehydrogenation of a fused saturated 6-ring (tetralin -> naphthalene type)
    for ring in m.GetRingInfo().AtomRings():
        if len(ring) != 6 or all(m.GetAtomWithIdx(j).GetIsAromatic() for j in ring):
            continue
        if not any(m.GetAtomWithIdx(j).GetIsAromatic() for j in ring):
            continue
        rw = Chem.RWMol(m)
        rs = set(ring)
        for b in rw.GetBonds():
            if b.GetBeginAtomIdx() in rs and b.GetEndAtomIdx() in rs:
                b.SetBondType(AR)
                b.SetIsAromatic(True)
        for j in rs:
            rw.GetAtomWithIdx(j).SetIsAromatic(True)
        p = canon(rw)
        if p in arom and p != s:
            R.add(("Z:HT", (me, "pool|[H]"), (f"aromatic|{p}",)))
    # polymethylbenzenes are hydrocarbon-pool members
    sk = [x for x in m.GetAtoms() if not x.GetIsAromatic()]
    if m.GetRingInfo().NumRings() == 1 and len(sk) >= 2 and all(x.GetDegree() == 1 for x in sk):
        R.add(("A:hcp", (me,), ("core|HCP",)))
    return list(R)


# ------------------------------------------------------------------ condensation / coupling (K: families, R:comb)
def carbonyl_c(m):
    """(carbonyl carbon, its =O) of an aldehyde / ketone / acid"""
    for b in m.GetBonds():
        if b.GetBondType() == D2:
            a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8:
                return a1.GetIdx(), a2.GetIdx()
            if a2.GetAtomicNum() == 6 and a1.GetAtomicNum() == 8:
                return a2.GetIdx(), a1.GetIdx()
    return None


def alpha_sites(m, c):
    """symmetry-distinct sp3 alpha carbons with H next to carbonyl carbon c"""
    rk = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    out, seen = [], set()
    for n in m.GetAtomWithIdx(c).GetNeighbors():
        if n.GetAtomicNum() != 6 or n.GetTotalNumHs() == 0 or any(b.GetBondType() != S1 for b in n.GetBonds()):
            continue
        if rk[n.GetIdx()] not in seen:
            seen.add(rk[n.GetIdx()])
            out.append(n.GetIdx())
    return out


def carbonyl_key(s):
    for k in ("aldehyde", "ketone"):
        if s in G[k]:
            return f"{k}|{s}"
    return None


def w_aldol(item):
    """K:aldol  donor (alpha C-H) + acceptor C=O  <=> beta-hydroxy carbonyl
       K:crot   beta-hydroxy carbonyl <=> enal / enone + H2O
       K:enH    enal / enone + H2 <=> saturated aldehyde / ketone"""
    dcls, ds = item
    md = mol(ds)
    nd = nC(md)
    cD = carbonyl_c(md)
    if cD is None:
        return []
    alphas = alpha_sites(md, cD[0])
    if not alphas:
        return []
    R = set()
    for n in range(1, CMAX - nd + 1):
        for acls, as_ in G["carbonyl_by_c"].get(n, ()):
            ma = mol(as_)
            cA = carbonyl_c(ma)
            if cA is None:
                continue
            off = md.GetNumAtoms()
            for al in alphas:
                rw = Chem.RWMol(Chem.CombineMols(md, ma))
                rw.AddBond(al, off + cA[0], S1)
                rw.GetBondBetweenAtoms(off + cA[0], off + cA[1]).SetBondType(S1)
                add = canon(rw)
                if add is None or not shape_ok(mol(add)):
                    continue
                hk = f"hydroxycarbonyl|{add}"
                R.add(("K:aldol", (f"{dcls}|{ds}", f"{acls}|{as_}"), (hk,)))
                m_add = mol(add)
                # dehydration needs an H left on the alpha carbon: work on the combined graph
                rw2 = Chem.RWMol(rw)
                try:
                    Chem.SanitizeMol(rw2)
                except Exception:
                    continue
                if rw2.GetAtomWithIdx(al).GetTotalNumHs() == 0:
                    continue
                rw2.GetBondBetweenAtoms(al, off + cA[0]).SetBondType(D2)
                rw2.RemoveAtom(off + cA[1])
                for a in rw2.GetAtoms():
                    a.SetNoImplicit(False)
                    a.SetNumExplicitHs(0)
                en = canon(rw2)
                if en is None:
                    continue
                R.add(("K:crot", (hk,), (f"enone|{en}", "pool|H2O")))
                me = mol(en)
                for b in me.GetBonds():                    # the new C=C (not the C=O)
                    if b.GetBondType() == D2 and b.GetBeginAtom().GetAtomicNum() == 6 and b.GetEndAtom().GetAtomicNum() == 6:
                        rw3 = Chem.RWMol(me)
                        rw3.GetBondWithIdx(b.GetIdx()).SetBondType(S1)
                        sk = carbonyl_key(canon(rw3) or "")
                        if sk:
                            R.add(("K:enH", (f"enone|{en}", "pool|H2"), (sk,)))
    return list(R)


def w_keto(a_s):
    """K:keto  R-COOH + R'-COOH <=> R-CO-R' + CO2 + H2O  (needs an alpha H on one acid; formic acid gives aldehydes)"""
    ma = mol(a_s)
    na = nC(ma)

    def parts(m):
        c, o2 = carbonyl_c(m)
        oh = next(n.GetIdx() for n in m.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() == 8 and n.GetIdx() != o2)
        r = next((n.GetIdx() for n in m.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() == 6), None)
        ah = r is not None and m.GetAtomWithIdx(r).GetTotalNumHs() > 0
        return c, o2, oh, r, ah
    cA, oA, ohA, rA, ahA = parts(ma)
    R = set()
    for n in range(1, CMAX - na + 2):
        for b_s in G["acid_by_c"].get(n, ()):
            if b_s < a_s:
                continue
            mb = mol(b_s)
            cB, oB, ohB, rB, ahB = parts(mb)
            if not (ahA or ahB):
                continue
            off = ma.GetNumAtoms()
            rw = Chem.RWMol(Chem.CombineMols(ma, mb))
            if rB is not None:
                rw.AddBond(cA, off + rB, S1)
            for i in sorted([ohA, off + cB, off + oB, off + ohB], reverse=True):
                rw.RemoveAtom(i)
            p = canon(rw)
            k = carbonyl_key(p or "")
            if k is None:
                continue
            R.add(("K:keto", (f"acid|{a_s}", f"acid|{b_s}"), (k, "core|CO2", "pool|H2O")))
    return list(R)


_CONJ = Chem.MolFromSmarts("[C;!a]=[C;!a]-[C;!a]=[C;!a]")


def w_da(d_s):
    """K:DA  conjugated diene + alkene <=> cyclohexene (both orientations of the dienophile)"""
    md = mol(d_s)
    if md.GetRingInfo().NumRings():
        return []
    nd = nC(md)
    sysm = {tuple(t) if t[0] < t[3] else tuple(reversed(t)) for t in md.GetSubstructMatches(_CONJ)}
    if not sysm:
        return []
    R = set()
    for n in range(2, CMAX - nd + 1):
        for e_s in G["alkene_by_c"].get(n, ()):
            me = mol(e_s)
            eb = next(b for b in me.GetBonds() if b.GetBondType() == D2)
            e1, e2 = eb.GetBeginAtomIdx(), eb.GetEndAtomIdx()
            off = md.GetNumAtoms()
            for d1, d2, d3, d4 in sysm:
                for x, y in ((e1, e2), (e2, e1)):
                    rw = Chem.RWMol(Chem.CombineMols(md, me))
                    rw.GetBondBetweenAtoms(d1, d2).SetBondType(S1)
                    rw.GetBondBetweenAtoms(d3, d4).SetBondType(S1)
                    rw.GetBondBetweenAtoms(d2, d3).SetBondType(D2)
                    rw.GetBondBetweenAtoms(off + x, off + y).SetBondType(S1)
                    rw.AddBond(d1, off + x, S1)
                    rw.AddBond(d4, off + y, S1)
                    p = canon(rw)
                    if p in G["cycloalkene"]:
                        R.add(("K:DA", (f"diene|{d_s}", f"alkene|{e_s}"), (f"cycloalkene|{p}",)))
    return list(R)


def w_comb(r_s):
    """R:comb  R. + R'. <=> R-R'  (every unordered pair of alkyl radicals with <= CMAX carbons)"""
    m1 = mol(r_s)
    n1 = nC(m1)
    i1 = cation_idx(m1)
    R = set()
    for n in range(1, CMAX - n1 + 1):
        for r2 in G["rad_by_c"].get(n, ()):
            if r2 < r_s:
                continue
            m2 = mol(r2)
            off = m1.GetNumAtoms()
            rw = Chem.RWMol(Chem.CombineMols(m1, m2))
            rw.AddBond(i1, off + cation_idx(m2), S1)
            for a in rw.GetAtoms():
                a.SetFormalCharge(0)
            p = canon(rw)
            if p in G["alkane"]:
                R.add(("R:comb", (f"rad|{r_s}", f"rad|{r2}"), (f"alkane|{p}",)))
    return list(R)


def w_info(k):
    cls, s = k.split("|", 1)
    if cls in ("core", "pool"):
        return k, (cls, 0 if cls == "pool" else 1, CORE.get(k, s), "")
    m = mol(s)
    sk = skeleton(m)
    return k, (cls, nC(m), s, Chem.MolToSmiles(sk) if sk is not None else "")


# ------------------------------------------------------------------ C1 core
CORE = {  # key: label
    "core|CO2": "CO₂", "core|CO2d": "CO₂^δ−*", "core|HCOO": "HCOO*", "core|H2COO": "H₂COO*", "core|H2CO": "H₂CO*",
    "core|CH3O": "CH₃O*", "core|CO*": "CO*", "core|HCO": "HCO*", "core|HCOH": "HCOH*", "core|CH": "CH*/C*",
    "core|CH2*": "CH₂*", "core|SMS": "CH₃OH/SMS", "core|DME": "DME", "core|HCP": "烃池 HCP",
    "pool|H2": "H₂", "pool|H*": "H*", "pool|H+": "H⁺", "pool|H2O": "H₂O", "pool|[H]": "[H]", "pool|H.": "H•",
}
CORE_RX = [
    ("C1", ("core|CO2", "pool|H*"), ("core|CO2d",)),            # CO2 activation at Vo / basic O
    ("C1", ("core|CO2d", "pool|H*"), ("core|HCOO",)),
    ("C1", ("core|HCOO", "pool|H*"), ("core|H2COO",)),
    ("C1", ("core|H2COO", "pool|H*"), ("core|H2CO", "pool|H2O")),
    ("C1", ("core|H2CO", "pool|H*"), ("core|CH3O",)),
    ("C1", ("core|CH3O", "pool|H*"), ("alcohol|CO",)),
    ("C1", ("core|HCOO", "pool|H*"), ("acid|O=CO",)),
    ("C1", ("core|H2CO",), ("aldehyde|C=O",)),
    ("C1", ("core|CO2", "pool|H2"), ("core|CO*", "pool|H2O")),  # RWGS
    ("C1", ("core|CO*", "pool|H*"), ("core|HCO",)),
    ("C1", ("core|HCO", "pool|H*"), ("core|HCOH",)),
    ("C1", ("core|HCO", "pool|H*"), ("aldehyde|C=O",)),
    ("C1", ("core|HCOH", "pool|H*"), ("core|CH", "pool|H2O")),
    ("C1", ("core|CH", "pool|H*"), ("core|CH2*",)),
    ("C1", ("core|CH2*", "pool|H*"), ("ads|[CH3+]",)),
    ("C1", ("alcohol|CO", "pool|H+"), ("core|SMS", "pool|H2O")),
    ("C1", ("core|SMS", "alcohol|CO"), ("core|DME", "pool|H+")),
    ("C1", ("core|SMS",), ("core|HCP",)),                        # methylation of the hydrocarbon pool
    ("C1", ("core|HCP",), ("alkene|C=C", "pool|H+")),           # side-chain / paring olefin release
    ("C1", ("core|HCP",), ("alkene|C=CC", "pool|H+")),
    ("C1", ("pool|H2",), ("pool|H*",)),
]


# ------------------------------------------------------------------ driver
def pmap(fn, items, chunksize=64):
    items = list(items)
    if not items:
        return []
    with mp.get_context("fork").Pool(PROCS) as pool:
        return list(pool.imap_unordered(fn, items, chunksize=chunksize))


def main():
    global CMAX, FULL, PROCS
    ap = argparse.ArgumentParser()
    ap.add_argument("--cmax", type=int, default=16)
    ap.add_argument("--full-upto", type=int, default=12)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=os.cpu_count() or 4)
    args = ap.parse_args()
    CMAX, FULL, PROCS = args.cmax, args.full_upto, args.procs
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tmp = gzip.open(out / "reactions_raw.tsv.gz", "wt")
    nrx = Counter()

    def emit(rxs):
        for fam, r, p in rxs:
            tmp.write(f"{fam}\t{','.join(r)}\t{','.join(p)}\n")
            nrx[fam] += 1

    # ---- skeleton families (alkanes, naphthenes, aromatics), layer by layer
    S = {k: defaultdict(set) for k in ("alkane", "naphthene", "aromatic")}
    S["alkane"][1] = {"C"}
    S["naphthene"][5] = {"C1CCCC1"}
    S["naphthene"][6] = {"C1CCCCC1"}
    S["aromatic"][6] = {"c1ccccc1"}
    if CMAX >= 9:
        S["aromatic"][9] = {"c1ccc2c(c1)CCC2"}
    if CMAX >= 10:
        S["aromatic"][10] = {"c1ccc2c(c1)CCCC2", "c1ccc2ccccc2c1"}
    for n in range(1, CMAX):
        for k in S:
            if S[k][n]:
                for kids in pmap(w_grow, S[k][n], 32):
                    S[k][n + 1].update(kids)
        log(f"C{n + 1}: alkanes {len(S['alkane'][n + 1])} naphthenes {len(S['naphthene'][n + 1])} aromatics {len(S['aromatic'][n + 1])}")
    for k in S:
        G[k] = set().union(*S[k].values())
    # ---- olefinic families
    G["alkene"] = {c for kids in pmap(w_double, G["alkane"]) for c in kids}
    log("alkenes", len(G["alkene"]))
    G["diene"] = {c for kids in pmap(w_double, G["alkene"]) for c in kids}
    log("dienes", len(G["diene"]))
    G["cycloalkene"] = {c for kids in pmap(w_double, G["naphthene"]) for c in kids}
    cyc2 = {c for kids in pmap(w_double, G["cycloalkene"]) for c in kids}
    G["cyclodiene"] = {s for s in cyc2 if all(b.IsInRing() for b in mol(s).GetBonds() if b.GetBondType() == D2)}
    log("cycloalkenes", len(G["cycloalkene"]), "cyclodienes", len(G["cyclodiene"]))
    # ---- carbenium ions + protonation
    olef = [(k, s) for k in ("alkene", "diene", "cycloalkene", "cyclodiene") for s in G[k]]
    kind = {s: k for k, s in olef}
    G["ions"] = set()
    for s, ions in pmap(w_prot, [s for _, s in olef], 128):
        for c in ions:
            G["ions"].add(c)
        emit([("Z:prot", (f"{kind[s]}|{s}", "pool|H+"), (f"{ion_cls(c)}|{c}",)) for c in ions])
    log("ions", len(G["ions"]))
    abc = defaultdict(list)
    for s in G["ions"]:
        if ion_cls(s) == "ion_alkyl":
            abc[nC(mol(s))].append(s)
    G["alkyl_ions_by_c"] = {k: sorted(v) for k, v in abc.items()}
    # ---- ion rearrangements, beta-scission, methylation, cyclisation
    done = 0
    ions_sorted = sorted(G["ions"])
    for start in range(0, len(ions_sorted), 200000):
        for rxs in pmap(w_ion, ions_sorted[start:start + 200000], 256):
            emit(rxs)
        done = min(len(ions_sorted), start + 200000)
        log(f"  ion moves {done}/{len(ions_sorted)}  reactions so far {sum(nrx.values())}")
    # ---- sites: metal, oxygenates, radicals
    site_spec = set()
    for spec, rxs in pmap(w_site, sorted(G["alkane"]), 16):
        site_spec.update(spec)
        emit(rxs)
    log("site species", len(site_spec))
    # ---- condensation / coupling (needs the site species)
    for k in ("aldehyde", "ketone", "acid", "rad"):
        G[k] = {x.split("|", 1)[1] for x in site_spec if x.startswith(k + "|")}
    by_c = defaultdict(list)
    for k in ("aldehyde", "ketone"):
        for x in G[k]:
            by_c[nC(mol(x))].append((k, x))
    G["carbonyl_by_c"] = {n: sorted(v) for n, v in by_c.items()}
    for k in ("acid", "alkene", "rad"):
        b2 = defaultdict(list)
        for x in G[k]:
            b2[nC(mol(x))].append(x)
        G[k + "_by_c"] = {n: sorted(v) for n, v in b2.items()}
    n0 = sum(nrx.values())
    for rxs in pmap(w_aldol, sorted(by for v in G["carbonyl_by_c"].values() for by in v), 8):
        emit(rxs)
    log("aldol / crotonisation / enone hydrogenation", sum(nrx.values()) - n0)
    n0 = sum(nrx.values())
    for rxs in pmap(w_keto, sorted(G["acid"]), 8):
        emit(rxs)
    log("ketonisation", sum(nrx.values()) - n0)
    n0 = sum(nrx.values())
    for rxs in pmap(w_da, sorted(G["diene"]), 16):
        emit(rxs)
    log("Diels-Alder", sum(nrx.values()) - n0)
    n0 = sum(nrx.values())
    for rxs in pmap(w_comb, sorted(G["rad"]), 16):
        emit(rxs)
    log("radical recombination", sum(nrx.values()) - n0)
    # ---- hydrogenation / hydrogen transfer
    for rxs in pmap(w_hyd, [(k, s) for k in ("alkene", "diene", "cycloalkene", "cyclodiene") for s in G[k]], 256):
        emit(rxs)
    # ---- aromatics
    for rxs in pmap(w_arom, sorted(G["aromatic"]), 4):
        emit(rxs)
    emit(CORE_RX)
    tmp.close()
    log("raw reactions", sum(nrx.values()), dict(nrx))

    # ---- species table with ids
    keys = set(CORE)
    for k in ("alkane", "naphthene", "aromatic", "alkene", "diene", "cycloalkene", "cyclodiene"):
        keys.update(f"{k}|{s}" for s in G[k])
    keys.update(f"{ion_cls(s)}|{s}" for s in G["ions"])
    keys.update(site_spec)
    # make sure every reaction participant exists
    with gzip.open(out / "reactions_raw.tsv.gz", "rt") as fh:
        for line in fh:
            _, r, p = line.rstrip("\n").split("\t")
            for x in r.split(",") + p.split(","):
                if x not in keys:
                    keys.add(x)
    CLASS_ORDER = ["pool", "core", "alcohol", "aldehyde", "ketone", "acid", "hydroxycarbonyl", "enone", "ads", "rad", "alkane", "alkene",
                   "ion_alkyl", "diene", "ion_alkenyl", "naphthene", "cycloalkene", "ion_cyc", "cyclodiene", "aromatic"]
    info = dict(pmap(w_info, sorted(keys), 512))
    order = sorted(keys, key=lambda k: (info[k][1], CLASS_ORDER.index(info[k][0]), info[k][3], info[k][2]))
    idx = {k: i for i, k in enumerate(order)}
    with gzip.open(out / "species.tsv.gz", "wt") as fh:
        fh.write("id\tclass\tC\tsmiles\tskeleton\n")
        for k in order:
            cls, n, s, sk = info[k]
            fh.write(f"{idx[k]}\t{cls}\t{n}\t{s}\t{sk}\n")
    log("species", len(order))

    # ---- reactions table (dedupe, ids, readable equation)
    def show(k):
        cls, s = k.split("|", 1)
        if cls in ("core", "pool"):
            return CORE.get(k, s)
        if cls == "rad":
            return s.replace("+]", "•]")
        if cls == "ads":
            return s.replace("+]", "*]")
        return s
    by_fam, by_c = Counter(), defaultdict(Counter)
    with gzip.open(out / "reactions_raw.tsv.gz", "rt") as fh, gzip.open(out / "reactions.tsv.gz", "wt") as fo:
        fo.write("id\tfamily\treactants\tproducts\tequation\n")
        rid = 0
        for line in fh:
            fam, r, p = line.rstrip("\n").split("\t")
            rs, ps = r.split(","), p.split(",")
            key = (fam, [idx[x] for x in rs], [idx[x] for x in ps])
            cmax = max(info[x][1] for x in rs + ps)
            fo.write(f"{rid}\t{fam}\t{','.join(map(str, key[1]))}\t{','.join(map(str, key[2]))}\t"
                     f"{' + '.join(show(x) for x in rs)} ⇌ {' + '.join(show(x) for x in ps)}\n")
            by_fam[fam] += 1
            by_c[cmax][fam] += 1
            rid += 1
    os.remove(out / "reactions_raw.tsv.gz")
    sp_by = defaultdict(Counter)
    for k in order:
        sp_by[info[k][0]][info[k][1]] += 1
    n_ion = sum(v for c, cnt in sp_by.items() if c.startswith("ion") for v in cnt.values())
    stats = {
        "cmax": CMAX, "full_upto": FULL, "n_species": len(order), "n_reactions": rid,
        "species_by_class_and_C": {c: {str(n): v for n, v in sorted(cnt.items())} for c, cnt in sp_by.items()},
        "reactions_by_family": dict(by_fam.most_common()),
        "reactions_by_max_C": {str(c): dict(v) for c, v in sorted(by_c.items())},
        "pair_families_not_listed": {
            "hydride transfer R+ + R'H (unordered ion pairs)": n_ion * (n_ion - 1) // 2,
        },
        "runtime_s": round(time.time() - T0, 1),
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))
    log("done", len(order), "species", rid, "reactions")


if __name__ == "__main__":
    main()
