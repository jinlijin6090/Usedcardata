"""Chemical space of C/H/O species for CO2 hydrogenation networks.

Two jobs:

1. Exhaustive generation of closed-shell C/H/O molecules (neutral, normal valences
   C=4, O=2, H implicit) by adding one heavy atom at a time. Every connected graph
   on n+1 heavy atoms has a non-cut vertex, and deleting it leaves a valid
   n-atom molecule, so "parent + one new atom bonded to any set of parent atoms"
   reaches every molecule. Duplicates are removed by RDKit canonical SMILES.
   Aromatic parents are expanded from all of their Kekule structures.

2. Per-molecule elementary-step analysis (radical sites, bond fissions,
   beta-scissions, H shifts, 4-centre eliminations, keto-enol shifts,
   retro-Diels-Alder). Each channel is identified at species level by the
   canonical SMILES of its products, so two symmetry-equivalent sites count once.

Filters
-------
raw tier ("valence"):   no O-O bond, no triple bond and no cumulated C (=C=) in a
                        ring smaller than 8. These are subgraph-monotone, so they
                        can prune during generation.
plausible tier:         raw + no atom in two rings of size <= 4, no bridged ring
                        systems (two SSSR rings sharing >= 3 atoms; this also
                        removes anti-Bredt alkenes), no two adjacent cumulated
                        centres (butatriene and longer cumulenes).
"""
from __future__ import annotations

from itertools import combinations

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

VAL = {"C": 4, "O": 2}
BT = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}


# --------------------------------------------------------------------------- utils
def bond_order(b) -> int:
    t = b.GetBondType()
    if t == Chem.BondType.SINGLE:
        return 1
    if t == Chem.BondType.DOUBLE:
        return 2
    if t == Chem.BondType.TRIPLE:
        return 3
    raise ValueError("non-integer bond order (kekulize first)")


def counts(mol) -> tuple[int, int]:
    c = o = 0
    for a in mol.GetAtoms():
        s = a.GetSymbol()
        if s == "C":
            c += 1
        elif s == "O":
            o += 1
    return c, o


def kekule_forms(mol):
    """All Kekule structures of an aromatic molecule, one copy otherwise."""
    if any(a.GetIsAromatic() for a in mol.GetAtoms()):
        sup = Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL)
        forms, seen = [], set()
        for m in sup:
            if m is None:
                continue
            m = Chem.Mol(m)
            Chem.Kekulize(m, clearAromaticFlags=True)
            key = tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), bond_order(b))
                               for b in m.GetBonds()))
            if key not in seen:
                seen.add(key)
                forms.append(m)
        if forms:
            return forms
    m = Chem.Mol(mol)
    Chem.Kekulize(m, clearAromaticFlags=True)
    return [m]


def explicit_h(mol):
    """Kekulized copy whose H counts are frozen as explicit counts."""
    m = Chem.RWMol(mol)
    Chem.Kekulize(m, clearAromaticFlags=True)
    for a in m.GetAtoms():
        h = a.GetTotalNumHs()
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(h)
    return m


def smiles_of(rw) -> str | None:
    m = Chem.Mol(rw)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return Chem.MolToSmiles(m)


# ------------------------------------------------------------------------- filters
def _ring_info(m):
    m.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(m)
    return m.GetRingInfo()


def check(m) -> tuple[bool, bool]:
    """(raw_ok, plausible_ok) for a kekulized molecule/radical with integer bonds."""
    ri = _ring_info(m)
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.TRIPLE and b.IsInRing():
            if any(b.IsInRingSize(k) for k in range(3, 8)):
                return False, False
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetSymbol() == "O" and a2.GetSymbol() == "O":
            return False, False
    cumul = set()
    for a in m.GetAtoms():
        nd = sum(1 for b in a.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)
        if nd >= 2:
            cumul.add(a.GetIdx())
            if a.IsInRing() and any(a.IsInRingSize(k) for k in range(3, 8)):
                return False, False
    plaus = True
    for b in m.GetBonds():
        if b.GetBeginAtomIdx() in cumul and b.GetEndAtomIdx() in cumul:
            plaus = False
            break
    if plaus and ri.NumRings() > 1:
        rings = [set(r) for r in ri.AtomRings()]
        small_count = {}
        for r in rings:
            if len(r) <= 4:
                for i in r:
                    small_count[i] = small_count.get(i, 0) + 1
        if any(v >= 2 for v in small_count.values()):
            plaus = False
        else:
            for r1, r2 in combinations(rings, 2):
                if len(r1 & r2) >= 3:
                    plaus = False
                    break
    return True, plaus


# ----------------------------------------------------------------------- generation
def children(smi: str, omax: int, cmax: int, sym_reduce: bool = True) -> set[str]:
    """All molecules obtained from `smi` by adding one C or O atom (raw tier)."""
    mol = Chem.MolFromSmiles(smi)
    nc, no = counts(mol)
    add = []
    if nc < cmax:
        add.append("C")
    if no < omax:
        add.append("O")
    out: set[str] = set()
    for km in kekule_forms(mol):
        n = km.GetNumAtoms()
        fv = [a.GetTotalNumHs() for a in km.GetAtoms()]
        sym = list(Chem.CanonicalRankAtoms(km, breakTies=False)) if sym_reduce else list(range(n))
        sym_list = [a.GetSymbol() for a in km.GetAtoms()]
        for X in add:
            v = VAL[X]
            cand = [i for i in range(n) if fv[i] > 0 and not (X == "O" and sym_list[i] == "O")]
            attachments = []

            def rec(start, left, cur):
                if cur:
                    attachments.append(list(cur))
                for j in range(start, len(cand)):
                    i = cand[j]
                    for o in range(1, min(fv[i], left, 3) + 1):
                        cur.append((i, o))
                        rec(j + 1, left - o, cur)
                        cur.pop()

            rec(0, v, [])
            seen_single = set()
            for att in attachments:
                if len(att) == 1:
                    key = (sym[att[0][0]], att[0][1])
                    if key in seen_single:
                        continue
                    seen_single.add(key)
                rw = Chem.RWMol(km)
                j = rw.AddAtom(Chem.Atom(X))
                for i, o in att:
                    rw.AddBond(i, j, BT[o])
                raw_ok, _ = check(rw)
                if not raw_ok:
                    continue
                s = smiles_of(rw)
                if s is not None:
                    out.add(s)
    return out


# ------------------------------------------------------------- per-molecule analysis
def _set_h(a, h):
    a.SetNoImplicit(True)
    a.SetNumExplicitHs(h)


def radical_from(m, i):
    """Radical obtained by removing one H from heavy atom i of explicit-H mol m."""
    rw = Chem.RWMol(m)
    a = rw.GetAtomWithIdx(i)
    _set_h(a, a.GetNumExplicitHs() - 1)
    a.SetNumRadicalElectrons(a.GetNumRadicalElectrons() + 1)
    return rw


def frag_smiles(rw) -> list[str]:
    m = Chem.Mol(rw)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return []
    return sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True))


def analyze(smi: str, want_products: bool = False, want_birad: bool = True) -> dict:
    """Elementary-step channels of one closed-shell molecule and its monoradicals.

    Returned counts (species-level distinct):
      rad            monoradicals R* (= distinct H-bearing sites)
      recomb         R1-R2 <-> R1* + R2* over non-ring heavy-heavy single bonds
      sub_orient     oriented (A-side, B-side) splits for SH2 substitution, with
                     (C, O) of the A-side fragment -> list in 'sub_hist'
      beta_heavy     R* + A=B <-> R-A-B*  (R heavy, intermolecular)
      beta_H         H* + A=B <-> H-A-B*  (also = beta-H channels for disproportionation)
      ring_open      cyclic radical <-> open-chain radical (intramolecular addition)
      co_add         R* + CO <-> R-C*=O
      hshift         intramolecular H migration pairs (all) / with 1,2..1,5 shifts
      elim4          X-A-B-H <-> A=B + HX  (concerted 4-centre, X = H or heavy group)
      taut13         X-H/A=B 1,3-H shift (keto-enol type), counted per direction
      rda            retro-Diels-Alder / Diels-Alder channels
      birad          biradicals/carbenes (2 H removed) - surface-intermediate proxy
    Each count has a raw and a plausible version (key + '_p').
    """
    mol0 = Chem.MolFromSmiles(smi)
    m = explicit_h(mol0)
    raw_ok, plaus_ok = check(Chem.RWMol(m))
    nc, no = counts(m)
    n = m.GetNumAtoms()
    H = [a.GetNumExplicitHs() for a in m.GetAtoms()]
    res = {"smiles": smi, "C": nc, "O": no, "plaus": plaus_ok}
    prod = {}

    # --- monoradicals
    rads = {}
    for i in range(n):
        if H[i] > 0:
            r = radical_from(m, i)
            s = smiles_of(r)
            if s and s not in rads:
                rads[s] = (i, r)
    res["rad"] = len(rads)
    prod["rad"] = sorted(rads)

    # --- distance matrix for H-shift span
    dm = Chem.GetDistanceMatrix(m)
    sites = [v[0] for v in rads.values()]
    hs_all = len(sites) * (len(sites) - 1) // 2
    hs_near = sum(1 for a, b in combinations(sites, 2) if dm[a][b] <= 4)
    res["hshift"], res["hshift_near"] = hs_all, hs_near

    # --- bond fission / recombination and substitution splits
    recomb = set()
    sub_orient = {}
    for b in m.GetBonds():
        if b.IsInRing() or bond_order(b) != 1:
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        rw = Chem.RWMol(m)
        rw.RemoveBond(i, j)
        for k in (i, j):
            rw.GetAtomWithIdx(k).SetNumRadicalElectrons(1)
        m2 = Chem.Mol(rw)
        try:
            Chem.SanitizeMol(m2)
        except Exception:
            continue
        frags_idx = Chem.GetMolFrags(m2)
        frags = Chem.GetMolFrags(m2, asMols=True)
        fs = {}
        for idxs, f in zip(frags_idx, frags):
            fs[i in idxs] = (Chem.MolToSmiles(f), counts(f))
        (si, ci), (sj, cj) = fs[True], fs[False]
        recomb.add(tuple(sorted((si, sj))))
        sub_orient[(si, sj)] = ci     # attack at A=i side, B=j leaves
        sub_orient[(sj, si)] = cj
    res["recomb"] = len(recomb)
    prod["recomb"] = sorted(recomb)
    hist = {}
    for (a, bsm), c in sub_orient.items():
        hist[c] = hist.get(c, 0) + 1
    res["sub_hist"] = hist
    prod["sub_orient"] = sorted(sub_orient.items())

    # --- beta scission from each monoradical
    beta_heavy, beta_H, ring_open, co_add = {}, {}, {}, {}
    for rs, (bi, r) in rads.items():
        rb = r.GetAtomWithIdx(bi)
        for bond in rb.GetBonds():
            ai = bond.GetOtherAtomIdx(bi)
            o_ab = bond_order(bond)
            if o_ab >= 3:
                continue
            aa = r.GetAtomWithIdx(ai)
            # (i) beta-H: H on a
            if aa.GetNumExplicitHs() > 0:
                rw = Chem.RWMol(r)
                _set_h(rw.GetAtomWithIdx(ai), aa.GetNumExplicitHs() - 1)
                rw.GetAtomWithIdx(bi).SetNumRadicalElectrons(0)
                rw.GetBondBetweenAtoms(ai, bi).SetBondType(BT[o_ab + 1])
                ok = check(rw)
                if ok[0]:
                    s = smiles_of(rw)
                    if s:
                        beta_H.setdefault((rs, s), ok[1])
            # (ii) heavy substituent R on a
            for b2 in aa.GetBonds():
                ri = b2.GetOtherAtomIdx(ai)
                if ri == bi or bond_order(b2) != 1:
                    continue
                rw = Chem.RWMol(r)
                rw.RemoveBond(ai, ri)
                rw.GetAtomWithIdx(ri).SetNumRadicalElectrons(
                    rw.GetAtomWithIdx(ri).GetNumRadicalElectrons() + 1)
                rw.GetAtomWithIdx(bi).SetNumRadicalElectrons(0)
                rw.GetBondBetweenAtoms(ai, bi).SetBondType(BT[o_ab + 1])
                ok = check(rw)
                if not ok[0]:
                    continue
                fr = frag_smiles(rw)
                if not fr:
                    continue
                if b2.IsInRing():
                    ring_open.setdefault((rs, fr[0]), ok[1])
                else:
                    beta_heavy.setdefault((rs, tuple(fr)), ok[1])
        # (iii) acyl-type radical R-C*=O  ->  R* + CO
        if rb.GetSymbol() == "C" and rb.GetNumExplicitHs() == 0 and rb.GetDegree() == 2:
            bonds = list(rb.GetBonds())
            dbl = [x for x in bonds if bond_order(x) == 2 and x.GetOtherAtom(rb).GetSymbol() == "O"]
            sgl = [x for x in bonds if bond_order(x) == 1]
            if len(dbl) == 1 and len(sgl) == 1:
                rw = Chem.RWMol(r)
                ri = sgl[0].GetOtherAtomIdx(bi)
                rw.RemoveBond(bi, ri)
                rw.GetAtomWithIdx(ri).SetNumRadicalElectrons(
                    rw.GetAtomWithIdx(ri).GetNumRadicalElectrons() + 1)
                fr = frag_smiles(rw)
                fr = [x for x in fr if x not in ("[C]=O",)]
                if fr:
                    co_add.setdefault((rs, fr[0]), True)
        elif rb.GetSymbol() == "C" and rb.GetNumExplicitHs() == 1 and rb.GetDegree() == 1:
            # formyl H-C*=O -> H* + CO
            bd = list(rb.GetBonds())[0]
            if bond_order(bd) == 2 and bd.GetOtherAtom(rb).GetSymbol() == "O":
                co_add.setdefault((rs, "[H]"), True)
    for key, d in (("beta_heavy", beta_heavy), ("beta_H", beta_H), ("ring_open", ring_open),
                   ("co_add", co_add)):
        res[key] = len(d)
        res[key + "_p"] = sum(1 for v in d.values() if v) if plaus_ok else 0
        prod[key] = sorted(d)

    # --- 4-centre elimination X-A-B-H -> A=B + HX ; 1,3-H shift X(H)-A=B -> X=A-B(H)
    elim, taut = {}, {}
    for bond in m.GetBonds():
        o = bond_order(bond)
        for A, B in ((bond.GetBeginAtom(), bond.GetEndAtom()), (bond.GetEndAtom(), bond.GetBeginAtom())):
            if o < 3 and B.GetNumExplicitHs() > 0:
                # X = H on A
                xs = [("H", None)] if A.GetNumExplicitHs() > 0 else []
                xs += [(None, x) for x in A.GetBonds()
                       if x.GetIdx() != bond.GetIdx() and bond_order(x) == 1 and not x.IsInRing()]
                for hx, xb in xs:
                    rw = Chem.RWMol(m)
                    bi, ai = B.GetIdx(), A.GetIdx()
                    _set_h(rw.GetAtomWithIdx(bi), B.GetNumExplicitHs() - 1)
                    if hx == "H":
                        _set_h(rw.GetAtomWithIdx(ai), A.GetNumExplicitHs() - 1)
                        rw.GetBondBetweenAtoms(ai, bi).SetBondType(BT[o + 1])
                        ok = check(rw)
                        if ok[0]:
                            s = smiles_of(rw)
                            if s:
                                elim.setdefault((s, "[H][H]"), ok[1])
                    else:
                        xi = xb.GetOtherAtomIdx(ai)
                        rw.RemoveBond(ai, xi)
                        xa = rw.GetAtomWithIdx(xi)
                        _set_h(xa, xa.GetNumExplicitHs() + 1)
                        rw.GetBondBetweenAtoms(ai, bi).SetBondType(BT[o + 1])
                        ok = check(rw)
                        if ok[0]:
                            fr = frag_smiles(rw)
                            if fr:
                                elim.setdefault(tuple(fr), ok[1])
            if o >= 2:
                # 1,3-H shift: H on X (neighbour of A, single bond), A=B -> X=A, B-H
                for xb in A.GetBonds():
                    if xb.GetIdx() == bond.GetIdx() or bond_order(xb) != 1:
                        continue
                    X = xb.GetOtherAtom(A)
                    if X.GetNumExplicitHs() == 0:
                        continue
                    rw = Chem.RWMol(m)
                    _set_h(rw.GetAtomWithIdx(X.GetIdx()), X.GetNumExplicitHs() - 1)
                    _set_h(rw.GetAtomWithIdx(B.GetIdx()), B.GetNumExplicitHs() + 1)
                    rw.GetBondBetweenAtoms(A.GetIdx(), B.GetIdx()).SetBondType(BT[o - 1])
                    rw.GetBondBetweenAtoms(A.GetIdx(), X.GetIdx()).SetBondType(BT[2])
                    ok = check(rw)
                    if ok[0]:
                        s = smiles_of(rw)
                        if s and s != smi:
                            taut.setdefault(s, ok[1])
    # --- 1,1-elimination to CO:  X-C(=O)-Y -> CO + X-Y  (HCOOH -> CO + H2O, RCHO -> RH + CO)
    elim11 = {}
    for a in m.GetAtoms():
        if a.GetSymbol() != "C":
            continue
        dbo = [x for x in a.GetBonds() if bond_order(x) == 2 and x.GetOtherAtom(a).GetSymbol() == "O"
               and x.GetOtherAtom(a).GetDegree() == 1]
        sg = [x for x in a.GetBonds() if bond_order(x) == 1]
        if len(dbo) != 1 or len(sg) + a.GetNumExplicitHs() != 2 or any(x.IsInRing() for x in sg):
            continue
        rw = Chem.RWMol(m)
        oi = dbo[0].GetOtherAtomIdx(a.GetIdx())
        nbrs = [x.GetOtherAtomIdx(a.GetIdx()) for x in sg]
        for k in nbrs:
            rw.RemoveBond(a.GetIdx(), k)
        if len(nbrs) == 2:
            rw.AddBond(nbrs[0], nbrs[1], BT[1])
        elif len(nbrs) == 1:
            x = rw.GetAtomWithIdx(nbrs[0])
            _set_h(x, x.GetNumExplicitHs() + 1)
        if len(nbrs) == 0:
            fr = ["[H][H]"]
            ok = (True, True)
        else:
            for k in sorted((a.GetIdx(), oi), reverse=True):
                rw.RemoveAtom(k)
            ok = check(rw)
            fr = frag_smiles(rw)
        if fr and ok[0]:
            elim11.setdefault(tuple(sorted(fr)), ok[1])
    for key, d in (("elim4", elim), ("taut13", taut), ("elim11", elim11)):
        res[key] = len(d)
        res[key + "_p"] = sum(1 for v in d.values() if v) if plaus_ok else 0
        prod[key] = sorted(d)

    # --- retro-Diels-Alder on six-membered rings containing a C=C
    rda = {}
    ri = m.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) != 6 or any(m.GetAtomWithIdx(x).GetSymbol() != "C" for x in ring):
            continue
        for s0 in range(6):
            c = [ring[(s0 + t) % 6] for t in range(6)]  # c[1]=c[2] is the new pi bond
            b12 = m.GetBondBetweenAtoms(c[1], c[2])
            b34 = m.GetBondBetweenAtoms(c[3], c[4])
            b50 = m.GetBondBetweenAtoms(c[5], c[0])
            b01 = m.GetBondBetweenAtoms(c[0], c[1])
            b23 = m.GetBondBetweenAtoms(c[2], c[3])
            b45 = m.GetBondBetweenAtoms(c[4], c[5])
            if bond_order(b12) != 2 or any(bond_order(x) != 1 for x in (b34, b50, b01, b23)):
                continue
            if bond_order(b45) == 3:
                continue
            rw = Chem.RWMol(m)
            rw.RemoveBond(c[3], c[4])
            rw.RemoveBond(c[5], c[0])
            rw.GetBondBetweenAtoms(c[1], c[2]).SetBondType(BT[1])
            rw.GetBondBetweenAtoms(c[0], c[1]).SetBondType(BT[2])
            rw.GetBondBetweenAtoms(c[2], c[3]).SetBondType(BT[2])
            rw.GetBondBetweenAtoms(c[4], c[5]).SetBondType(BT[bond_order(b45) + 1])
            ok = check(rw)
            if ok[0]:
                fr = frag_smiles(rw)
                if fr:
                    rda.setdefault(tuple(fr), ok[1])
    res["rda"] = len(rda)
    res["rda_p"] = sum(1 for v in rda.values() if v) if plaus_ok else 0
    prod["rda"] = sorted(rda)

    # --- biradicals / carbenes (two H removed)
    if want_birad:
        bir = set()
        for i in range(n):
            if H[i] >= 2:
                rw = Chem.RWMol(m)
                a = rw.GetAtomWithIdx(i)
                _set_h(a, H[i] - 2)
                a.SetNumRadicalElectrons(2)
                s = smiles_of(rw)
                if s:
                    bir.add(s)
        for i, j in combinations([k for k in range(n) if H[k] > 0], 2):
            rw = Chem.RWMol(m)
            for k in (i, j):
                a = rw.GetAtomWithIdx(k)
                _set_h(a, H[k] - 1)
                a.SetNumRadicalElectrons(1)
            s = smiles_of(rw)
            if s:
                bir.add(s)
        res["birad"] = len(bir)
    # radical/recombination counts in the plausible tier simply follow the molecule
    for key in ("rad", "recomb", "hshift", "hshift_near"):
        res[key + "_p"] = res[key] if plaus_ok else 0
    if want_birad:
        res["birad_p"] = res["birad"] if plaus_ok else 0
    if want_products:
        res["products"] = prod
    return res
