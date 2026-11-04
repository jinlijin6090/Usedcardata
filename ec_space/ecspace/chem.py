"""Chemistry kernel: exhaustive closed-shell C/H/O generation and per-parent derivation of
radicals, ions, Li+ adducts and elementary reactions.

Closed-shell generation: add one heavy atom (C or O) at a time, bonded to any subset of atoms of
the parent with free valence (bond orders 1-3). Every connected (n+1)-atom graph has a non-cut
vertex whose removal leaves a valid n-atom molecule, so this reaches every molecule; aromatic
parents are expanded from all Kekule structures. Canonical SMILES (RDKit) identify molecules.

Removed only when chemically impossible ("valence" tier):
  * normal valences C 4, O 2 (H implicit); a few charge-separated specials (CO) are added by hand
  * triple bond or cumulated (=C=) centre inside a ring smaller than 8 (cannot close geometrically)
  * O-O bonds: allowed as a single peroxide link, but no O-O-O chains (allow_peroxide=yes);
    or no O-O at all (allow_peroxide=no)
Optional "plausible" tier also removes bicyclobutane-type fused small rings, bridged ring systems
(anti-Bredt) and butatriene-or-longer cumulenes.

Derived species keep their parent recoverable (put the removed H / electron back on the marked atom),
so they are unique per parent and need no global de-duplication. Delocalised radical ions M.+/M.-
are written as "<parent SMILES>|+1" / "|-1".
"""
from __future__ import annotations

from itertools import combinations

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
VAL = {"C": 4, "O": 2}
BT = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
S1, D2 = Chem.BondType.SINGLE, Chem.BondType.DOUBLE


# ------------------------------------------------------------------ helpers
def bond_order(b) -> int:
    t = b.GetBondType()
    return 1 if t == S1 else 2 if t == D2 else 3 if t == Chem.BondType.TRIPLE else 0


def counts(m):
    c = o = 0
    for a in m.GetAtoms():
        s = a.GetSymbol()
        if s == "C":
            c += 1
        elif s == "O":
            o += 1
    return c, o


def kekule_forms(mol):
    if any(a.GetIsAromatic() for a in mol.GetAtoms()):
        forms, seen = [], set()
        for m in Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL):
            if m is None:
                continue
            m = Chem.Mol(m)
            Chem.Kekulize(m, clearAromaticFlags=True)
            key = tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), bond_order(b)) for b in m.GetBonds()))
            if key not in seen:
                seen.add(key)
                forms.append(m)
        if forms:
            return forms
    m = Chem.Mol(mol)
    Chem.Kekulize(m, clearAromaticFlags=True)
    return [m]


def stable_smiles(m):
    """Canonical SMILES made robust by round-tripping until stable (guards against
    aromatic/Kekule ambiguities that can give one molecule two canonical strings)."""
    s = Chem.MolToSmiles(m)
    for _ in range(3):
        m2 = Chem.MolFromSmiles(s)
        if m2 is None:
            return s
        s2 = Chem.MolToSmiles(m2)
        if s2 == s:
            return s
        s = s2
    return s


def smiles_of(rw):
    m = Chem.Mol(rw)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return stable_smiles(m)


class Filters:
    def __init__(self, allow_peroxide=True, strain="valence", max_C=99, max_O=99):
        self.allow_peroxide = allow_peroxide
        self.plausible = strain == "plausible"
        self.max_C, self.max_O = max_C, max_O

    def ok(self, m) -> bool:
        m.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(m)
        ri = m.GetRingInfo()
        for b in m.GetBonds():
            if b.GetBondType() == Chem.BondType.TRIPLE and b.IsInRing() and any(b.IsInRingSize(k) for k in range(3, 8)):
                return False
        for a in m.GetAtoms():
            if a.GetSymbol() == "O":
                on = sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "O")
                if on and not self.allow_peroxide:
                    return False
                if on >= 2:
                    return False                       # O-O-O chain
            nd = sum(1 for b in a.GetBonds() if b.GetBondType() == D2)
            if nd >= 2 and a.IsInRing() and any(a.IsInRingSize(k) for k in range(3, 8)):
                return False
        if self.plausible:
            cum = {a.GetIdx() for a in m.GetAtoms() if sum(1 for b in a.GetBonds() if b.GetBondType() == D2) >= 2}
            if any(b.GetBeginAtomIdx() in cum and b.GetEndAtomIdx() in cum for b in m.GetBonds()):
                return False
            if ri.NumRings() > 1:
                rings = [set(r) for r in ri.AtomRings()]
                small = {}
                for r in rings:
                    if len(r) <= 4:
                        for i in r:
                            small[i] = small.get(i, 0) + 1
                if any(v >= 2 for v in small.values()):
                    return False
                if any(len(r1 & r2) >= 3 for r1, r2 in combinations(rings, 2)):
                    return False
        return True


def children(smi: str, F: Filters) -> list[str]:
    mol = Chem.MolFromSmiles(smi)
    nc, no = counts(mol)
    add = ([("C")] if nc < F.max_C else []) + (["O"] if no < F.max_O else [])
    out = set()
    for km in kekule_forms(mol):
        n = km.GetNumAtoms()
        fv = [a.GetTotalNumHs() for a in km.GetAtoms()]
        sym = list(Chem.CanonicalRankAtoms(km, breakTies=False))
        el = [a.GetSymbol() for a in km.GetAtoms()]
        for X in add:
            v = VAL[X]
            cand = [i for i in range(n) if fv[i] > 0]
            atts = []

            def rec(start, left, cur):
                if cur:
                    atts.append(list(cur))
                for j in range(start, len(cand)):
                    i = cand[j]
                    for o in range(1, min(fv[i], left, 3) + 1):
                        cur.append((i, o))
                        rec(j + 1, left - o, cur)
                        cur.pop()

            rec(0, v, [])
            seen1 = set()
            for att in atts:
                if len(att) == 1:
                    key = (sym[att[0][0]], att[0][1])
                    if key in seen1:
                        continue
                    seen1.add(key)
                if X == "O" and not F.allow_peroxide and any(el[i] == "O" for i, _ in att):
                    continue
                rw = Chem.RWMol(km)
                j = rw.AddAtom(Chem.Atom(X))
                for i, o in att:
                    rw.AddBond(i, j, BT[o])
                if not F.ok(rw):
                    continue
                s = smiles_of(rw)
                if s is not None:
                    out.add(s)
    return list(out)


def formula(m):
    c = h = o = li = 0
    for a in m.GetAtoms():
        s = a.GetSymbol()
        h += a.GetTotalNumHs()
        if s == "C":
            c += 1
        elif s == "O":
            o += 1
        elif s == "Li":
            li += 1
        elif s == "H":
            h += 1
    return c, h, o, li


def frozen(mol):
    """Kekulized copy with H counts frozen as explicit counts (so edits never re-derive H)."""
    m = Chem.RWMol(mol)
    Chem.Kekulize(m, clearAromaticFlags=True)
    for a in m.GetAtoms():
        h = a.GetTotalNumHs()
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(h)
    return m


def canon(rw):
    m = Chem.Mol(rw)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return stable_smiles(m)


def info_of(smi):
    """(charge, multiplicity, C, H, O, Li) of a key."""
    if "|" in smi:
        base, q = smi.split("|")
        m = Chem.MolFromSmiles(base)
        c, h, o, li = formula(m)
        return int(q) + sum(a.GetFormalCharge() for a in m.GetAtoms()), 2, c, h, o, li
    m = Chem.MolFromSmiles(smi)
    c, h, o, li = formula(m)
    q = sum(a.GetFormalCharge() for a in m.GetAtoms())
    ne = sum(a.GetNumRadicalElectrons() for a in m.GetAtoms() if a.GetSymbol() != "Li")   # RDKit books phantom radicals on dative Li+
    return q, ne + 1, c, h, o, li


def centres(smi):
    """(radical centres, anion centres, cation centres [(element, n_neighbours_incl_H)], Li count) of a localised key."""
    m = Chem.MolFromSmiles(smi.split("|")[0])
    r = an = 0
    cat = []
    li = 0
    for a in m.GetAtoms():
        s = a.GetSymbol()
        if s == "Li":
            li += 1
            continue
        r += a.GetNumRadicalElectrons()
        q = a.GetFormalCharge()
        if q < 0:
            an += -q
        elif q > 0:
            nb = sum(1 for n in a.GetNeighbors() if n.GetSymbol() != "Li") + a.GetTotalNumHs()
            cat.append((s, nb))
    return r, an, cat, li


POOL = {"H.": "[H]", "H+": "[H+]", "H-": "[H-]", "e-": "[e-]", "Li+": "[Li+]"}


# ------------------------------------------------------------------ derivation
class Deriver:
    """All derived species and elementary reactions that are 'owned' by one closed-shell parent."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.charges = set(cfg.charges)
        self.fam = set(cfg.families)

    def allowed(self, key):
        """Would this localised species be part of the enumerated space under the current configuration?"""
        cfg = self.cfg
        if info_of(key)[0] not in self.charges:
            return False
        r, an, cat, li = centres(key)
        if li > cfg.li_max:
            return False
        prot = [c for c in cat if c[0] == "O" and c[1] == 3]
        other = [c for c in cat if not (c[0] == "O" and c[1] == 3)]
        if prot and not cfg.protonation:
            return False
        if any(c[0] not in cfg.cation_centers for c in other):
            return False
        ncat = len(cat)
        pattern = (r, an, ncat)
        if pattern == (0, 0, 0):
            return True
        if pattern == (1, 0, 0):
            return cfg.radicals
        if pattern == (0, 1, 0):
            return cfg.anions
        if pattern == (0, 0, 1):
            return True
        P = set(cfg.pairs)
        if pattern == (1, 1, 0):
            return "distonic" in P
        if pattern == (0, 2, 0):
            m = Chem.MolFromSmiles(key.split("|")[0])
            return "dianion" in P and all(a.GetSymbol() == "O" for a in m.GetAtoms() if a.GetFormalCharge() < 0)
        if pattern == (2, 0, 0):
            return "biradical" in P
        if pattern == (0, 1, 1):
            return "zwitterion" in P
        if pattern == (1, 0, 1):
            return "radical_cation_pair" in P
        return False

    # site operations on a frozen kekulized copy --------------------------------
    @staticmethod
    def _op(m, i, kind):
        a = m.GetAtomWithIdx(i)
        h = a.GetNumExplicitHs()
        if kind == "rad":
            if h < 1:
                return False
            a.SetNumExplicitHs(h - 1)
            a.SetNumRadicalElectrons(a.GetNumRadicalElectrons() + 1)
        elif kind == "an":
            if h < 1 or a.GetFormalCharge() != 0:
                return False
            a.SetNumExplicitHs(h - 1)
            a.SetFormalCharge(-1)
        elif kind == "cat":
            if h < 1 or a.GetFormalCharge() != 0:
                return False
            a.SetNumExplicitHs(h - 1)
            a.SetFormalCharge(1)
        elif kind == "prot":
            if a.GetSymbol() != "O" or a.GetFormalCharge() != 0 or a.GetNumRadicalElectrons():
                return False
            a.SetNumExplicitHs(h + 1)
            a.SetFormalCharge(1)
        return True

    def _site_ok(self, el, kind):
        if kind == "cat":
            return el in self.cfg.cation_centers
        if kind == "prot":
            return el == "O"
        return el in ("C", "O")

    def derive(self, psmi, write=True):
        """Return (species: dict key -> type, reactions: list of (family, reactants, products))."""
        cfg = self.cfg
        P = Chem.MolFromSmiles(psmi)
        base = frozen(P)
        n = base.GetNumAtoms()
        el = [a.GetSymbol() for a in base.GetAtoms()]
        rank = list(Chem.CanonicalRankAtoms(P, breakTies=False))
        reps = {}
        for i, r in enumerate(rank):
            reps.setdefault(r, i)
        rep_atoms = sorted(reps.values())
        species = {}          # key -> type
        rx = []               # (family, reactants, products)

        def add_sp(key, typ):
            if key and key not in species:
                q, mult, *_ = info_of(key)
                if q in self.charges:
                    species[key] = typ
                    return True
                return False
            return key in species

        add_sp(psmi, "closed_shell")
        kinds = []
        if cfg.radicals:
            kinds.append("rad")
        if cfg.anions:
            kinds.append("an")
        if cfg.cation_centers:
            kinds.append("cat")
        if cfg.protonation:
            kinds.append("prot")
        TYPE1 = {"rad": "radical", "an": "anion", "cat": "cation", "prot": "protonated"}
        LINK = {"rad": ("h_loss", None, "H."), "an": ("deprotonation", None, "H+"), "cat": ("hydride_loss", None, "H-"),
                "prot": ("protonation", "H+", None)}
        singles = {}          # (kind, atom) -> key
        for kind in kinds:
            for i in rep_atoms:
                if not self._site_ok(el[i], kind):
                    continue
                m = Chem.RWMol(base)
                if not self._op(m, i, kind):
                    continue
                k = canon(m)
                if not k or not add_sp(k, TYPE1[kind]):
                    continue
                singles[(kind, i)] = k
                fam, pin, pout = LINK[kind]
                if fam in self.fam:
                    if pin:
                        rx.append((fam, (psmi, POOL[pin]), (k,)))
                    else:
                        rx.append((fam, (psmi,), (k, POOL[pout])))
        if "site_electron_transfer" in self.fam:
            for (kind, i), k in list(singles.items()):
                if kind == "an" and ("rad", i) in singles:
                    rx.append(("site_electron_transfer", (singles[("rad", i)], POOL["e-"]), (k,)))
                if kind == "cat" and ("rad", i) in singles:
                    rx.append(("site_electron_transfer", (singles[("rad", i)],), (k, POOL["e-"])))
        # radical ions (delocalised)
        if cfg.radical_ions:
            for q, fam, typ in ((-1, "electron_capture", "radical_anion"), (1, "electron_loss", "radical_cation")):
                key = f"{psmi}|{q:+d}"
                if q + 0 in self.charges and add_sp(key, typ) and fam in self.fam:
                    if q < 0:
                        rx.append((fam, (psmi, POOL["e-"]), (key,)))
                    else:
                        rx.append((fam, (psmi,), (key, POOL["e-"])))
        # pairs
        PAIRS = {"distonic": ("rad", "an", "distonic_radical_anion"), "dianion": ("an", "an", "dianion"),
                 "biradical": ("rad", "rad", "biradical"), "zwitterion": ("an", "cat", "zwitterion"),
                 "radical_cation_pair": ("rad", "cat", "distonic_radical_cation")}
        pair_keys = []
        for pname in cfg.pairs:
            if pname not in PAIRS:
                continue
            k1, k2, typ = PAIRS[pname]
            for i in range(n):
                if not self._site_ok(el[i], k1) or (pname == "dianion" and el[i] != "O"):
                    continue
                for j in range(n):
                    if (k1 == k2 and j < i) or (j == i and not (pname == "biradical" and cfg.carbenes)):
                        continue
                    if not self._site_ok(el[j], k2) or (pname == "dianion" and el[j] != "O"):
                        continue
                    m = Chem.RWMol(base)
                    if not self._op(m, i, k1) or not self._op(m, j, k2):
                        continue
                    k = canon(m)
                    if not k or k in species:
                        continue
                    if not add_sp(k, typ):
                        continue
                    pair_keys.append(k)
                    if "pair_steps" in self.fam:
                        # the two one-step predecessors (each a single-site species of this parent)
                        for (ka, ia), (kb, ib) in (((k1, i), (k2, j)), ((k2, j), (k1, i))):
                            mm = Chem.RWMol(base)
                            if not self._op(mm, ia, ka):
                                continue
                            pre = canon(mm)
                            if pre not in species:
                                continue
                            fam, pin, pout = LINK[kb]
                            rx.append(("pair_steps", (pre,) + ((POOL[pin],) if pin else ()), (k,) + ((POOL[pout],) if pout else ())))
        # Li+ adducts on O atoms (dative O->Li)
        if cfg.li_max > 0:
            for key, typ in list(species.items()):
                self._li(key, typ, species, rx)
        # bond cleavage of the parent and of its radical ions
        bonds = []
        seenb = set()
        for b in P.GetBonds():
            if b.GetBondType() != S1:
                continue
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            kb = tuple(sorted((rank[i], rank[j]))) + (b.IsInRing(),)
            if kb in seenb:
                continue
            seenb.add(kb)
            bonds.append((i, j, b.IsInRing()))
        for i, j, ring in bonds:
            if "homolysis" in self.fam:
                self._cleave(base, i, j, ("rad", "rad"), psmi, "homolysis", rx, species, ring)
            if "heterolysis" in self.fam:
                for pol in (("cat", "an"), ("an", "cat")):
                    if self._site_ok(el[i], pol[0]) and self._site_ok(el[j], pol[1]):
                        self._cleave(base, i, j, pol, psmi, "heterolysis", rx, species, ring)
            if cfg.radical_ions and "radical_anion_fragmentation" in self.fam:
                for pol in (("rad", "an"), ("an", "rad")):
                    self._cleave(base, i, j, pol, f"{psmi}|-1", "radical_anion_fragmentation", rx, species, ring)
            if cfg.radical_ions and "radical_cation_fragmentation" in self.fam:
                for pol in (("rad", "cat"), ("cat", "rad")):
                    if self._site_ok(el[i], pol[0]) and self._site_ok(el[j], pol[1]):
                        self._cleave(base, i, j, pol, f"{psmi}|+1", "radical_cation_fragmentation", rx, species, ring)
        # beta-scission of every species of this parent that has exactly one localised radical centre
        if "beta_scission" in self.fam:
            for key in list(species):
                if "|" in key or "Li" in key:
                    continue
                self._beta(key, rx)
        return species, list(dict.fromkeys(rx))

    def _cleave(self, base, i, j, pol, src, fam, rx, species, ring):
        m = Chem.RWMol(base)
        m.RemoveBond(i, j)
        for idx, kind in ((i, pol[0]), (j, pol[1])):
            a = m.GetAtomWithIdx(idx)
            if kind == "rad":
                a.SetNumRadicalElectrons(a.GetNumRadicalElectrons() + 1)
            elif kind == "an":
                a.SetFormalCharge(a.GetFormalCharge() - 1)
            elif kind == "cat":
                a.SetFormalCharge(a.GetFormalCharge() + 1)
        mm = m.GetMol()
        try:
            Chem.SanitizeMol(mm)
        except Exception:
            return
        frags = [stable_smiles(f) for f in Chem.GetMolFrags(mm, asMols=True)]
        if not all(self.allowed(f) for f in frags):
            return
        rx.append((fam, (src,), tuple(frags)))

    def _beta(self, key, rx):
        m0 = Chem.MolFromSmiles(key)
        rad = [a.GetIdx() for a in m0.GetAtoms() if a.GetNumRadicalElectrons() == 1]
        if len(rad) != 1 or sum(a.GetNumRadicalElectrons() for a in m0.GetAtoms()) != 1:
            return
        base = frozen(m0)
        a = rad[0]
        seen = set()
        for B in base.GetAtomWithIdx(a).GetNeighbors():
            b = B.GetIdx()
            if base.GetBondBetweenAtoms(a, b).GetBondType() != S1:
                continue
            for C in B.GetNeighbors():
                c = C.GetIdx()
                if c == a or base.GetBondBetweenAtoms(b, c).GetBondType() != S1:
                    continue
                m = Chem.RWMol(base)
                m.RemoveBond(b, c)
                m.GetBondBetweenAtoms(a, b).SetBondType(D2)
                m.GetAtomWithIdx(a).SetNumRadicalElectrons(0)
                m.GetAtomWithIdx(c).SetNumRadicalElectrons(m.GetAtomWithIdx(c).GetNumRadicalElectrons() + 1)
                mm = m.GetMol()
                try:
                    Chem.SanitizeMol(mm)
                except Exception:
                    continue
                frags = tuple(sorted(stable_smiles(f) for f in Chem.GetMolFrags(mm, asMols=True)))
                if frags in seen:
                    continue
                seen.add(frags)
                if not all(self.allowed(f) for f in frags):
                    continue
                rx.append(("beta_scission", (key,), frags))

    def _li(self, key, typ, species, rx):
        cfg = self.cfg
        delocal = "|" in key
        base_smi, q0 = (key.split("|")[0], int(key.split("|")[1])) if delocal else (key, 0)
        m0 = Chem.MolFromSmiles(base_smi)
        qmol = sum(a.GetFormalCharge() for a in m0.GetAtoms()) + q0
        rank = list(Chem.CanonicalRankAtoms(m0, breakTies=False))
        Os = [a.GetIdx() for a in m0.GetAtoms() if a.GetSymbol() == "O" and a.GetFormalCharge() <= 0]
        done = set()
        for k in range(1, cfg.li_max + 1):
            if qmol + k not in self.charges:
                continue
            for combo in combinations(Os, k):
                sig = tuple(sorted(rank[i] for i in combo))
                if sig in done:
                    continue
                done.add(sig)
                m = Chem.RWMol(m0)
                for i in combo:
                    li = m.AddAtom(Chem.Atom(3))
                    m.GetAtomWithIdx(li).SetFormalCharge(1)
                    m.GetAtomWithIdx(li).SetNoImplicit(True)
                    m.AddBond(i, li, Chem.BondType.DATIVE)
                s = canon(m)
                if not s:
                    continue
                nk = f"{s}|{q0:+d}" if delocal else s
                if nk in species:
                    continue
                species[nk] = f"li{k}_" + typ
                if "li_binding" in self.fam:
                    # one Li at a time: previous adduct (k-1 Li) + Li+ -> this adduct
                    if k == 1:
                        rx.append(("li_binding", (key, POOL["Li+"]), (nk,)))
                    else:
                        m1 = Chem.RWMol(m0)
                        for i in combo[:-1]:
                            li = m1.AddAtom(Chem.Atom(3))
                            m1.GetAtomWithIdx(li).SetFormalCharge(1)
                            m1.GetAtomWithIdx(li).SetNoImplicit(True)
                            m1.AddBond(i, li, Chem.BondType.DATIVE)
                        s1 = canon(m1)
                        if s1:
                            rx.append(("li_binding", ((f"{s1}|{q0:+d}" if delocal else s1), POOL["Li+"]), (nk,)))
