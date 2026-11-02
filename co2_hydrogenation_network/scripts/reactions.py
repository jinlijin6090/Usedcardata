"""Explicit elementary-reaction network over a set of closed-shell species.

Species set S (closed shell) -> monoradicals (all H-removals) -> reactions of the
families below. Reactions are stored reversibly (A + B = C + D counted once) and
deduplicated at species level.

Unimolecular / association families (number grows linearly with |S|):
  recomb      R1* + R2* = R1-R2            (bond fission / radical recombination)
  H_recomb    R*  + H*  = RH
  add_beta    R*  + A=B = R-A-B*           (radical addition / beta-scission, R heavy)
  add_H       H*  + A=B = H-A-B*
  ring_open   cyclic R* = open-chain R*    (intramolecular addition / ring opening)
  co_add      R*  + CO  = R-C*=O
  h_shift     R*  = R'*  (same parent, intramolecular H migration)
  elim4       X-A-B-H = A=B + HX           (concerted 4-centre, incl. dehydration)
  elim11      X-C(=O)-Y = CO + X-Y         (decarbonylation, HCOOH -> CO + H2O)
  taut13      keto-enol-type 1,3-H shift
  diels_alder diene + dienophile = cyclohexene
Bimolecular H/atom-transfer families (grow as |radicals| x |partners|):
  h_abstr     RH + Q* = R* + QH
  dispro      R* + Q* = A=B + QH
  sh2         Q* + A-B = Q-A + B*          (radical substitution)
With `pool` given, Q is restricted to that set of small radicals.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from multiprocessing import Pool

from rdkit import Chem

from chem_space import analyze, counts

H_ATOM, H2, CO = "[H]", "[H][H]", "[C-]#[O+]"
QUADRATIC = ("h_abstr", "dispro", "sh2")


def _an(smi):
    return analyze(smi, want_products=True, want_birad=False)


@lru_cache(maxsize=2_000_000)
def join_radicals(r1: str, r2: str) -> str | None:
    """Bond the radical centres of two monoradicals."""
    if r1 == H_ATOM:
        r1, r2 = r2, r1
    m1 = Chem.MolFromSmiles(r1)
    if r2 == H_ATOM:
        rw = Chem.RWMol(m1)
        for a in rw.GetAtoms():
            if a.GetNumRadicalElectrons():
                a.SetNumRadicalElectrons(0)
                a.SetNumExplicitHs(a.GetNumExplicitHs() + 1)
                break
        m = rw.GetMol()
    else:
        m2 = Chem.MolFromSmiles(r2)
        comb = Chem.RWMol(Chem.CombineMols(m1, m2))
        idx = [a.GetIdx() for a in comb.GetAtoms() if a.GetNumRadicalElectrons()]
        if len(idx) != 2:
            return None
        for i in idx:
            comb.GetAtomWithIdx(i).SetNumRadicalElectrons(0)
        comb.AddBond(idx[0], idx[1], Chem.BondType.SINGLE)
        m = comb.GetMol()
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return Chem.MolToSmiles(m)


def _fix(s: str) -> str:
    return CO if s == "[C]=O" else s


class Network:
    def __init__(self, store_quadratic=True):
        self.rxn = {}          # key -> family
        self.fam = defaultdict(int)
        self.store_q = store_quadratic
        self._qhash = set()    # dedupe of unstored quadratic-family reactions

    def add(self, fam, lhs, rhs):
        lhs = tuple(sorted(_fix(x) for x in lhs))
        rhs = tuple(sorted(_fix(x) for x in rhs))
        if lhs == rhs:
            return
        key = min((lhs, rhs), (rhs, lhs))
        if not self.store_q and fam in QUADRATIC:
            h = hash(key)
            if h not in self._qhash:
                self._qhash.add(h)
                self.fam[fam] += 1
            return
        if key not in self.rxn:
            self.rxn[key] = fam
            self.fam[fam] += 1

    def __len__(self):
        return len(self.rxn) + len(self._qhash)

    def species(self):
        s = set()
        for (l, r) in self.rxn:
            s.update(l)
            s.update(r)
        return s


def build(species, pool=None, quadratic=True, near_hshift=False, procs=4, cmax=12, omax=2,
          extra_species=(H2, "O", CO), store_quadratic=True):
    """Build the network on `species` (closed-shell SMILES). Products of a step must be
    in the species set (or its radicals). pool=None: all radicals are H-transfer partners."""
    S = set(species) | set(extra_species)
    S.discard(H2)
    S.discard(CO)
    with Pool(procs) as p:
        res = p.map(_an, sorted(S), chunksize=32)
    an = {r["smiles"]: r for r in res}
    closed = set(an) | {H2, CO}
    parent = {H_ATOM: H2}
    for m, r in an.items():
        for rad in r["products"]["rad"]:
            parent[rad] = m
    rads = set(parent)
    allsp = closed | rads
    net = Network(store_quadratic)
    ok = allsp.__contains__

    for m, r in an.items():
        pr = r["products"]
        for rad in pr["rad"]:
            net.add("H_recomb", (rad, H_ATOM), (m,))
        for a, b in pr["recomb"]:
            if ok(a) and ok(b):
                net.add("recomb", (a, b), (m,))
        for rad, frs in pr["beta_heavy"]:
            if all(ok(_fix(x)) for x in frs):
                net.add("add_beta", tuple(frs), (rad,))
        for rad, u in pr["beta_H"]:
            if ok(u):
                net.add("add_H", (H_ATOM, u), (rad,))
        for rad, op in pr["ring_open"]:
            if ok(op):
                net.add("ring_open", (rad,), (op,))
        for rad, rr in pr["co_add"]:
            if ok(rr):
                net.add("co_add", (rr, CO), (rad,))
        sites = pr["rad"]
        if not near_hshift:
            for i in range(len(sites)):
                for j in range(i + 1, len(sites)):
                    net.add("h_shift", (sites[i],), (sites[j],))
        for prod in pr["elim4"]:
            if all(ok(x) for x in prod):
                net.add("elim4", (m,), tuple(prod))
        for prod in pr["elim11"]:
            prod = tuple(prod) + (CO,)
            if all(ok(x) for x in prod):
                net.add("elim11", (m,), prod)
        for t in pr["taut13"]:
            if ok(t):
                net.add("taut13", (m,), (t,))
        for prod in pr["rda"]:
            if all(ok(x) for x in prod):
                net.add("diels_alder", tuple(prod), (m,))

    partners = sorted(rads if pool is None else (set(pool) & rads) | {H_ATOM})
    pc = {q: ((0, 0) if q == H_ATOM else counts(Chem.MolFromSmiles(q))) for q in partners}
    if quadratic:
        rad_list = sorted(rads)
        for r1 in rad_list:
            for q in partners:
                if q != r1:
                    net.add("h_abstr", (parent[r1], q), (r1, parent[q]))
        for m, r in an.items():
            for rad, u in r["products"]["beta_H"]:
                if not ok(u):
                    continue
                for q in partners:
                    net.add("dispro", (rad, q), (u, parent[q]))
            for (A, B), (ca, oa) in r["products"]["sub_orient"]:
                for q in partners:
                    cq, oq = pc[q]
                    if cq + ca > cmax or oq + oa > omax:
                        continue
                    qa = join_radicals(q, A)
                    if qa and ok(qa) and ok(B):
                        net.add("sh2", (q, m), (qa, B))
    nbh = sum(1 for r in an.values() for _, u in r["products"]["beta_H"] if ok(u))
    return net, {"closed": len(closed), "radicals": len(rads), "partners": len(partners),
                 "beta_H_total": nbh}


def write_csv(net, path):
    import csv
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["family", "reactants", "products"])
        for (l, r), fam in sorted(net.rxn.items(), key=lambda kv: (kv[1], kv[0])):
            w.writerow([fam, " + ".join(l), " + ".join(r)])
