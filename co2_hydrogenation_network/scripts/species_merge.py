"""Species-level merge of a network_gen.py network: one node per chemical species, no ions.

    python species_merge.py ../data/net_C16_full  ../data/net_C16_species        # most complete
    python species_merge.py ../data/net_C16_shape ../data/net_C16_shape_species

What is merged into one node
  1. R+ (acid site: carbenium ion = R-O-Z alkoxide / ion pair), R* (metal site) and R. (gas radical) that share the
     same skeleton and the same site: they are the same fragment "R" (one H missing at that carbon).  Written as a
     neutral open-shell SMILES.
  2. Resonance forms of allylic / pentadienyl fragments (localised writings of one species).
  3. Adsorbed and gas forms of the same molecule / fragment of the C1 core: H+, H*, H. -> H atom;
     CO2^delta-* -> CO2; H2CO* -> formaldehyde; SMS (surface methoxy CH3-O-Z) -> methyl; CH3O* -> methoxy etc.
No pseudo-species: the generator's [H] (hydrogen-transfer acceptor) and HCP (lumped hydrocarbon pool) are removed.
  Z:HT   A + [H] -> B          becomes  A <=> B + n H2   (n = 1 or 2; the H2 stands for the hydrogen handed to another
                                         olefin in a hydrogen-transfer step)
  HCP    CH3 -> HCP, HCP -> C2H4/C3H6 + H, polymethylbenzene -> HCP   are dropped; the hydrocarbon-pool chemistry is
         written explicitly: ring methylation (A:meth), dealkylation (A:FC reverse), beta-H elimination, plus the new
         family A:scm  side-chain methylation  ArC-H + CH3 <=> ArC-CH3 + H  at every benzylic C-H.

Corrections of shorthand steps of the generator (so that every step is element-balanced):
  M:CO   R + CO -> RCHO (one H short)  becomes  R + CO <=> RC(=O) (acyl, new nodes)  and  RC(=O) + H <=> RCHO (M:acylH)
  C1     H2COO + H -> H2CO + H2O        becomes  H2COO + H <=> H2COOH  and  H2COOH + H <=> H2CO + H2O
  C1     H2 -> H                         becomes  H2 <=> H + H

Reactions are remapped to nodes; equations that become identical (in either direction) are merged into one row that
lists every mechanism family that gives it (acid / metal / radical / oxygenate ...) and how many original elementary
reactions it stands for.  Steps in which one side only adds a species to the other (e.g. CO2 + H -> CO2) are dropped.

Outputs in <out>/
  nodes.tsv.gz      node_id  smiles  label  kind  formula  C  H  O  states  n_forms  skeleton
  aliases.tsv.gz    old_id  old_class  old_key  node_id  form_smiles      (every species of the input, -> its node)
  reactions.tsv.gz  id  reactants  products  families  n_original  equation      (split into parts if large)
  stats.json, README.md
"""
from __future__ import annotations

import gzip
import json
import multiprocessing as mp
import os
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
T0 = time.time()
S1, D2 = Chem.BondType.SINGLE, Chem.BondType.DOUBLE


def log(*a):
    print(f"[{time.time() - T0:7.0f}s]", *a, flush=True)


SITE_CLASSES = {"ion_alkyl": "酸位 R⁺（R–O–Z 烷氧基 / 离子对）", "ion_alkenyl": "酸位 R⁺（R–O–Z 烷氧基 / 离子对）",
                "ion_cyc": "酸位 R⁺（R–O–Z 烷氧基 / 离子对）", "ads": "金属位 R*", "rad": "气相自由基 R•"}
MOL_KIND = {"alkane": "烷烃", "alkene": "烯烃", "diene": "二烯", "naphthene": "环烷", "cycloalkene": "环烯",
            "cyclodiene": "环二烯", "aromatic": "芳烃", "alcohol": "醇", "aldehyde": "醛", "ketone": "酮", "acid": "羧酸"}
# C1 core and pool (keys are the labels used in species.tsv of network_gen.py)
CORE = {
    "CO₂": ("O=C=O", "CO₂", "C1 分子", "气相"),
    "CO₂^δ−*": ("O=C=O", "CO₂", "C1 分子", "活化吸附 CO₂^δ−*"),
    "CO*": ("[C-]#[O+]", "CO", "C1 分子", "吸附 CO*"),
    "H₂CO*": ("C=O", "甲醛", "醛", "吸附 H₂CO*"),
    "DME": ("COC", "二甲醚 DME", "醚", "分子"),
    "HCOO*": ("O=C[O]", "甲酸根基团 HCOO", "C1 碎片", "吸附 HCOO*（甲酸盐）"),
    "H₂COO*": ("[O]C[O]", "二氧亚甲基 H₂COO", "C1 碎片", "吸附 H₂COO*"),
    "CH₃O*": ("C[O]", "甲氧基 CH₃O", "C1 碎片", "吸附 CH₃O*"),
    "HCO*": ("[CH]=O", "甲酰基 HCO", "C1 碎片", "吸附 HCO*"),
    "HCOH*": ("O[CH]", "羟基亚甲基 HCOH", "C1 碎片", "吸附 HCOH*"),
    "CH*/C*": ("[CH]", "次甲基 CH（含 C*）", "C1 碎片", "吸附 CH*/C*"),
    "CH₂*": ("[CH2]", "亚甲基 CH₂", "C1 碎片", "吸附 CH₂*"),
    "CH₃OH/SMS": ("[CH3]", "甲基 CH₃", "R 基团", "酸位 R⁺（R–O–Z 烷氧基 / 离子对）"),   # surface methoxy = CH3 on Z
    "烃池 HCP": ("pseudo:HCP", "", "", ""),
    "H*": ("[H]", "H 原子", "氢", "金属位 H*"),
    "H•": ("[H]", "H 原子", "氢", "气相 H•"),
    "H⁺": ("[H]", "H 原子", "氢", "酸位 H⁺（H–Z）"),
    "H₂": ("[H][H]", "H₂", "氢", "分子"),
    "H₂O": ("O", "H₂O", "水", "分子"),
    "[H]": ("pseudo:HT", "", "", ""),
}


def site_to_radical(smi):
    """'+'-marked site key (R+, R*, R.) -> neutral radical mol with H counts frozen"""
    m = Chem.MolFromSmiles(smi)
    rw = Chem.RWMol(m)
    Chem.Kekulize(rw, clearAromaticFlags=True)
    for a in rw.GetAtoms():
        h = a.GetTotalNumHs()
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(h)
        if a.GetFormalCharge() == 1:
            a.SetFormalCharge(0)
            a.SetNumRadicalElectrons(1)
    return rw


def smiles(rw):
    m = Chem.Mol(rw)
    Chem.SanitizeMol(m)
    s = Chem.MolToSmiles(m)
    m2 = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m2) if m2 is not None else s


def resonance_forms(rw):
    """all localised writings reachable by allylic shifts a.-b=c -> a=b-c. (closure)"""
    seen = {}
    todo = [Chem.RWMol(rw)]
    while todo:
        cur = todo.pop()
        s = smiles(cur)
        if s in seen:
            continue
        seen[s] = True
        for a in cur.GetAtoms():
            if a.GetNumRadicalElectrons() != 1:
                continue
            for bd in a.GetBonds():
                if bd.GetBondType() != S1:
                    continue
                b = bd.GetOtherAtom(a)
                for bd2 in b.GetBonds():
                    if bd2.GetBondType() != D2:
                        continue
                    c = bd2.GetOtherAtom(b)
                    if c.GetIdx() == a.GetIdx():
                        continue
                    nw = Chem.RWMol(cur)
                    nw.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx()).SetBondType(D2)
                    nw.GetBondBetweenAtoms(b.GetIdx(), c.GetIdx()).SetBondType(S1)
                    nw.GetAtomWithIdx(a.GetIdx()).SetNumRadicalElectrons(0)
                    nw.GetAtomWithIdx(c.GetIdx()).SetNumRadicalElectrons(1)
                    todo.append(nw)
    return sorted(seen)


def w_site(smi):
    rw = site_to_radical(smi)
    forms = resonance_forms(rw)
    return smi, smiles(rw), forms[0], len(forms)


def formula(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    c = Counter(a.GetSymbol() for a in m.GetAtoms())
    return c.get("C", 0), c.get("H", 0), c.get("O", 0)


def w_acyl(smi):
    """aldehyde R-CH=O -> acyl R-C(.)=O"""
    m = Chem.MolFromSmiles(smi)
    rw = Chem.RWMol(m)
    for a in rw.GetAtoms():
        h = a.GetTotalNumHs()
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(h)
    for a in rw.GetAtoms():
        if a.GetSymbol() == "C" and a.GetNumExplicitHs() > 0 and any(
                b.GetBondType() == D2 and b.GetOtherAtom(a).GetSymbol() == "O" for b in a.GetBonds()):
            a.SetNumExplicitHs(a.GetNumExplicitHs() - 1)
            a.SetNumRadicalElectrons(1)
            return smi, smiles(rw)
    return smi, None


def w_scm(smi):
    """aromatic: methyl added to every symmetry-distinct benzylic C-H -> canonical product SMILES"""
    m = Chem.MolFromSmiles(smi)
    rank = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    seen, out = set(), []
    for a in m.GetAtoms():
        if a.GetSymbol() != "C" or a.GetIsAromatic() or a.GetTotalNumHs() == 0:
            continue
        if not any(n.GetIsAromatic() for n in a.GetNeighbors()) or rank[a.GetIdx()] in seen:
            continue
        seen.add(rank[a.GetIdx()])
        rw = Chem.RWMol(m)
        j = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(a.GetIdx(), j, S1)
        try:
            out.append(smiles(rw))
        except Exception:
            pass
    return smi, out


def w_canon(smi):
    return smiles(Chem.MolFromSmiles(smi)), smi


def w_formula(smi):
    return smi, formula(smi)


def fstr(c, h, o):
    return "".join(f"{el}{n if n > 1 else ''}" for el, n in (("C", c), ("H", h), ("O", o)) if n)


def r_kind(smi):
    m = Chem.MolFromSmiles(smi)
    if m.GetRingInfo().NumRings():
        return "R 基团 · 环状"
    if any(b.GetBondType() == D2 for b in m.GetBonds()):
        return "R 基团 · 烯基 / 烯丙基"
    return "R 基团 · 烷基"


def main():
    src, out = Path(sys.argv[1]), Path(sys.argv[2])
    procs = int(os.environ.get("PROCS", os.cpu_count() or 4))
    out.mkdir(parents=True, exist_ok=True)
    # ---------------------------------------------------------------- species
    old = []                                     # (old_id, class, smiles/label, skeleton)
    with gzip.open(src / "species.tsv.gz", "rt") as fh:
        next(fh)
        for line in fh:
            i, cls, _c, s, sk = line.rstrip("\n").split("\t")
            old.append((int(i), cls, s, sk))
    log(f"input species {len(old):,}")
    site_keys = sorted({s for _, cls, s, _ in old if cls in SITE_CLASSES})
    conv = {}
    with mp.get_context("fork").Pool(procs) as pool:
        for k, (smi, form, rep, nf) in enumerate(pool.imap_unordered(w_site, site_keys, chunksize=2000)):
            conv[smi] = (form, rep)
            if k % 500000 == 0:
                log(f"  R fragments {k:,}/{len(site_keys):,}")
    log(f"R-type keys {len(site_keys):,} converted")
    node_of_key, node_info = {}, {}              # node key -> dict
    alias = []                                   # (old_id, cls, key, node key, form)
    for i, cls, s, sk in old:
        if cls in ("core", "pool"):
            nk, label, kind, state = CORE[s]
            form = nk
            if nk.startswith("pseudo:"):
                alias.append((i, cls, s, nk, ""))
                continue
        elif cls in SITE_CLASSES:
            form, nk = conv[s]
            label, kind, state = "", None, SITE_CLASSES[cls]
        else:
            nk, form, label, kind, state = s, s, "", MOL_KIND[cls], "分子"
        d = node_info.setdefault(nk, {"label": label, "kind": kind, "states": set(), "forms": set(), "sk": sk})
        d["states"].add(state)
        d["forms"].add(form)
        if label and not d["label"]:
            d["label"] = label
        if kind and d["kind"] is None:
            d["kind"] = kind
        if sk and not d["sk"]:
            d["sk"] = sk
        alias.append((i, cls, s, nk, form))
    # acyl fragment for every aldehyde (CO-insertion product), and H2COOH
    alds = sorted(k for k, d in node_info.items() if d["kind"] == "醛")
    with mp.get_context("fork").Pool(procs) as pool:
        acyl_of_key = dict(pool.imap_unordered(w_acyl, alds, chunksize=2000))
    for ald, ac in acyl_of_key.items():
        if ac is None:
            continue
        d = node_info.setdefault(ac, {"label": "", "kind": "酰基 RC(=O)", "states": set(), "forms": {ac}, "sk": node_info[ald]["sk"]})
        d["states"].add("金属位 RC(=O)*（CO 插入产物）")
    node_info.setdefault("[O]CO", {"label": "H₂COOH（羟基甲氧基）", "kind": "C1 碎片", "states": {"吸附 H₂COOH*"},
                                   "forms": {"[O]CO"}, "sk": ""})
    log(f"nodes {len(node_info):,} (acyl fragments added for {len(alds):,} aldehydes)")
    # order: pool / pseudo / C1 first, then by carbon number, kind, smiles
    real = list(node_info)
    with mp.get_context("fork").Pool(procs) as pool:
        fm = dict(pool.imap_unordered(w_formula, real, chunksize=5000))
    for k in node_info:
        node_info[k]["f"] = fm.get(k, (0, 0, 0))
        if node_info[k]["kind"] in (None, "R 基团"):
            node_info[k]["kind"] = r_kind(k)
    order = sorted(node_info, key=lambda k: (node_info[k]["f"][0], node_info[k]["kind"], len(k), k))
    nid = {k: i for i, k in enumerate(order)}
    nid["pseudo:HCP"], nid["pseudo:HT"] = -2, -3
    with gzip.open(out / "nodes.tsv.gz", "wt", compresslevel=6) as fh:
        fh.write("node_id\tsmiles\tlabel\tkind\tformula\tC\tH\tO\tstates\tn_forms\tskeleton\n")
        for k in order:
            d = node_info[k]
            c, h, o = d["f"]
            fh.write(f"{nid[k]}\t{k}\t{d['label']}\t{d['kind']}\t{fstr(c, h, o)}\t{c}\t{h}\t{o}\t"
                     f"{'；'.join(sorted(d['states']))}\t{len(d['forms'])}\t{d['sk']}\n")
    with gzip.open(out / "aliases.tsv.gz", "wt", compresslevel=6) as fh:
        fh.write("old_id\told_class\told_key\tnode_id\tform_smiles\n")
        for i, cls, s, nk, form in alias:     # pseudo-species: node_id empty (removed)
            fh.write(f"{i}\t{cls}\t{s}\t{nid[nk] if nid[nk] >= 0 else ''}\t{form}\n")
    o2n = np.full(len(old), -1, dtype=np.int64)
    for i, _cls, _s, nk, _f in alias:
        o2n[i] = nid[nk]
    log("nodes / aliases written")
    # node statistics now, then free what the reaction stage does not need
    kinds = Counter(node_info[k]["kind"] for k in order)
    by_kind_c = defaultdict(Counter)
    states = Counter()
    for k in order:
        di = node_info[k]
        by_kind_c[di["kind"]][di["f"][0]] += 1
        st = di["states"]
        if any(x.startswith(("酸位 R", "金属位 R*", "气相自由基")) for x in st):
            states[" + ".join(sorted(x.split(" ")[0] for x in st))] += 1
    old_cls = Counter(cls for _, cls, _, _ in old)
    res_merged = sum(1 for k in order if len(node_info[k]["forms"]) > 1)
    n_old = len(old)
    F = np.array([node_info[k]["f"] for k in order], dtype=np.int64)
    aro = [k for k in order if node_info[k]["kind"] == "芳烃"]
    lab = [node_info[k]["label"] if k in ("[H]", "[H][H]", "O", "O=C=O", "[C-]#[O+]") else k for k in order]
    del alias, old, conv, node_info, site_keys
    import gc
    gc.collect()
    # ---------------------------------------------------------------- reactions
    W = 3                                        # species per side
    rfiles = sorted(src.glob("reactions.part*.tsv.gz")) or [src / "reactions.tsv.gz"]
    fams, fi = [], {}
    rows, blocks = [], []
    n_in = 0
    for f in rfiles:
        with gzip.open(f, "rt") as fh:
            for line in fh:
                p = line.split("\t", 4)
                if p[0] == "id":
                    continue
                fam = p[1]
                if fam not in fi:
                    fi[fam] = len(fams)
                    fams.append(fam)
                r = [int(x) for x in p[2].split(",")]
                q = [int(x) for x in p[3].split(",")]
                rows.append((fi[fam], *(r + [-1] * (W - len(r))), *(q + [-1] * (W - len(q)))))
                n_in += 1
                if len(rows) >= 1_000_000:
                    blocks.append(np.array(rows, dtype=np.int64))
                    rows = []
        log(f"  read {f.name}: total {n_in:,}")
    if rows:
        blocks.append(np.array(rows, dtype=np.int64))
    del rows
    A = np.concatenate(blocks)
    del blocks
    famcol = A[:, 0].copy()
    X = A[:, 1:].copy()
    del A
    m = X >= 0
    X[m] = o2n[X[m]]                             # pseudo-species become -2 (HCP) / -3 (HT acceptor)

    def addfam(name):
        if name not in fi:
            fi[name] = len(fams)
            fams.append(name)
        return fi[name]
    Hn, H2n, H2On, CH3n = nid["[H]"], nid["[H][H]"], nid["O"], nid["[CH3]"]
    extra_f, extra = [], []
    # ---- M:CO -> acyl, + M:acylH
    acyl_arr = np.full(len(order), -1, dtype=np.int64)
    for ald, ac in acyl_of_key.items():
        if ac is not None:
            acyl_arr[nid[ald]] = nid[ac]
    if "M:CO" in fi:
        sel = famcol == fi["M:CO"]
        alds_used = np.unique(X[sel, W])
        X[sel, W] = acyl_arr[X[sel, W]]
        assert (X[sel, W] >= 0).all()
        fa = addfam("M:acylH")
        for a_ in alds_used.tolist():
            extra_f.append(fa)
            extra.append((acyl_arr[a_], Hn, -1, a_, -1, -1))
    # ---- C1 shorthand steps
    if "C1" in fi:
        h2coo, h2cooh, h2co = nid["[O]C[O]"], nid["[O]CO"], nid["C=O"]
        for r in np.flatnonzero(famcol == fi["C1"]).tolist():
            row = X[r].tolist()
            L_, R_ = sorted(x for x in row[:W] if x != -1), sorted(x for x in row[W:] if x != -1)
            if L_ == sorted([h2coo, Hn]) and R_ == sorted([h2co, H2On]):
                X[r] = (h2coo, Hn, -1, h2cooh, -1, -1)
                extra_f.append(fi["C1"])
                extra.append((h2cooh, Hn, -1, h2co, H2On, -1))
            elif L_ == [H2n] and R_ == [Hn]:
                X[r] = (H2n, -1, -1, Hn, Hn, -1)
    # ---- hydrogen transfer: drop the [H] acceptor, release n H2
    ht = (X == -3).any(axis=1)
    n_ht = int(ht.sum())
    for r in np.flatnonzero(ht).tolist():
        row = X[r].tolist()
        L_ = [x for x in row[:W] if x >= 0]
        R_ = [x for x in row[W:] if x >= 0]
        dh = int(F[L_, 1].sum() - F[R_, 1].sum())
        assert dh in (2, 4) and len(L_) == 1 and len(R_) == 1 and "-3" not in str(R_), (row, dh)
        R_ += [H2n] * (dh // 2)
        X[r] = (L_ + [-1] * (W - len(L_))) + (R_ + [-1] * (W - len(R_)))
    # ---- hydrocarbon pool: drop lumped steps, add explicit side-chain methylation
    hcp = (X == -2).any(axis=1)
    n_hcp = int(hcp.sum())
    X, famcol = X[~hcp], famcol[~hcp]
    with mp.get_context("fork").Pool(procs) as pool:
        can = dict(pool.imap_unordered(w_canon, aro, chunksize=500))
        scm = list(pool.imap_unordered(w_scm, aro, chunksize=200))
    fs = addfam("A:scm")
    n_scm = 0
    for smi, prods in scm:
        for pr in prods:
            k = can.get(pr)
            if k is not None:
                extra_f.append(fs)
                extra.append((nid[smi], CH3n, -1, nid[k], Hn, -1))
                n_scm += 1
    if extra:
        X = np.vstack([X, np.array(extra, dtype=np.int64)])
        famcol = np.concatenate([famcol, np.array(extra_f, dtype=famcol.dtype)])
    log(f"corrections: {len(extra) - n_scm:,} extra steps; hydrogen transfer rewritten {n_ht:,}; "
        f"hydrocarbon-pool steps dropped {n_hcp:,}; side-chain methylation added {n_scm:,}")
    assert (X >= -1).all()
    # ---- canonical form: sort each side (-1 last), orient so that the smaller side comes first
    BIG = np.iinfo(np.int64).max
    L = np.sort(np.where(X[:, :W] < 0, BIG, X[:, :W]), axis=1)
    R = np.sort(np.where(X[:, W:] < 0, BIG, X[:, W:]), axis=1)
    del X
    swap = np.zeros(len(L), dtype=bool)
    undecided = np.ones(len(L), dtype=bool)
    for c in range(W):
        swap |= undecided & (L[:, c] > R[:, c])
        undecided &= L[:, c] == R[:, c]
    L2 = np.where(swap[:, None], R, L)
    R2 = np.where(swap[:, None], L, R)
    del L, R
    # drop trivial steps: identical sides, or one side contained in the other (it only adds species)

    def contained(P, Q):                         # multiset P (shorter) inside Q, both sorted, BIG padded
        nP, nQ = (P != BIG).sum(axis=1), (Q != BIG).sum(axis=1)
        res = np.zeros(len(P), dtype=bool)
        one = nP == 1
        res |= one & (nQ > 1) & (Q == P[:, :1]).any(axis=1)
        two = (nP == 2) & (nQ == 3)
        for j in range(W):
            Qr = np.delete(Q, j, axis=1)
            res |= two & (Qr[:, 0] == P[:, 0]) & (Qr[:, 1] == P[:, 1])
        return res
    same = (L2 == R2).all(axis=1)
    drop = same | contained(L2, R2) | contained(R2, L2)
    log(f"reactions in {n_in:,}; dropped trivial {int(drop.sum()):,} (identical sides {int(same.sum()):,})")
    K = np.concatenate([L2, R2], axis=1)[~drop]
    fc = famcol[~drop]
    del L2, R2
    o = np.lexsort(tuple(K[:, c] for c in range(2 * W - 1, -1, -1)))
    K, fc = K[o], fc[o]
    new = np.ones(len(K), dtype=bool)
    new[1:] = (K[1:] != K[:-1]).any(axis=1)
    starts = np.flatnonzero(new)
    mask = np.bitwise_or.reduceat((np.int64(1) << fc), starts)
    nor = np.diff(np.append(starts, len(K)))
    U = K[starts]
    U[U == BIG] = -1
    NR = len(U)
    log(f"unique node-level reactions {NR:,}")
    # element balance
    tot = np.zeros((NR, 3), dtype=np.int64)
    for c in range(2 * W):
        v = U[:, c]
        ok = v >= 0
        tot[ok] += (F[v[ok]] if c < W else -F[v[ok]])
    unbal = (tot != 0).any(axis=1)
    log(f"element-unbalanced reactions: {int(unbal.sum()):,}")
    # write reactions
    by_fam = Counter()
    by_nfam = Counter()
    part, rows_in_part, PART_ROWS = 1, 0, 4_000_000
    big = NR > PART_ROWS
    path = lambda p: out / (f"reactions.part{p:02d}.tsv.gz" if big else "reactions.tsv.gz")  # noqa: E731
    fo = gzip.open(path(part), "wt", compresslevel=6)
    hdr = "id\treactants\tproducts\tfamilies\tn_original\tequation\n"
    fo.write(hdr)
    del K, fc, tot
    famlists = {}
    CH = 500_000
    for rid in range(NR):
        if rid % CH == 0:
            Ul = U[rid:rid + CH].tolist()
            ml = mask[rid:rid + CH].tolist()
            nl = nor[rid:rid + CH].tolist()
        j_ = rid % CH
        row = Ul[j_]
        mk = ml[j_]
        fl = famlists.get(mk)
        if fl is None:
            fl = famlists[mk] = [fams[j] for j in range(len(fams)) if mk >> j & 1]
        for x in fl:
            by_fam[x] += 1
        by_nfam[len(fl)] += 1
        rs = [x for x in row[:W] if x >= 0]
        ps = [x for x in row[W:] if x >= 0]
        fo.write(f"{rid}\t{','.join(map(str, rs))}\t{','.join(map(str, ps))}\t{';'.join(fl)}\t{nl[j_]}\t"
                 f"{' + '.join(lab[x] for x in rs)} ⇌ {' + '.join(lab[x] for x in ps)}\n")
        rows_in_part += 1
        if big and rows_in_part >= PART_ROWS and rid + 1 < NR:
            fo.close()
            part += 1
            rows_in_part = 0
            fo = gzip.open(path(part), "wt", compresslevel=6)
            fo.write(hdr)
    fo.close()
    log("reactions written")
    # ---------------------------------------------------------------- stats
    stats = {
        "source": str(src), "input_species": n_old, "input_reactions": n_in,
        "nodes": len(order), "reactions": NR, "dropped_trivial": int(drop.sum()),
        "element_unbalanced": int(unbal.sum()), "hydrogen_transfer_rewritten": n_ht,
        "hydrocarbon_pool_steps_dropped": n_hcp, "side_chain_methylation_added": n_scm,
        "input_species_by_class": dict(old_cls),
        "nodes_by_kind": dict(kinds.most_common()),
        "nodes_by_kind_and_C": {k: {str(c): v for c, v in sorted(cnt.items())} for k, cnt in by_kind_c.items()},
        "R_nodes_by_site_combination": dict(states.most_common()),
        "nodes_with_several_resonance_forms": res_merged,
        "reactions_by_family": dict(by_fam.most_common()),
        "reactions_by_number_of_families": {str(k): v for k, v in sorted(by_nfam.items())},
        "runtime_s": round(time.time() - T0, 1),
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))
    log(json.dumps({k: stats[k] for k in ("nodes", "reactions", "dropped_trivial", "element_unbalanced",
                                            "nodes_with_several_resonance_forms")}))


if __name__ == "__main__":
    main()
