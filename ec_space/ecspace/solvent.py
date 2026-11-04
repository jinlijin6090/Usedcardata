"""Solvent-partner reactions: bimolecular steps between every enumerated (localised, Li-free) radical / anion and the
configured solvent molecules (EC, DMC, EA, DME, ... or any SMILES).  Each reaction is owned by the open-shell /
anionic reactant X, so it is written once; written as X + S -> products and, like all reactions here, reversible.

  solvent_h_abstraction   X•  + S-H        ->  X-H + S•           radical H abstraction from the solvent
  solvent_deprotonation   X⁻  + S-H        ->  X-H + S⁻           (C or O) anion takes a proton from the solvent
  solvent_acyl_transfer   X⁻  + S C(=O)-OR ->  X-C(=O)... + RO⁻   nucleophilic acyl substitution (BAc2): transesterification,
                                                                    Claisen-type C-acylation by carbanions; ring opening of EC
  solvent_alkyl_transfer  X⁻  + S R-OR'    ->  X-R + R'O⁻         SN2 at an sp3 carbon bound to O (BAl2): methylation by DMC,
                                                                    ether cleavage of DME, alkyl-O ring opening of EC
Products must belong to the enumerated space (allowed(), <= max_heavy).  When the reverse step would also be
generated (a product is itself a solvent), only one direction is kept; degenerate exchanges are dropped.

Stage runner: run(cfg) -> out/solvent/<TAG>/{reactions,species}/L{n}/part-*.tsv.gz and L{n}.json, TAG = "EC+DMC" etc.
Layers <= derive_max_heavy: every species of the base database is a candidate owner (exact, written).
Larger layers: the same crc32 sample of parents as in derive, re-derived on the fly, counts only.
"""
from __future__ import annotations

import gzip
import json
import multiprocessing as mp
import shutil
import time
import zlib
from collections import Counter
from pathlib import Path

from rdkit import Chem

from . import chem

LIBRARY = {   # abbreviation -> SMILES (battery-electrolyte meaning)
    "EC": "O=C1OCCO1",            # ethylene carbonate
    "PC": "CC1COC(=O)O1",         # propylene carbonate
    "VC": "O=C1OC=CO1",           # vinylene carbonate
    "FEC": None,                  # needs F - not in the C/H/O space
    "DMC": "COC(=O)OC",           # dimethyl carbonate
    "DEC": "CCOC(=O)OCC",         # diethyl carbonate
    "EMC": "CCOC(=O)OC",          # ethyl methyl carbonate
    "EA": "CCOC(C)=O",            # ethyl acetate
    "MA": "COC(C)=O",             # methyl acetate
    "DME": "COCCOC",              # 1,2-dimethoxyethane (monoglyme) - the battery "DME"
    "DIMETHYL_ETHER": "COC",      # dimethyl ether - the CO2-chemistry "DME"
    "THF": "C1CCOC1",
    "DOL": "C1COCO1",             # 1,3-dioxolane
}
NAME_CN = {"EC": "碳酸乙烯酯", "PC": "碳酸丙烯酯", "VC": "碳酸亚乙烯酯", "DMC": "碳酸二甲酯", "DEC": "碳酸二乙酯", "EMC": "碳酸甲乙酯",
           "EA": "乙酸乙酯", "MA": "乙酸甲酯", "DME": "乙二醇二甲醚", "DIMETHYL_ETHER": "二甲醚", "THF": "四氢呋喃", "DOL": "1,3-二氧戊环"}
FAMILIES = ["solvent_h_abstraction", "solvent_deprotonation", "solvent_acyl_transfer", "solvent_alkyl_transfer"]


def resolve(names):
    """[(label, canonical key)] for abbreviations or raw SMILES."""
    out = []
    for n in names:
        u = n.strip().upper()
        if u in LIBRARY:
            if LIBRARY[u] is None:
                raise ValueError(f"{n}: contains elements outside C/H/O, not supported")
            smi, lab = LIBRARY[u], u
        else:
            smi, lab = n.strip(), n.strip()
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise ValueError(f"unknown solvent / bad SMILES: {n}")
        out.append((lab, chem.stable_smiles(m)))
    return out


def tag_of(cfg):
    return "+".join(lab for lab, _ in resolve(cfg.solvents)) if cfg.solvents else ""


def heavy_of(m):
    return sum(1 for a in m.GetAtoms() if a.GetSymbol() in ("C", "O"))


def type_of(key):
    """Species class of a localised key (same names as the derive stage)."""
    if "|" in key:
        return "radical_anion" if key.endswith("-1") else "radical_cation"
    r, an, cat, li = chem.centres(key)
    pat = (r, an, len(cat))
    base = {(0, 0, 0): "closed_shell", (1, 0, 0): "radical", (0, 1, 0): "anion", (1, 1, 0): "distonic_radical_anion",
            (0, 2, 0): "dianion", (2, 0, 0): "biradical", (0, 1, 1): "zwitterion", (1, 0, 1): "distonic_radical_cation"}.get(pat)
    if base is None:
        base = "cation" if pat == (0, 0, 1) else "other"
        if pat == (0, 0, 1) and cat[0][0] == "O" and cat[0][1] == 3:
            base = "protonated"
    return f"li{li}_{base}" if li else base


class SolventChem:
    def __init__(self, cfg):
        self.cfg = cfg
        self.D = chem.Deriver(cfg)
        self.fam = set(cfg.solvent_families)
        self.solvents = resolve(cfg.solvents)
        self.keys = {k for _, k in self.solvents}
        self.info = []
        for name, smi in self.solvents:
            m = Chem.MolFromSmiles(smi)
            base = chem.frozen(m)
            rank = list(Chem.CanonicalRankAtoms(m, breakTies=False))
            reps = {}
            for i, r in enumerate(rank):
                reps.setdefault(r, i)
            rad_sites, an_sites = set(), set()
            for i in sorted(reps.values()):
                for kind, store in (("rad", rad_sites), ("an", an_sites)):
                    mm = Chem.RWMol(base)
                    if chem.Deriver._op(mm, i, kind):
                        k = chem.canon(mm)
                        if k and self.D.allowed(k):
                            store.add(k)
            acyl, alkyl, seen = [], [], set()
            for b in base.GetBonds():
                if b.GetBondType() != chem.BT[1]:
                    continue
                for c, o in ((b.GetBeginAtom(), b.GetEndAtom()), (b.GetEndAtom(), b.GetBeginAtom())):
                    if c.GetSymbol() != "C" or o.GetSymbol() != "O" or o.GetDegree() < 2:
                        continue                                   # leaving group must be an O-R (ether / ester O)
                    key = (rank[c.GetIdx()], rank[o.GetIdx()])
                    if key in seen:
                        continue
                    seen.add(key)
                    carbonyl = any(x.GetBondType() == chem.BT[2] and x.GetOtherAtom(c).GetSymbol() == "O" for x in c.GetBonds())
                    sp3 = all(x.GetBondType() == chem.BT[1] for x in c.GetBonds())
                    if carbonyl:
                        acyl.append((c.GetIdx(), o.GetIdx()))
                    elif sp3:
                        alkyl.append((c.GetIdx(), o.GetIdx()))
            self.info.append({"name": name, "smi": smi, "mol": base, "rad": sorted(rad_sites), "an": sorted(an_sites),
                              "acyl": acyl, "alkyl": alkyl, "heavy": heavy_of(m)})

    def describe(self):
        return [{"name": S["name"], "smiles": S["smi"], "radical_sites": S["rad"], "anion_sites": S["an"],
                 "acyl_sites": len(S["acyl"]), "alkyl_sites": len(S["alkyl"])} for S in self.info]

    def _ok(self, k):
        if not k or not self.D.allowed(k):
            return False
        return heavy_of(Chem.MolFromSmiles(k)) <= self.cfg.max_heavy

    def _keep(self, fam, R, P, out):
        rs, ps = tuple(sorted(R)), tuple(sorted(P))
        if rs == ps:
            return                                                 # degenerate exchange
        if any(p in self.keys for p in P) and rs > ps:
            return                                                 # reverse direction is generated by the other owner
        out.append((fam, R, P))

    def reactions(self, key):
        """Solvent-partner reactions owned by species `key` (localised and Li-free, else none)."""
        if "|" in key or "Li" in key:
            return []
        m0 = Chem.MolFromSmiles(key)
        if m0 is None:
            return []
        base = chem.frozen(m0)
        rad = [a.GetIdx() for a in base.GetAtoms() if a.GetNumRadicalElectrons() > 0]
        an = [a.GetIdx() for a in base.GetAtoms() if a.GetFormalCharge() < 0]
        if not rad and not an:
            return []
        hk = heavy_of(m0)
        out = []

        def with_h(i, kind):
            m = Chem.RWMol(base)
            a = m.GetAtomWithIdx(i)
            a.SetNumExplicitHs(a.GetNumExplicitHs() + 1)
            if kind == "rad":
                a.SetNumRadicalElectrons(a.GetNumRadicalElectrons() - 1)
            else:
                a.SetFormalCharge(a.GetFormalCharge() + 1)
            k = chem.canon(m)
            return k if self._ok(k) else None

        xh_rad = {i: with_h(i, "rad") for i in rad} if "solvent_h_abstraction" in self.fam else {}
        xh_an = {i: with_h(i, "an") for i in an} if "solvent_deprotonation" in self.fam else {}
        for S in self.info:
            for i, xh in xh_rad.items():
                if xh:
                    for s in S["rad"]:
                        self._keep("solvent_h_abstraction", (key, S["smi"]), (xh, s), out)
            for i, xh in xh_an.items():
                if xh:
                    for s in S["an"]:
                        self._keep("solvent_deprotonation", (key, S["smi"]), (xh, s), out)
            if hk + S["heavy"] > self.cfg.max_heavy:
                continue
            for fam, sites in (("solvent_acyl_transfer", S["acyl"]), ("solvent_alkyl_transfer", S["alkyl"])):
                if fam not in self.fam:
                    continue
                for i in an:
                    for c_idx, o_idx in sites:
                        comb = Chem.RWMol(Chem.CombineMols(base, S["mol"]))
                        off = base.GetNumAtoms()
                        nu = comb.GetAtomWithIdx(i)
                        nu.SetFormalCharge(nu.GetFormalCharge() + 1)
                        comb.AddBond(i, off + c_idx, chem.BT[1])
                        comb.RemoveBond(off + c_idx, off + o_idx)
                        lg = comb.GetAtomWithIdx(off + o_idx)
                        lg.SetFormalCharge(lg.GetFormalCharge() - 1)
                        mm = comb.GetMol()
                        try:
                            Chem.SanitizeMol(mm)
                        except Exception:
                            continue
                        frags = tuple(chem.stable_smiles(f) for f in Chem.GetMolFrags(mm, asMols=True))
                        if all(self._ok(f) for f in frags):
                            self._keep(fam, (key, S["smi"]), frags, out)
        return list(dict.fromkeys(out))


# ------------------------------------------------------------------ stage
_SC = None
_CFG = None


def _init(cfg):
    global _SC, _CFG
    _CFG = cfg
    _SC = SolventChem(cfg)


def _task(args):
    n, tid, items, rx_dir, sp_dir, write, rederive = args
    rx_c, sp_c = Counter(), Counter()
    rx_lines, sp_lines, new_sp = [], [], set()
    owners = items
    if rederive:            # items are closed-shell parents: re-derive their species
        D = chem.Deriver(_CFG)
        owners = []
        for p in items:
            try:
                owners.extend(D.derive(p)[0].keys())
            except Exception:
                rx_c["_errors"] += 1
    for k in owners:
        try:
            rx = _SC.reactions(k)
        except Exception:
            rx_c["_errors"] += 1
            continue
        for fam, R, P in rx:
            rx_c[fam] += 1
            if write:
                rx_lines.append(f"{fam}\t{' + '.join(R)}\t{' + '.join(P)}\t{n}")
                for x in P:
                    if x not in new_sp and heavy_of(Chem.MolFromSmiles(x)) > _CFG.derive_max_heavy:
                        new_sp.add(x)
    for x in sorted(new_sp):
        q, mult, c, h, o, li = chem.info_of(x)
        t = type_of(x)
        sp_c[t] += 1
        sp_lines.append(f"{x}\t{t}\t{q}\t{mult}\t{c}\t{h}\t{o}\t{li}\t{c + o}\t")
    if write:
        lvl = _CFG.gzip_level
        with gzip.open(Path(rx_dir) / f"part-{tid:06d}.tsv.gz", "wt", compresslevel=lvl) as fh:
            if rx_lines:
                fh.write("\n".join(rx_lines) + "\n")
        if sp_lines:
            with gzip.open(Path(sp_dir) / f"part-{tid:06d}.tsv.gz", "wt", compresslevel=lvl) as fh:
                fh.write("\n".join(sp_lines) + "\n")
    return len(items), rx_c, sp_c


def run(cfg):
    from .pipeline import MP_START, layer_files, log, read_lines
    if not cfg.solvents:
        log("solvent: no [solvent] solvents configured, stage skipped")
        return
    tag = tag_of(cfg)
    root = cfg.out / "solvent" / tag
    root.mkdir(parents=True, exist_ok=True)
    SC = SolventChem(cfg)
    (root / "SOLVENTS.json").write_text(json.dumps({"tag": tag, "families": sorted(SC.fam), "solvents": SC.describe()},
                                                   indent=1, ensure_ascii=False))
    for S in SC.describe():
        log(f"solvent {S['name']} = {S['smiles']}: H-abstraction sites {len(S['radical_sites'])}, deprotonation sites "
            f"{len(S['anion_sites'])}, acyl sites {S['acyl_sites']}, alkyl sites {S['alkyl_sites']}")
    W = max(1, cfg.workers)
    for n in range(1, cfg.max_heavy + 1):
        dj = root / f"L{n}.json"
        if dj.exists():
            log(f"solvent[{tag}] L{n}: already done")
            continue
        write = n <= cfg.derive_max_heavy
        frac = 1.0 if write else max(1e-6, cfg.count_sample)
        rx_dir, sp_dir = root / "reactions" / f"L{n}", root / "species" / f"L{n}"
        for d in (rx_dir, sp_dir):
            if d.exists():
                shutil.rmtree(d)
            if write:
                d.mkdir(parents=True)
        thr = int(frac * 1_000_000)

        def tasks():
            tid, batch = 0, []
            if write:     # owners = every written species of this layer of the base database
                B = 4000
                for f in sorted((cfg.out / "species" / f"L{n}").glob("part-*.tsv.gz")):
                    for line in read_lines(f):
                        p = line.split("\t", 2)
                        if p[1] == "closed_shell" or p[1].startswith("li") or "|" in p[0]:
                            continue
                        batch.append(p[0])
                        if len(batch) >= B:
                            yield (n, tid, batch, str(rx_dir), str(sp_dir), True, False)
                            tid, batch = tid + 1, []
            else:         # same parent sample as derive, re-derived here
                B = 50
                for f in layer_files(cfg.out, n):
                    for s in read_lines(f):
                        if zlib.crc32(s.encode()) % 1_000_000 >= thr:
                            continue
                        batch.append(s)
                        if len(batch) >= B:
                            yield (n, tid, batch, str(rx_dir), str(sp_dir), False, True)
                            tid, batch = tid + 1, []
            if batch:
                yield (n, tid, batch, str(rx_dir), str(sp_dir), write, not write)

        t1 = time.time()
        rx_c, sp_c = Counter(), Counter()
        done, last = 0, time.time()
        ctx = mp.get_context(MP_START)
        with ctx.Pool(W, initializer=_init, initargs=(cfg,)) as pool:
            for k, a, b in pool.imap_unordered(_task, tasks(), chunksize=1):
                done += k
                rx_c.update(a)
                sp_c.update(b)
                if time.time() - last > 60:
                    last = time.time()
                    log(f"  solvent[{tag}] L{n}: {done:,} owners/parents, reactions {sum(rx_c.values()):,}")
        scale = 1.0 / frac
        res = {"n": n, "tag": tag, "exact": write, "processed": done,
               "reactions": {k: round(v * scale) for k, v in rx_c.items()},
               "new_species": dict(sp_c), "seconds": round(time.time() - t1, 1)}
        dj.write_text(json.dumps(res, indent=1))
        log(f"solvent[{tag}] L{n}: reactions {sum(v for k, v in res['reactions'].items() if not k.startswith('_')):,}"
            f"{'' if write else ' (estimated)'}  ({time.time() - t1:.0f}s)")
