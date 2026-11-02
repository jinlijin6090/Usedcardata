"""Extrapolate exact enumeration statistics to C12 and count reactions per family.

Inputs:  data/stats_CHO_O2.json  (C+O <= 9, O <= 2)
         data/stats_CH.json      (hydrocarbons, C <= 10)
Output:  data/estimate_C12.json, printed tables

Model
-----
N(C,O): exact where enumerated; beyond, the growth ratio r(C) = N(C)/N(C-1) is fitted
linearly in C on the last 4 exact ratios; the central estimate keeps the ratio
rising with that slope from the last exact value. Bounds: low =
last exact ratio held constant, high = the ratio rising 1.5x faster than the fit.
Per-molecule channel averages (radicals/molecule, beta channels/molecule, ...) are
fitted linearly in C on the last 3 exact points. The plausible fraction is fitted
log-linearly.
Reaction counts (reversible, species-level):
  linear families  = sum over species of channel counts
  h_abstr          = Nrad (Nrad - 1) / 2
  dispro           = (sum of beta-H channels) * Nrad
  sh2              = 1/2 sum_{oriented splits A|B} #radicals Q with C(Q)+C(A) <= Cmax, O(Q)+O(A) <= Omax
"""
from __future__ import annotations

import json
import os
from math import comb

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CMAX, OMAX = 12, 2
LIN = ["recomb", "rad", "beta_heavy", "beta_H", "ring_open", "co_add", "hshift", "elim4",
       "elim11", "taut13", "rda"]
CH_KEYS = LIN + ["birad", "hshift_near"]


def load():
    st = {}
    for fn in ("stats_CHO_O2.json", "stats_CH.json", "stats_CHcount.json"):
        p = os.path.join(DATA, fn)
        if not os.path.exists(p):
            continue
        for k, v in json.load(open(p))["stats"].items():
            c, o = map(int, k.split(","))
            if "rad" not in v:      # counted but not analysed
                if (c, o) not in st:
                    st[(c, o)] = {"mol": v["mol"], "_counts_only": True}
                continue
            if (c, o) not in st or st[(c, o)].get("_counts_only"):
                st[(c, o)] = v
    return st


def series(st, o, key):
    xs = sorted(c for (c, oo) in st if oo == o and key in st[(c, oo)] and c >= 1)
    return xs, [st[(c, o)][key] for c in xs]


def extrapolate_counts(xs, ys, cmax, mode="central"):
    xs, ys = list(xs), list(ys)
    r = [ys[i] / ys[i - 1] for i in range(1, len(ys))]
    rc = list(range(xs[1], xs[-1] + 1))
    k = min(4, len(r))
    slope = max(np.polyfit(rc[-k:], r[-k:], 1)[0], 0.0)
    out = dict(zip(xs, ys))
    n = ys[-1]
    for c in range(xs[-1] + 1, cmax + 1):
        if mode == "low":          # growth ratio frozen at the last exact value
            ratio = r[-1]
        elif mode == "high":       # ratio keeps rising, 1.5x the fitted slope
            ratio = r[-1] + 1.5 * slope * (c - xs[-1])
        else:                      # ratio keeps rising with the fitted slope
            ratio = r[-1] + slope * (c - xs[-1])
        n *= ratio
        out[c] = n
    return out


def extrapolate_linear(xs, ys, cmax, k=3):
    out = dict(zip(xs, ys))
    slope, icpt = np.polyfit(xs[-k:], ys[-k:], 1)
    for c in range(xs[-1] + 1, cmax + 1):
        out[c] = slope * c + icpt
    return out


def build_table(st, mode="central", cmax=CMAX, omax=OMAX):
    tab = {}
    for o in range(omax + 1):
        xs, ys = series(st, o, "mol")
        N = extrapolate_counts(xs, ys, cmax, mode)
        xa = [c for c in xs if "rad" in st[(c, o)]]
        fp = [st[(c, o)]["mol_p"] / st[(c, o)]["mol"] for c in xa]
        lp = np.polyfit(xa[-3:], np.log(fp[-3:]), 1)
        per = {}
        for key in CH_KEYS:
            per[key] = extrapolate_linear(xa, [st[(c, o)][key] / st[(c, o)]["mol"] for c in xa], cmax)
            per[key + "_p"] = extrapolate_linear(
                xa, [st[(c, o)][key + "_p"] / max(st[(c, o)]["mol_p"], 1) for c in xa], cmax)
        for c in range(1, cmax + 1):
            exact = (c, o) in st and "rad" in st[(c, o)]
            row = {"exact": (c, o) in st, "mol": N[c]}
            row["frac_p"] = (st[(c, o)]["mol_p"] / st[(c, o)]["mol"]) if exact else float(np.exp(np.polyval(lp, c)))
            row["mol_p"] = row["mol"] * row["frac_p"]
            for key in CH_KEYS:
                if exact:
                    row[key] = st[(c, o)][key]
                    row[key + "_p"] = st[(c, o)][key + "_p"]
                else:
                    row[key] = per[key][c] * row["mol"]
                    row[key + "_p"] = per[key + "_p"][c] * row["mol_p"]
            if exact:
                row["subA"] = {k[5:]: v for k, v in st[(c, o)].items() if k.startswith("subA:")}
                row["subA_p"] = {k[7:]: v for k, v in st[(c, o)].items() if k.startswith("subA_p:")}
            tab[(c, o)] = row
    # water (C=0)
    tab[(0, 1)] = {"exact": True, "mol": 1, "mol_p": 1, "frac_p": 1, **{k: 0 for k in CH_KEYS},
                   **{k + "_p": 0 for k in CH_KEYS}}
    tab[(0, 1)]["rad"] = tab[(0, 1)]["rad_p"] = 1
    return tab


def sh2_count(tab, cmax, omax, tier=""):
    # radicals by (C,O)
    nr = {k: v["rad" + tier] for k, v in tab.items() if k[0] <= cmax}
    nr[(0, 0)] = 1  # H atom

    def nrad_fit(cl, ol):
        return sum(v for (c, o), v in nr.items() if c <= cl and o <= ol)

    total = 0.0
    for (c, o), row in tab.items():
        if c > cmax or c == 0:
            continue
        hist = row.get("subA" + tier)
        if not hist:
            # model: oriented splits = 2*recomb, A-side carbon uniform on 0..c, O shared pro rata
            splits = 2 * row["recomb" + tier]
            hist = {}
            for ca in range(0, c + 1):
                oa = round(o * ca / c) if c else 0
                hist[f"{ca},{oa}"] = hist.get(f"{ca},{oa}", 0) + splits / (c + 1)
        for k, v in hist.items():
            ca, oa = map(int, k.split(","))
            total += v * nrad_fit(cmax - ca, omax - oa)
    return total / 2


def totals(tab, cmax=CMAX, omax=OMAX, tier=""):
    rows = [v for (c, o), v in tab.items() if c <= cmax and o <= omax]
    g = lambda k: sum(r[k + tier] for r in rows)
    n_closed = sum(r["mol" + tier] for r in rows) + 2         # + H2, CO
    n_rad = g("rad") + 1                                       # + H atom
    fam = {
        "recomb (R+R'=RR')": g("recomb"),
        "H_recomb (R+H=RH)": g("rad"),
        "add_beta (R+A=B)": g("beta_heavy"),
        "add_H (H+A=B)": g("beta_H"),
        "ring_open/cyclisation": g("ring_open"),
        "co_add (R+CO)": g("co_add"),
        "h_shift (all)": g("hshift"),
        "elim4 (A=B+HX)": g("elim4"),
        "elim11 (CO extrusion)": g("elim11"),
        "taut13 (keto-enol)": g("taut13") / 2,
        "diels_alder": g("rda"),
    }
    lin = sum(fam.values())
    q = {
        "h_abstr (RH+Q=R+QH)": n_rad * (n_rad - 1) / 2,
        "dispro (R+Q=A=B+QH)": g("beta_H") * n_rad,
        "sh2 (Q+AB=QA+B)": sh2_count(tab, cmax, omax, tier),
    }
    return {"closed_shell": n_closed, "radicals": n_rad, "species": n_closed + n_rad,
            "biradicals": g("birad"), "hshift_near (1,2-1,5)": g("hshift_near"),
            "linear_families": fam, "linear_total": lin, "quadratic_families": q,
            "quadratic_total": sum(q.values()), "reactions_total": lin + sum(q.values())}


def stereo_factor(sample=1500, seed=1):
    """Average number of stereoisomers per constitutional isomer, sampled per heavy-atom level."""
    import gzip
    import random
    from rdkit import Chem
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
    opt = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=4096)
    out = {}
    random.seed(seed)
    for h in range(4, 9):
        p = os.path.join(DATA, "species", f"CHO_O2_h{h}.smi.gz")
        if not os.path.exists(p):
            continue
        L = gzip.open(p, "rt").read().split()
        L = random.sample(L, min(sample, len(L)))
        tot = 0
        for s in L:
            m = Chem.MolFromSmiles(s)
            tot += sum(1 for _ in EnumerateStereoisomers(m, options=opt))
        out[h] = tot / len(L)
    return out


def main():
    st = load()
    res = {}
    for mode in ("low", "central", "high"):
        tab = build_table(st, mode)
        res[mode] = {"full": totals(tab), "plausible": totals(tab, tier="_p")}
        if mode == "central":
            res["by_cmax"] = {cm: {"full": totals(tab, cm), "plausible": totals(tab, cm, tier="_p")}
                              for cm in range(1, CMAX + 1)}
            res["by_carbon"] = {f"{c},{o}": {k: v for k, v in row.items() if not k.startswith("subA")}
                                for (c, o), row in sorted(tab.items())}
    res["stereo_factor_by_heavy_atoms"] = stereo_factor()
    with open(os.path.join(DATA, "estimate_C12.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    for mode in ("low", "central", "high"):
        for tier in ("full", "plausible"):
            t = res[mode][tier]
            print(f"{mode:8s} {tier:9s} closed {t['closed_shell']:.3g}  rad {t['radicals']:.3g}  "
                  f"lin {t['linear_total']:.3g}  quad {t['quadratic_total']:.3g}  total {t['reactions_total']:.3g}")
    print("cmax  closed  radicals  linear  quadratic  total")
    for cm, t in res["by_cmax"].items():
        t = t["full"]
        print(cm, f"{t['closed_shell']:.3g} {t['radicals']:.3g} {t['linear_total']:.3g} "
                  f"{t['quadratic_total']:.3g} {t['reactions_total']:.3g}")
    print("stereo", res["stereo_factor_by_heavy_atoms"])


if __name__ == "__main__":
    main()
