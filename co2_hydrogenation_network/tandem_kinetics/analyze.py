"""Kinetics-only analyses of the OX-ZEO model.

1. identifiability: ensemble of fits to the G1 targets from randomized priors -> which
   outputs are fixed by the abstract data and which are not
2. operating sweeps (T, P, GHSV, H2/CO2, oxide fraction) over ensemble members
3. degree of rate control (d ln Y / d ln k) at the G1 central fit
4. pore-size (shape cap) sweep for the C8-C16 cut
5. hydride-transfer ceiling for aromatics selectivity (closed form)
Writes results/*.json and figures/*.png
"""
from __future__ import annotations

import json
import os
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
from scipy.optimize import least_squares

import calibrate as C
from model import Conditions, Model

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
KEYS = ["X_CO2", "S_CO", "S_MeOH+DME", "CH4", "C2-C4=", "C2-C4o", "C5+", "aromatics",
        "olefins_C5+", "paraffins_C5+", "C8-C16", "C5-C11", "C17+", "SAF_usable", "Y_C8-C16", "nA_mean"]
G1 = "G1_GaZrOx_GaZSM5"


def _summ(params, cond, shape=None):
    try:
        s = Model(params, shape=shape).run(cond).summary()
        return {k: float(s[k]) for k in KEYS} | {"jet_cut_comp": {k: float(v) for k, v in s["jet_cut_comp"].items()}}
    except Exception:
        return None


# ------------------------------------------------------------------ 1. identifiability
def _ens_fit(seed):
    rng = np.random.default_rng(seed)
    shift = rng.uniform(-1, 1, len(C.FIT))          # random prior, +-1 decade
    prior = {k: (C.PRIOR[k][0] * 10 ** s, C.PRIOR[k][1]) for k, s in zip(C.FIT, shift)}
    case = C.CASES[G1]

    def resid(x):
        p = {k: (prior[k][0] * 10 ** xi, prior[k][1]) for k, xi in zip(C.FIT, x)}
        try:
            s = Model(p).run(case["cond"]).summary()
        except Exception:
            return np.full(len(case["targets"]) + len(x), 10.0)
        r = [(s[k] - v) / sd for k, (v, sd) in case["targets"].items()]
        return np.array(r + list(C.LAMBDA * np.asarray(x)))

    sol = least_squares(resid, np.zeros(len(C.FIT)), bounds=(-3, 3), diff_step=0.05, max_nfev=60)
    p = {k: (prior[k][0] * 10 ** xi, prior[k][1]) for k, xi in zip(C.FIT, sol.x)}
    s = _summ(p, case["cond"])
    if s is None:
        return None
    tgt_ok = all(abs(s[k] - v) <= 2 * sd for k, (v, sd) in case["targets"].items())
    return {"seed": seed, "ok": tgt_ok, "params": {k: list(v) for k, v in p.items()}, "pred": s}


# ------------------------------------------------------------------ 2. sweeps
SWEEPS = {
    "T_C": [280, 300, 320, 340, 360, 380, 400],
    "P_MPa": [1, 2, 3, 4, 5, 6],
    "GHSV": [1000, 2000, 4000, 8000, 16000],
    "H2_CO2": [2, 3, 4, 5],
    "phi_ox": [0.2, 0.33, 0.5, 0.67, 0.8],
}


def _sweep_member(params):
    base = C.CASES[G1]["cond"]
    out = {}
    for var, vals in SWEEPS.items():
        out[var] = [_summ(params, replace(base, **{var: v})) for v in vals]
    return out


# ------------------------------------------------------------------ 3. DRC
def drc(params, cond, h=0.05):
    base = _summ(params, cond)
    outs = ["X_CO2", "S_CO", "C5+", "aromatics", "C8-C16", "Y_C8-C16", "SAF_usable"]
    res = {}
    for k in C.FIT + ["K_w"]:
        up = dict(params)
        up[k] = (params[k][0] * np.exp(h), params[k][1])
        dn = dict(params)
        dn[k] = (params[k][0] * np.exp(-h), params[k][1])
        su, sd = _summ(up, cond), _summ(dn, cond)
        res[k] = {o: (np.log(max(su[o], 1e-12)) - np.log(max(sd[o], 1e-12))) / (2 * h) for o in outs}
    return base, res


# ------------------------------------------------------------------ 4. pore cap
def _cap(args):
    params, ncap, ncapA = args
    s = _summ(params, C.CASES[G1]["cond"], shape={"n_cap": ncap, "n_capA": ncapA})
    return {"n_cap": ncap, "n_capA": ncapA, "pred": s}


def ht_ceiling():
    nA = np.arange(7, 15)
    nP = np.array([1.5, 2, 2.5, 3, 3.5, 4])
    return {"nA": nA.tolist(), "nP": nP.tolist(),
            "S_A_max": [[float(a / (a + 3 * p)) for a in nA] for p in nP]}


def main():
    os.makedirs(RES, exist_ok=True)
    cal = json.load(open(os.path.join(HERE, "calibration.json")))
    pcen = {k: tuple(v) for k, v in cal[G1]["params"].items()}
    with Pool(os.cpu_count()) as pool:
        ens = [e for e in pool.map(_ens_fit, range(32)) if e]
        good = [e for e in ens if e["ok"]]
        print(f"ensemble: {len(good)}/{len(ens)} fits reproduce the G1 targets within 2 sigma")
        json.dump(ens, open(os.path.join(RES, "ensemble_G1.json"), "w"), indent=1)
        members = [pcen] + [{k: tuple(v) for k, v in e["params"].items()} for e in good[:11]]
        sw = pool.map(_sweep_member, members)
        json.dump({"values": SWEEPS, "members": sw}, open(os.path.join(RES, "sweeps_G1.json"), "w"))
        caps = pool.map(_cap, [(pcen, nc, nca) for nc in (9, 10, 11, 12, 13, 14, 16)
                               for nca in (10, 12)])
        json.dump(caps, open(os.path.join(RES, "pore_cap_G1.json"), "w"), indent=1)
        caps_e = pool.map(_cap, [(m, nc, 10 if nc <= 11 else 12) for m in members[:8]
                                 for nc in (9, 11, 13, 16)])
        json.dump(caps_e, open(os.path.join(RES, "pore_cap_ensemble.json"), "w"), indent=1)
    base, d = drc(pcen, C.CASES[G1]["cond"])
    json.dump({"base": base, "drc": d}, open(os.path.join(RES, "drc_G1.json"), "w"), indent=1)
    json.dump(ht_ceiling(), open(os.path.join(RES, "ht_ceiling.json"), "w"), indent=1)
    # print compact summaries
    for k in KEYS:
        v = np.array([e["pred"][k] for e in good])
        print(f"{k:14s} min {v.min():.3f}  median {np.median(v):.3f}  max {v.max():.3f}")


if __name__ == "__main__":
    main()
