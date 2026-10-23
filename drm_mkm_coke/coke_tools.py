"""Shared helpers for the coking analyses."""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mkmc import *
from mkmc import network as net
from mkmc.network import REACTIONS
from mkmc.solver import steady_state
from mkmc.rates import rate_vector, net_rates

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 130, "savefig.bbox": "tight"})

_CACHE = {}

def solve(rc, pp, theta0=None):
    """Steady state -> dict of coking-relevant metrics. rc: RateConstants (40 k)."""
    t, th, geom = steady_state(rc.k, pp, theta0=theta0)
    ss = th[-1]
    r = rate_vector(ss, rc.k, pp.p, geom)
    rn = net_rates(r)
    th_d = {s: ss[i] for i, s in enumerate(net.SPECIES)}
    out = dict(
        r_DRM=rn[5],                 # net H2 formation (R6)  [mol/gcat/s]
        r_CH4=rn[0],                 # CH4 consumption (R1)
        r_CO2=rn[11] + rn[13],       # CO2 consumption (R12 + R14)
        r_CO=rn[9] + rn[12],         # CO formation (R10 + R13)
        r_H2O=rn[16],
        r_coke=2.0 * rn[19],         # carbon atoms deposited (R20 x 2)  [mol C/gcat/s]
        S_coke=2.0 * rn[19] / rn[0] if rn[0] > 0 else np.nan,   # fraction of CH4 carbon -> coke
        r_cc_f=r[34], r_cc_r=r[35], r_cc_net=rn[17],
        r_c2ox=rn[18],
        r_Cform=r[8], r_Cox=r[16], r_Otrans=rn[10],
        ratio_Cform_Otrans=r[8] / r[20] if r[20] > 0 else np.nan,
        H2_CO=rn[5] / (rn[9] + rn[12]) if (rn[9] + rn[12]) != 0 else np.nan,
    )
    out.update({f"th_{s}": v for s, v in th_d.items()})
    out["_ss"] = ss; out["_r"] = r; out["_rn"] = rn
    return out

def objective(rc, pp, key, theta0=None):
    return solve(rc, pp, theta0)[key]

def drc(rc, pp, key="r_DRM", dk=0.01, reactions=None):
    """Campbell degree of rate control of objective `key` w.r.t. each reaction
    (forward k perturbed, reverse rescaled to hold Keq)."""
    base = objective(rc, pp, key)
    idx = reactions if reactions is not None else range(len(REACTIONS))
    X = {}
    for i in idx:
        f = 2 * i
        if rc.k[f] == 0:
            X[REACTIONS[i]] = 0.0; continue
        kp = rc.k.copy(); d = dk * kp[f]; kp[f] += d
        if np.isfinite(rc.Keq[i]) and rc.Keq[i] > 0:
            kp[f + 1] = kp[f] / rc.Keq[i]
        rc2 = RateConstants(k=kp, Keq=rc.Keq, Ea=rc.Ea, k0=rc.k0)
        rp = objective(rc2, pp, key)
        X[REACTIONS[i]] = (rc.k[f] / base) * (rp - base) / d
    return pd.Series(X)

def pp_with(**kw):
    """ProcessParams with partial-pressure keywords pCH4,pCO2,pH2,pCO,pH2O."""
    pp = ProcessParams()
    p = pp.p.copy()
    for i, name in enumerate(["pCH4", "pCO2", "pH2", "pCO", "pH2O"]):
        if name in kw: p[i] = kw.pop(name)
    pp.p = p
    for k, v in kw.items(): setattr(pp, k, v)
    return pp

def sweep(var_values, make_pp, make_rc, keys=None):
    rows = []
    for v in var_values:
        pp = make_pp(v); rc = make_rc(v)
        try:
            o = solve(rc, pp)
        except Exception as e:
            o = {"error": str(e)}
        o = {k: val for k, val in o.items() if not k.startswith("_")}
        o["x"] = v
        rows.append(o)
    return pd.DataFrame(rows)
