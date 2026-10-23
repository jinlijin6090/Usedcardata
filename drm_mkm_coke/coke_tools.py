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


# ---------------------------------------------------------------------------
# Robust steady state for strongly perturbed parameter sets (Sobol samples):
# 1) the engine's solver, 2) site-normalised ODE, 3) Newton root-finding from
# the base-case steady state with the site balances enforced.
# ---------------------------------------------------------------------------
from scipy.integrate import solve_ivp as _solve_ivp
from scipy.optimize import root as _root
from mkmc.sites import calculate_sites as _calc_sites
from mkmc.config import NA as _NA

def _site_norm(pp):
    s = _calc_sites(*pp.as_site_args())
    n_metal = net.N_METAL
    return np.array([s.NSm / _NA] * n_metal + [s.NSs / _NA] * (len(net.SPECIES) - n_metal))

def _theta0():
    y = np.zeros(len(net.SPECIES)); y[net.SPECIES.index("*")] = 1; y[net.SPECIES.index("Ov")] = 0.7; y[net.SPECIES.index("v")] = 0.3
    return y

def steady_state_robust(rc, pp, theta_guess=None):
    """Steady state reached from the clean surface: engine solver first, then the
    site-normalised ODE (Radau/BDF, tiny first step).  No root-finding fallback,
    because the network can have several steady states and only the one reached
    from the fresh catalyst is meaningful."""
    from mkmc.solver import steady_state as _ss
    from mkmc.rates import make_geom
    geom = make_geom(pp)
    norm = _site_norm(pp)
    rhs = lambda t, y: (net.STOICH @ net_rates(rate_vector(np.clip(y, 0, None), rc.k, pp.p, geom))) / norm
    y = _theta0(); t_end = 1e3
    def imbalance(y):
        r = rate_vector(np.clip(y, 0, None), rc.k, pp.p, geom)
        return np.max(np.abs(net.STOICH @ net_rates(r))) / max(np.max(np.abs(r)), 1e-300)
    try:
        t, th, g = _ss(rc.k, pp)
        if imbalance(th[-1]) < 1e-5:                    # net production of every species << largest rate
            return th[-1], geom
        y = np.clip(th[-1], 0, None)                    # otherwise continue from there in real time
    except Exception:
        pass
    for attempt in range(6):
        res = _solve_ivp(rhs, (0, t_end), y, method="BDF", rtol=1e-7, atol=1e-12)
        if not res.success or res.y.shape[1] < 2:
            res = _solve_ivp(rhs, (0, t_end), y, method="Radau", rtol=1e-7, atol=1e-12, first_step=1e-12)
        if not res.success or res.y.shape[1] < 2:
            raise RuntimeError("normalised ODE failed")
        y = np.clip(res.y[:, -1], 0, None)
        if imbalance(y) < 1e-5:
            return y, geom
        t_end *= 100
    raise RuntimeError("robust steady state did not converge")

def solve_robust(rc, pp, theta_guess=None):
    ss, geom = steady_state_robust(rc, pp, theta_guess)
    r = rate_vector(ss, rc.k, pp.p, geom); rn = net_rates(r)
    return dict(r_DRM=rn[5], r_CH4=rn[0], r_coke=2.0 * rn[19], S_coke=2.0 * rn[19] / rn[0] if rn[0] > 0 else np.nan,
                th_C=ss[net.SPECIES.index("C*")], th_O=ss[net.SPECIES.index("O*")], _ss=ss)
