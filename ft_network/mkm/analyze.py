#!/usr/bin/env python3
"""对拟合结果做产物分布与速率控制度 (DRC) 分析。"""
import json
import math
from pathlib import Path

import numpy as np

from co_fts_mkm import DEFAULT, CoMKM
from params_lit import P_REF, R_REF, T_REF, dft_params

HERE = Path(__file__).resolve().parent
import os
BEST = os.environ.get("BEST", "D")
STEPS = ["dir", "hco", "hcod", "ch", "ch2", "ch3", "ch4", "oh", "h2o", "ohoh", "co2", "grow", "par", "beta"]


def model(variant):
    p = dict(DEFAULT, **dft_params())
    p["nsite_dir"] = 2 if variant in ("C", "D") else 1
    p["h_sep"] = variant == "D"
    if variant != "DFT":
        fit = json.loads((HERE / f"fit_result_{variant}.json").read_text())["fit"]
        for k, v in fit.items():
            p[k] = (v, p[k][1]) if isinstance(p[k], tuple) else v
    return CoMKM(p, N=80)


def drc(m, d=0.01):
    """X_i = ∂ln(y)/∂ln(k_i), 以同时降低正逆反应活化能 d (eV) 实现 (平衡常数不变)。"""
    base = m.solve(T_REF, P_REF, R_REF)
    y0 = dict(TOF=math.log(base.tof_co), S_CH4=math.log(base.s_ch4), alpha=math.log(base.alpha()))
    lnf = d / (8.617333e-5 * T_REF)
    out = {}
    for s in STEPS:
        Ea, dE = m.p[s]
        vals = {}
        for sign in (+1, -1):
            m.p[s] = (Ea - sign * d, dE)
            r = m.solve(T_REF, P_REF, R_REF)
            vals[sign] = dict(TOF=math.log(r.tof_co), S_CH4=math.log(r.s_ch4), alpha=math.log(r.alpha()))
        m.p[s] = (Ea, dE)
        out[s] = {k: (vals[1][k] - vals[-1][k]) / (2 * lnf) for k in y0}
    m.solve(T_REF, P_REF, R_REF)
    return out


if __name__ == "__main__":
    res = {}
    for v in ("DFT", "A", "B", "C", "D"):
        m = model(v)
        r = m.solve(T_REF, P_REF, R_REF)
        n = np.arange(len(r.paraffin))
        tot = r.paraffin + r.olefin
        a_n = tot[2:] / tot[1:-1] if False else None
        # 链长 n 的增长概率 α_n = Σ_{k>n} / Σ_{k≥n} (按摩尔)
        cum = np.cumsum(tot[::-1])[::-1]
        alpha_n = cum[2:] / cum[1:-1]
        res[v] = dict(W=r.weight_frac().tolist(), mol=(tot / tot[1:].sum()).tolist(), alpha_n=alpha_n.tolist(),
                      op=[(r.olefin[k] / r.paraffin[k]) if r.paraffin[k] > 0 else None for k in range(len(n))],
                      tof=r.tof_co, s_ch4=r.s_ch4, alpha=r.alpha(), theta=r.theta)
        print(v, "TOF", f"{r.tof_co:.3g}", "S_CH4", f"{r.s_ch4:.3f}", "alpha10-20", f"{r.alpha():.3f}",
              "alpha_n(1,2,3,5,10,20)", [round(alpha_n[k - 1], 3) for k in (1, 2, 3, 5, 10, 20)],
              "O/P(2,3,5,10)", [round(r.olefin[k] / r.paraffin[k], 2) for k in (2, 3, 5, 10)])
        if v == BEST:
            dr = drc(m)
            res[v]["drc"] = dr
            for s, x in dr.items():
                if max(abs(t) for t in x.values()) > 0.03:
                    print(f"  DRC {s:5s} TOF {x['TOF']:+.2f}  S_CH4 {x['S_CH4']:+.2f}  alpha {x['alpha']:+.3f}")
    (HERE / "analysis.json").write_text(json.dumps(res, ensure_ascii=False))
