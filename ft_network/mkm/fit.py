#!/usr/bin/env python3
"""在 DFT 先验附近拟合 MKM 参数, 使其同时复现实验 TOF、CH4 选择性、α、表观活化能、反应级数和 C3 烯烷比。

目标函数 = Σ[(预测 − 实验)/容差]² + Σ[(参数 − DFT 值)/σ]²  (贝叶斯式正则化, 防止参数跑到不物理的位置)
只调整活化能 (ΔE 固定) 和吸附参数。
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from co_fts_mkm import DEFAULT, CoMKM, kinetics
from params_lit import DFT_STEP, P_REF, R_REF, T_REF, TARGETS, dft_params

HERE = Path(__file__).resolve().parent
NAMES = [k for k, *_ in DFT_STEP]
X0 = np.array([v[0] if isinstance(v, tuple) else v for _, v, _, _ in DFT_STEP])
SIG = np.array([s for _, _, s, _ in DFT_STEP])
LO = np.where(np.char.startswith(np.array(NAMES), "w_") | np.isin(NAMES, ["c_des", "d_C2"]), 0.0, 0.1)
HI = X0 + 4 * SIG


def to_params(x):
    p = dict(DEFAULT, **dft_params())
    p["nsite_dir"] = 2 if os.environ.get("FIT_VARIANT", "A") == "C" else 1
    for (k, v, _, _), xi in zip(DFT_STEP, x):
        p[k] = (float(xi), v[1]) if isinstance(v, tuple) else float(xi)
    return p


def predict(m, T=T_REF, P=P_REF, R=R_REF):
    r = m.solve(T, P, R)
    Ea, nH, nCO = kinetics(m, T, P, R)
    return dict(log10_TOF=math.log10(r.tof_co), S_CH4=r.s_ch4, alpha=r.alpha(), Ea=Ea, n_H2=nH, n_CO=nCO,
                log10_OP3=math.log10(r.op_ratio(3)), theta_CO=r.theta["CO"],
                S_C5=float(r.carbon_sel()[5:].sum()), tail=float(r.carbon_sel()[-1])), r


# 变体 B: 额外要求 CO 覆盖度约 0.4–0.5 ML (反应条件下 Co 表面 CO 覆盖度高; 软约束)
VARIANT = os.environ.get("FIT_VARIANT", "A")
if VARIANT in ("B", "C"):
    TARGETS = TARGETS + [("theta_CO", 0.45, 0.08, "软约束: 反应条件下 CO 覆盖度约 0.4–0.5 ML")]


NCHAIN = 80                                   # 链长上限; α = 0.9 时 C80 以上的碳 < 0.5 %
MODEL = CoMKM(N=NCHAIN)


def residuals(x):
    MODEL.p = to_params(x)
    try:
        pred, r = predict(MODEL)
        if r.resid > 1e-5 or not all(np.isfinite(list(pred.values()))):
            raise RuntimeError
    except Exception:
        MODEL._last = None
        return np.full(len(TARGETS) + 1 + len(x), 30.0)
    res = [(pred[k] - t) / tol for k, t, tol, _ in TARGETS]
    res.append(max(0.0, pred["tail"] - 0.01) / 0.003)          # 防护: 截断处 (C_N) 的碳不得超过 1 %
    return np.concatenate([res, (x - X0) / SIG])


def main():
    MODEL.p = to_params(X0)
    pred0, _ = predict(MODEL)
    start = X0.copy()
    if VARIANT == "C":                          # 变体 C: CO 解离需要两个空位
        MODEL.p["nsite_dir"] = 2
    if VARIANT in ("B", "C"):                   # 从 CO 吸附更强的起点出发
        start[NAMES.index("E_CO")] = 1.40
        start[NAMES.index("w_CO")] = 1.00
        MODEL._last = None
    fit = least_squares(residuals, start, bounds=(LO, HI), x_scale=SIG, diff_step=1e-3, max_nfev=400, verbose=1)
    MODEL._last = None
    MODEL.p = to_params(fit.x)
    pred1, r1 = predict(MODEL)
    out = dict(
        dft=dict(zip(NAMES, X0.tolist())), fit=dict(zip(NAMES, fit.x.tolist())), sigma=dict(zip(NAMES, SIG.tolist())),
        pred_dft=pred0, pred_fit=pred1, targets={k: [t, tol] for k, t, tol, _ in TARGETS},
        cost=float(fit.cost), theta_fit=r1.theta,
    )
    (HERE / f"fit_result_{VARIANT}.json").write_text(json.dumps(out, ensure_ascii=False, indent=1))
    print(f"{'指标':12s}{'实验':>10s}{'DFT':>10s}{'拟合':>10s}")
    for k, t, tol, _ in TARGETS:
        print(f"{k:12s}{t:10.3f}{pred0[k]:10.3f}{pred1[k]:10.3f}")
    for k in ("theta_CO", "tail"):
        print(f"{k:12s}{'':10s}{pred0[k]:10.3f}{pred1[k]:10.3f}")
    print("\n参数变化 (eV, 以 σ 计)")
    for k, a, b, s in zip(NAMES, X0, fit.x, SIG):
        print(f"{k:8s} {a:6.3f} → {b:6.3f}   Δ={b - a:+.3f}  ({(b - a) / s:+.1f}σ)")
    print("覆盖度", {k: f"{v:.3g}" for k, v in r1.theta.items()})
    return 0


if __name__ == "__main__":
    sys.exit(main())
