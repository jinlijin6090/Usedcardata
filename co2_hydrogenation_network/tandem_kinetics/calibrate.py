"""Calibrate the lumped OX-ZEO model to literature data points (abstract-level data).

Targets (carbon selectivities of hydrocarbons exclude CO, as in the papers):
  G1  Guo group, GaZrOx/Ga-ZSM-5 (Adv. Sci. 2026): 320 C, 3 MPa, X=17.2 %, S_CO=11.7 %,
      C5+ = 83.2 % of HC, CH4 < 1 %.          [H2/CO2, GHSV, mass ratio not in abstract]
  G2  Guo/Song, ZnO/ZrO2-ZSM-5 (J. CO2 Util. 2019): 340 C, X = 9 %, aromatics = 70 % of HC,
      CH4 < 1 %.                               [P, H2/CO2, GHSV, CO selectivity not in abstract]
  B1  benchmark, not Guo group: In2O3/HZSM-5 (Gao et al., Nat. Chem. 2017, 10.1038/nchem.2794):
      340 C, 3 MPa, H2/CO2 = 3, GHSV 9000 mL/(g h), X = 13.1 %, S_CO = 44.8 %, C5+ = 78.6 %,
      CH4 = 1 %.
Missing conditions are filled with the assumptions in ASSUMED (typical OX-ZEO values) and
are part of the uncertainty. The fit is regularized toward the prior (DEFAULT) in log space,
because each data point has only 3-5 numbers for 12 rate constants.
"""
from __future__ import annotations

import json
import os

import numpy as np
from scipy.optimize import least_squares

from model import DEFAULT, Conditions, Model

HERE = os.path.dirname(os.path.abspath(__file__))
FIT = ["k_meoh", "k_rwgs", "k_hyd", "k_init", "k_meth", "k_olig", "k_crack", "k_cyc", "k_ht",
       "k_dh", "k_ameth", "k_dealk", "k_ch4"]

# zeolite prior: scaled-up from DEFAULT (light-off needs ~10-100x, see README)
PRIOR = dict(DEFAULT)
for k, f in {"k_init": 300, "k_meth": 30, "k_olig": 30, "k_crack": 3, "k_cyc": 10, "k_ht": 10,
             "k_ameth": 10, "k_dealk": 3, "k_ch4": 0.3}.items():
    PRIOR[k] = (DEFAULT[k][0] * f, DEFAULT[k][1])

ASSUMED = {"H2_CO2": 3.0, "GHSV": 4000.0, "phi_ox": 0.5, "P_MPa": 3.0}

CASES = {
    "G1_GaZrOx_GaZSM5": {
        "cond": Conditions(T_C=320, P_MPa=3.0, H2_CO2=3.0, GHSV=4000, phi_ox=0.5),
        "targets": {"X_CO2": (0.172, 0.01), "S_CO": (0.117, 0.01), "C5+": (0.832, 0.02),
                    "CH4": (0.008, 0.005), "S_MeOH+DME": (0.01, 0.01)},
        "source": "Guo group, Adv. Sci. 2026 (abstract)"},
    "G2_ZnZrO_ZSM5": {
        "cond": Conditions(T_C=340, P_MPa=3.0, H2_CO2=3.0, GHSV=4000, phi_ox=0.5),
        "targets": {"X_CO2": (0.09, 0.01), "aromatics": (0.70, 0.03), "CH4": (0.008, 0.005),
                    "S_MeOH+DME": (0.01, 0.01)},
        "source": "Guo/Song, J. CO2 Util. 2019 (abstract)"},
    "B1_In2O3_HZSM5": {
        "cond": Conditions(T_C=340, P_MPa=3.0, H2_CO2=3.0, GHSV=9000, phi_ox=0.5),
        "targets": {"X_CO2": (0.131, 0.01), "S_CO": (0.448, 0.02), "C5+": (0.786, 0.02),
                    "CH4": (0.01, 0.005), "S_MeOH+DME": (0.01, 0.01)},
        "source": "benchmark (not Guo): Gao et al., Nat. Chem. 2017"},
}
# Soft assumption, not from the three abstracts: OX-ZEO aromatics are mostly C8-C10
# (xylenes, trimethylbenzenes), so the molar-mean aromatic carbon number is set to 8.7 +- 0.5.
for _c in CASES.values():
    _c["targets"]["nA_mean"] = (8.7, 0.5)
LAMBDA = 0.25   # ridge weight per decade of deviation from prior


def params_from(x):
    return {k: (PRIOR[k][0] * 10 ** xi, PRIOR[k][1]) for k, xi in zip(FIT, x)}


def residuals(x, case):
    try:
        s = Model(params_from(x)).run(case["cond"]).summary()
    except Exception:
        return np.full(len(case["targets"]) + len(x), 10.0)
    res = [(s[k] - v) / sd for k, (v, sd) in case["targets"].items()]
    return np.array(res + list(LAMBDA * np.asarray(x)))


def fit(name):
    case = CASES[name]
    x0 = np.zeros(len(FIT))
    sol = least_squares(residuals, x0, args=(case,), bounds=(-3, 3), diff_step=0.05,
                        x_scale=1.0, max_nfev=80)
    p = params_from(sol.x)
    s = Model(p).run(case["cond"]).summary()
    return {"case": name, "source": case["source"], "log10_shift": dict(zip(FIT, sol.x.round(3))),
            "params": {k: [v[0], v[1]] for k, v in p.items()},
            "fit": {k: s[k] for k in case["targets"]},
            "targets": {k: v[0] for k, v in case["targets"].items()},
            "prediction": {k: (v if not isinstance(v, dict) else v) for k, v in s.items()},
            "cond": vars(case["cond"]), "cost": float(sol.cost)}


if __name__ == "__main__":
    from multiprocessing import Pool
    with Pool(len(CASES)) as pool:
        results = pool.map(fit, list(CASES))
    out = {}
    for name, r in zip(CASES, results):
        out[name] = r
        print(name, "cost", round(r["cost"], 3))
        for k in r["targets"]:
            print(f"   {k:12s} target {r['targets'][k]:.3f}  fit {r['fit'][k]:.3f}")
        pr = r["prediction"]
        print("   predicted:", {k: round(pr[k], 3) for k in
                                ("C2-C4=", "C2-C4o", "C5+", "aromatics", "C8-C16", "SAF_usable")})
    with open(os.path.join(HERE, "calibration.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
