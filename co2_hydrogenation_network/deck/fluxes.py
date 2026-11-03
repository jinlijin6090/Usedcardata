"""Integrated reaction fluxes through the calibrated G1 lumped network (for the deck)."""
import json
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

import os
TK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tandem_kinetics")
sys.path.insert(0, TK)
import calibrate as C  # noqa: E402
from model import Model, k_dme, k_meoh, k_rwgs  # noqa: E402

cal = json.load(open(os.path.join(TK, "calibration.json")))
OUT = {}
for case in ("G1_GaZrOx_GaZSM5", "G2_ZnZrO_ZSM5", "B1_In2O3_HZSM5"):
    p = {k: tuple(v) for k, v in cal[case]["params"].items()}
    m = Model(p)
    net = m.net
    c = C.CASES[case]["cond"]
    T, Ptot = c.T_C + 273.15, c.P_MPa * 10
    F0 = m.feed(c)
    sol = solve_ivp(lambda w, f: m.rhs(w, f, c, T, Ptot), (0, c.W), F0, method="BDF",
                    rtol=1e-6, atol=1e-14 * F0.sum(), dense_output=True)
    I = net.idx
    W = np.linspace(0, c.W, 801)
    acc = np.zeros(len(net.rxns))
    ox = np.zeros(3)
    ht_partner = np.zeros(len(net.O_idx))
    for w0, w1 in zip(W[:-1], W[1:]):
        F = sol.sol(0.5 * (w0 + w1))
        Ft = F.sum()
        P = np.maximum(F, 0) / Ft * Ptot
        pco2, ph2, pco, pw, pm, pd = (P[I[x]] for x in ("CO2", "H2", "CO", "H2O", "MeOH", "DME"))
        inh = 1 + m.kT("K_w", T) * pw
        q1 = pm * pw / (pco2 * ph2 ** 3) / k_meoh(T)
        q2 = pco * pw / (pco2 * ph2) / k_rwgs(T)
        r1 = m.kT("k_meoh", T) * pco2 * ph2 * (1 - q1) / inh * c.phi_ox
        r2 = m.kT("k_rwgs", T) * pco2 * math.sqrt(ph2) * (1 - q2) / inh * c.phi_ox
        rd = m.kT("k_dme", T) * (pm * pm - pd * pw / k_dme(T)) * (1 - c.phi_ox)
        kv = np.array([m.kT(k, T) for k in net.knames])[net.kidx]
        pi = P[net.r_i]
        pj = np.where(net.r_j >= 0, P[np.maximum(net.r_j, 0)], 1.0)
        rate = kv * net.fac * pi * pj
        po = P[net.O_idx]
        rate[net.is_ht] *= po.sum()
        rate *= np.where(net.site_ox, c.phi_ox, 1 - c.phi_ox)
        dw = w1 - w0
        acc += rate * dw
        ox += np.array([r1, r2, rd]) * dw
        if po.sum() > 0:
            ht_partner += 3 * rate[net.is_ht].sum() * po / po.sum() * dw
    dco2 = F0[I["CO2"]] - sol.y[I["CO2"], -1]
    rx = []
    for r, (fam, reac, st, fac, kn, site) in enumerate(net.rxns):
        rx.append({"family": fam,
                   "reactants": {net.names[i]: o for i, o in reac.items()},
                   "stoich": {net.names[i]: v for i, v in st.items()},
                   "flux_rel": float(acc[r] / dco2)})   # mol per mol CO2 converted
    fam_tot = {}
    for x in rx:
        fam_tot[x["family"]] = fam_tot.get(x["family"], 0) + x["flux_rel"]
    OUT[case] = {"oxide": {"MeOH_synthesis": float(ox[0] / dco2), "RWGS": float(ox[1] / dco2),
                           "DME_net": float(ox[2] / dco2)},
                 "family_flux": fam_tot, "reactions": rx,
                 "ht_partner_flux": {f"ole{n}": float(v / dco2) for n, v in zip(range(2, net.N + 1), ht_partner)},
                 "outlet": {n: float(sol.y[i, -1] / dco2) for n, i in I.items()}}
    print(case, {k: round(v, 3) for k, v in OUT[case]["oxide"].items()},
          {k: round(v, 3) for k, v in sorted(fam_tot.items(), key=lambda kv: -kv[1])})
json.dump(OUT, open("fluxes.json", "w"), indent=1)
