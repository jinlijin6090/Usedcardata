"""
Scenario C: growth-type carbon (whisker / filament and surface carbon islands).

Abstraction of the classical Ni whisker mechanism plus surface encapsulation:

    surface:   C* + C*  -> C2* (R18)  -> 2 Cg*  (R20b)   island nucleation (slow)
               C* + Cg* -> 2 Cg*      (R24)               island edge growth (autocatalytic)
               Cg* + O* -> CO* + *    (R21)               island gasification
    bulk:      C* + (bulk site) <-> C_b + *   (R22)       dissolution into Ni (barrier Ea_sol)
               C_b  -> C_fil         (R23)                diffusion-limited precipitation on the
                                                          rear-side filament, k = D_C / r_m^2,
                                                          driving force (x - x_eq), only once a
                                                          filament nucleus exists
    nucleation: df/dt = (S-1)^2 / tau_nuc * (1-f),  S = x / x_eq (supersaturation)
                f = fraction of Ni particles carrying a growing filament

Cg* blocks Ni sites (encapsulation / carbon particles on the surface); C_fil does not
(whisker carbon grows away from the active surface). The ODE here is normalised by the
site counts of each site type (the published engine integrates un-normalised rates,
which is fine for steady states but not for time-resolved behaviour).
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp
from mkmc import network as net
from mkmc.coke import CokeParams, load_model, KBT_H
from mkmc.config import ProcessParams, R_EV, NA
from mkmc.sites import calculate_sites
from mkmc.rates import rate_vector, net_rates


@dataclass
class FilamentParams:
    Ea_sol: float = 1.20     # C* -> C_b barrier [eV]  (lit. 1.0-1.6 eV surface->subsurface on Ni)
    dH_sol: float = 0.10     # C* -> C_b enthalpy [eV]
    dG_gr: float = 0.50      # graphite -> C_b dissolution free energy [eV] -> x_eq = exp(-dG/kT) (~0.26% at 973 K)
    D0: float = 1.2e-5       # C diffusion in Ni, pre-exponential [m^2/s]  (Lander et al. 1952: 0.12 cm^2/s)
    E_D: float = 1.41        # C diffusion in Ni, activation energy [eV]
    E_nuc: float = 2.50      # filament nucleation barrier [eV]: J = A_nuc exp(-E_nuc/kT) (S-1)^2 per particle
    A_nuc: float = KBT_H
    Ea_edge: float = 1.00    # C* + Cg(edge)* -> 2 Cg* barrier [eV] (edge growth of surface islands)
    A_edge: float = KBT_H
    A_sol: float = KBT_H
    c_edge: float = 2 * np.sqrt(np.pi)   # edge atoms of a disc island of N atoms ~ c_edge*sqrt(N)
    whisker: bool = True     # switch the bulk route on/off
    islands: bool = True     # switch the surface-island route on/off

    def with_(self, **kw): return replace(self, **kw)


class FilamentModel:
    def __init__(self, pp: ProcessParams, cp: CokeParams | None = None, fp: FilamentParams | None = None):
        net.use_variant("blocking")
        self.pp = pp
        self.cp = cp or CokeParams(Ea_cc=2.2, Ea_gas=1.8, dH_g=-0.40)   # C* -> Cg* = dH_sol + (-dG_gr) = -0.40
        self.fp = fp or FilamentParams()
        self.rc = load_model(self.cp)
        T = pp.T; kT = R_EV * T
        s = calculate_sites(*pp.as_site_args())
        self.geom = net.build_geom(s, pp.Zsm, pp.Zms)
        self.nsm = s.NSm / NA; self.nss = s.NSs / NA            # mol sites / g
        # bulk interstitial sites: one per Ni atom
        self.nb = pp.m_load * pp.cat_wt / 58.693
        self.Zmm = s.Zmm
        self.kf22 = self.fp.A_sol * np.exp(-self.fp.Ea_sol / kT)
        self.kr22 = self.kf22 / np.exp(-self.fp.dH_sol / kT)
        self.x_eq = np.exp(-self.fp.dG_gr / kT)
        D = self.fp.D0 * np.exp(-self.fp.E_D / kT)
        self.k_diff = D / (pp.rm * 1e-9) ** 2
        self.k24 = self.fp.A_edge * np.exp(-self.fp.Ea_edge / kT)
        self.knuc = self.fp.A_nuc * np.exp(-self.fp.E_nuc / kT)
        # island gasification R21 is re-expressed on edge sites below: switch the lattice form off
        self.k21f, self.k21r = self.rc.k[40], self.rc.k[41]
        self.rc.k[40] = 0.0; self.rc.k[41] = 0.0
        S = net.SPECIES; self.i = {k: j for j, k in enumerate(S)}
        self.n_sp = len(S); self.n_metal = net.N_METAL
        self.norm = np.array([self.nsm] * self.n_metal + [self.nss] * (self.n_sp - self.n_metal))

    def edge(self, th, n_isl):
        """Edge-site coverage of the surface carbon islands (perimeter scaling)."""
        thCg = max(th[self.i["Cg*"]], 0.0); n = max(n_isl, 0.0)
        return min(self.fp.c_edge * np.sqrt(n * thCg), thCg)

    def extra_rates(self, th, x, f, n_isl):
        i = self.i; fp = self.fp
        r22 = (self.kf22 * self.nsm * th[i["C*"]] * (1 - x) - self.kr22 * self.nsm * th[i["*"]] * x) if fp.whisker else 0.0
        r23 = f * self.k_diff * self.nb * max(x - self.x_eq, 0.0) if fp.whisker else 0.0
        th_e = self.edge(th, n_isl)
        r24 = self.k24 * self.Zmm * self.nsm * th[i["C*"]] * th_e if fp.islands else 0.0
        r21 = self.Zmm * self.nsm * (self.k21f * th_e * th[i["O*"]] - self.k21r * th[i["CO*"]] * th[i["*"]] * th_e / max(th[i["Cg*"]], 1e-30)) if fp.islands else 0.0
        S = x / self.x_eq
        jnuc = self.knuc * max(S - 1, 0.0) ** 2 if fp.whisker else 0.0
        return r22, r23, r24, jnuc, r21

    def rhs(self, t, y):
        th = y[:self.n_sp]; x, f, cfil, n_isl = y[self.n_sp:]
        r = rate_vector(th, self.rc.k, self.pp.p, self.geom)
        rn = net_rates(r)
        if not self.fp.islands:
            rn[17:] = 0.0                      # no C2*/Cg* chemistry at all
        d = (net.STOICH @ rn) / self.norm
        r22, r23, r24, jnuc, r21 = self.extra_rates(th, x, f, n_isl)
        i = self.i
        d[i["C*"]] += (-r22 - r24) / self.nsm
        d[i["*"]] += (r22 + r21) / self.nsm
        d[i["Cg*"]] += (r24 - r21) / self.nsm
        d[i["O*"]] += -r21 / self.nsm
        d[i["CO*"]] += r21 / self.nsm
        dx = (r22 - r23) / self.nb
        df = jnuc * (1 - f)
        dCfil = r23
        dn = rn[19] / self.nsm                 # one island per C2* -> 2Cg* event
        return np.concatenate([d, [dx, df, dCfil, dn]])

    def y0(self):
        y = np.zeros(self.n_sp + 4); y[self.i["*"]] = 1; y[self.i["Ov"]] = 0.7; y[self.i["v"]] = 0.3
        return y

    def run(self, t_end=1e5, n=400):
        t_eval = np.logspace(-3, np.log10(t_end), n); t_eval[-1] = t_end
        res = solve_ivp(self.rhs, (0, t_end), self.y0(), method="BDF", t_eval=t_eval, rtol=1e-7, atol=1e-12)
        return self.summarise(res)

    def summarise(self, res):
        t = res.t; Y = res.y.T; i = self.i
        rows = []
        for y in Y:
            th = y[:self.n_sp]; x, f, cfil, n_isl = y[self.n_sp:]
            r = rate_vector(th, self.rc.k, self.pp.p, self.geom); rn = net_rates(r)
            if not self.fp.islands: rn[17:] = 0.0
            r22, r23, r24, jnuc, r21 = self.extra_rates(th, x, f, n_isl)
            rows.append(dict(r_DRM=rn[5], r_CH4=rn[0], r_CO=rn[9] + rn[12] + r21, th_C=th[i["C*"]], th_Cg=th[i["Cg*"]], th_free=th[i["*"]],
                             th_O=th[i["O*"]], th_OHv=th[i["OHv"]], x_bulk=x, S=x / self.x_eq, f_nuc=f, C_fil=cfil, n_isl=n_isl,
                             island_size=th[i["Cg*"]] / n_isl if n_isl > 0 else 0.0, th_edge=self.edge(th, n_isl),
                             r_dissolve=r22, r_fil=r23, r_island=r24 + 2 * rn[19], r_gasify=r21))
        import pandas as pd
        df = pd.DataFrame(rows); df.insert(0, "t", t)
        df["C_island_mgC_g"] = df.th_Cg * self.nsm * 12.011e3
        df["C_fil_mgC_g"] = df.C_fil * 12.011e3
        df["C_bulk_mgC_g"] = df.x_bulk * self.nb * 12.011e3
        df["fil_rate_mgC_g_h"] = df.r_fil * 12.011e3 * 3600
        return df
