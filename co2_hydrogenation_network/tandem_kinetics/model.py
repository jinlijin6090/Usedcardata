"""Lumped kinetic model of one-step CO2 hydrogenation over oxide/zeolite (OX-ZEO)
bifunctional catalysts, aimed at the jet-fuel (C8-C16) cut.

Mechanism (methanol-mediated tandem route; see README for the evidence):

  Oxide (ZnZrOx / GaZrOx / In2O3 ...), rates per g of oxide
    OX1  CO2 + 3 H2 <=> CH3OH + H2O          formate/methoxy route, reversible (thermo)
    OX2  CO2 + H2   <=> CO + H2O             RWGS, reversible (thermo)
    OX3  CnH2n + H2  -> CnH2n+2              olefin hydrogenation (H spill-over)
  Zeolite (H-ZSM-5 type acid sites), rates per g of zeolite
    Z1   2 CH3OH <=> CH3OCH3 + H2O           fast, reversible
    Z2   2 CH3OH  -> C2H4 + 2 H2O            first C-C bond (slow initiation)
    Z3   CnH2n + CH3OH/DME -> Cn+1H2n+2      olefin methylation (olefin cycle)
    Z4   CiH2i + CjH2j -> Ci+jH2(i+j)        oligomerization
    Z5   CnH2n -> CiH2i + Cn-iH2(n-i)        beta-scission cracking
    Z6   CnH2n -> cyclic CnH2n (naphthene)   cyclization (n >= 6)
    Z7   CnH2n(cyc) + 3 olefin -> ArCn + 3 paraffin   hydride transfer / aromatization
    Z7b  cyclic CnH2n-x -> ArCn + H2          dehydrogenative aromatization on Lewis (Ga/Zn)
                                             sites; the cyclic lump stands for cyclo-olefinic
                                             intermediates (cyclohexene -> benzene + 2 H2 has
                                             K(593 K) ~ 5e3 bar^2), written here as +3 H2
    Z8   ArCn + CH3OH/DME -> ArCn+1          aromatic methylation (aromatic cycle)
    Z9   ArCn -> ArCn-k + CkH2k (k = 2,3)    dealkylation / paring (aromatic cycle)
    Z10  CH3OH + H2 -> CH4 + H2O             methane by-product
  Shape selectivity: formation of a product with n carbons in the pores is
  multiplied by s(n) = 1 / (1 + exp((n - n_cap)/w)); separate caps for chains and
  aromatics (MFI: chains ~C11, aromatics ~C10).

Species are lumped by class and carbon number: olefins O2..ON, paraffins P2..PN,
naphthenes N6..NN, aromatics A6..NA. Isomers inside a lump are assumed to share
the lump's reactivity (single-class lumping); the iso/normal split of paraffins
is estimated afterwards from an equilibrium branching fraction.

Units: partial pressure bar, molar flow mol/h, catalyst mass g, rates mol/(g h).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

R = 8.314462618e-3   # kJ/(mol K)
TREF = 593.15        # K, 320 C

# --------------------------------------------------------------------------- thermo
def k_meoh(T):   # CO2 + 3H2 <=> MeOH + H2O, bar^-2  (Bussche & Froment 1996)
    return 10 ** (3066.0 / T - 10.592)


def k_rwgs(T):   # CO2 + H2 <=> CO + H2O  (inverse of Bussche-Froment WGS)
    return 10 ** (-2073.0 / T + 2.029)


def k_dme(T):    # 2 MeOH <=> DME + H2O (fit to Diep & Wainwright; ~20 at 593 K)
    return math.exp(2835.0 / T - 1.675)


# ----------------------------------------------------------------------- parameters
DEFAULT = {
    # name: (value at TREF, Ea kJ/mol)
    "k_meoh": (2.0e-4, 103.0),   # OX1, mol/(g h bar^2)      Ea ~103 (apparent, lit.)
    "k_rwgs": (6.0e-5, 117.0),   # OX2, mol/(g h bar^1.5)    Ea ~117 (apparent, lit.)
    "K_w": (0.05, 0.0),          # water inhibition on oxide, 1/bar
    "k_hyd": (2.0e-3, 60.0),     # OX3 olefin hydrogenation, mol/(g h bar^2)
    "k_dme": (5.0, 60.0),        # Z1, mol/(g h bar^2), fast
    "k_init": (1.0e-5, 140.0),   # Z2, mol/(g h bar^2)
    "k_meth": (0.30, 80.0),      # Z3 (n>=4), mol/(g h bar^2)
    "k_olig": (0.02, 70.0),      # Z4, mol/(g h bar^2)
    "k_crack": (0.05, 120.0),    # Z5 for C6, mol/(g h bar)
    "k_cyc": (0.01, 110.0),      # Z6, mol/(g h bar)
    "k_ht": (0.10, 90.0),        # Z7, mol/(g h bar^2)
    "k_dh": (0.05, 120.0),       # Z7b dehydrogenative aromatization, mol/(g h bar)
    "k_ameth": (0.05, 80.0),     # Z8, mol/(g h bar^2)
    "k_dealk": (0.01, 130.0),    # Z9 (from C10 aromatics), mol/(g h bar)
    "k_ch4": (2.0e-4, 120.0),    # Z10, mol/(g h bar^2)
}
SHAPE = {"n_cap": 11.0, "w": 0.7, "n_capA": 10.0, "wA": 0.5,
         "alpha_crack": 0.9,      # cracking rate x e^(alpha (n-6))
         "eps_eth": 0.05, "eps_prop": 0.5,   # relative methylation of C2=, C3=
         "alpha_dealk": 0.8}


@dataclass
class Conditions:
    T_C: float = 320.0
    P_MPa: float = 3.0
    H2_CO2: float = 3.0
    GHSV: float = 4000.0          # mL(STP)/(g_cat h), on total catalyst
    phi_ox: float = 0.5           # oxide mass fraction
    inert: float = 0.0            # mol inert per mol CO2 (e.g. Ar/N2 internal std)
    W: float = 1.0                # g catalyst (results scale-free in GHSV)


@dataclass
class Network:
    N: int = 20
    NA: int = 16
    names: list = field(default_factory=list)
    idx: dict = field(default_factory=dict)
    nC: np.ndarray = None
    cls: list = field(default_factory=list)
    rxns: list = field(default_factory=list)   # (family, reactants{i:order}, stoich{i:nu}, factor, kname, site)

    def __post_init__(self):
        base = [("CO2", 1, "CO2"), ("H2", 0, "H2"), ("CO", 1, "CO"), ("H2O", 0, "H2O"),
                ("MeOH", 1, "oxy"), ("DME", 2, "oxy"), ("CH4", 1, "P"), ("inert", 0, "inert")]
        for n, c, k in base:
            self._add(n, c, k)
        for n in range(2, self.N + 1):
            self._add(f"ole{n}", n, "O")
        for n in range(2, self.N + 1):
            self._add(f"par{n}", n, "P")
        for n in range(6, self.N + 1):
            self._add(f"nap{n}", n, "N")
        for n in range(6, self.NA + 1):
            self._add(f"aro{n}", n, "A")
        self.nC = np.array(self.nC_list)

    def _add(self, name, c, k):
        if not hasattr(self, "nC_list"):
            self.nC_list = []
        self.idx[name] = len(self.names)
        self.names.append(name)
        self.nC_list.append(c)
        self.cls.append(k)


def s_shape(n, cap, w):
    return 1.0 / (1.0 + math.exp((n - cap) / w))


def build_network(N=20, NA=16, shape=None):
    """Generate the lumped reaction list from the family rules above."""
    sh = dict(SHAPE, **(shape or {}))
    net = Network(N=N, NA=NA)
    I = net.idx
    add = net.rxns.append
    Ox = lambda n: I[f"ole{n}"]
    Pa = lambda n: I["CH4"] if n == 1 else I[f"par{n}"]
    # OX3 hydrogenation
    for n in range(2, N + 1):
        add(("OX3_hyd", {Ox(n): 1, I["H2"]: 1}, {Ox(n): -1, I["H2"]: -1, Pa(n): 1}, 1.0, "k_hyd", "ox"))
    # Z2 initiation
    add(("Z2_init", {I["MeOH"]: 2}, {I["MeOH"]: -2, Ox(2): 1, I["H2O"]: 2}, 1.0, "k_init", "zeo"))
    # Z10 methane
    add(("Z10_ch4", {I["MeOH"]: 1, I["H2"]: 1}, {I["MeOH"]: -1, I["H2"]: -1, I["CH4"]: 1, I["H2O"]: 1},
         1.0, "k_ch4", "zeo"))
    # Z3 methylation by MeOH and DME
    for n in range(2, N):
        f = sh["eps_eth"] if n == 2 else sh["eps_prop"] if n == 3 else 1.0
        f *= s_shape(n + 1, sh["n_cap"], sh["w"])
        add(("Z3_meth", {Ox(n): 1, I["MeOH"]: 1}, {Ox(n): -1, I["MeOH"]: -1, Ox(n + 1): 1, I["H2O"]: 1},
             f, "k_meth", "zeo"))
        add(("Z3_meth", {Ox(n): 1, I["DME"]: 1}, {Ox(n): -1, I["DME"]: -1, Ox(n + 1): 1, I["MeOH"]: 1},
             f, "k_meth", "zeo"))
    # Z4 oligomerization (i <= j, ethene sluggish)
    for i in range(2, N):
        for j in range(i, N + 1 - i):
            f = s_shape(i + j, sh["n_cap"], sh["w"])
            f *= (sh["eps_eth"] if i == 2 else 1.0) * (sh["eps_eth"] if j == 2 else 1.0)
            st = {Ox(i): -1, Ox(j): -1, Ox(i + j): 1} if i != j else {Ox(i): -2, Ox(i + j): 1}
            add(("Z4_olig", {Ox(i): 1, Ox(j): 1} if i != j else {Ox(i): 2}, st, f, "k_olig", "zeo"))
    # Z5 cracking: CnH2n -> CiH2i + C(n-i)H2(n-i), fragments >= C2, secondary/tertiary preferred
    for n in range(5, N + 1):
        splits = [(i, n - i) for i in range(2, n // 2 + 1) if n - i >= 2]
        wts = [1.0 if (a >= 3 and b >= 3) else 0.1 for a, b in splits]
        tot = sum(wts)
        base = math.exp(sh["alpha_crack"] * (n - 6))
        for (a, b), wt in zip(splits, wts):
            st = {Ox(n): -1}
            st[Ox(a)] = st.get(Ox(a), 0) + 1
            st[Ox(b)] = st.get(Ox(b), 0) + 1
            add(("Z5_crack", {Ox(n): 1}, st, base * wt / tot, "k_crack", "zeo"))
    # Z6 cyclization
    for n in range(6, N + 1):
        add(("Z6_cyc", {Ox(n): 1}, {Ox(n): -1, I[f"nap{n}"]: 1},
             s_shape(n, sh["n_cap"], sh["w"]), "k_cyc", "zeo"))
    # Z7 hydride transfer (special stoichiometry, handled in rates()) - one entry per naphthene
    for n in range(6, N + 1):
        tgt = I[f"aro{n}"] if n <= NA else None
        add(("Z7_HT", {I[f"nap{n}"]: 1}, {I[f"nap{n}"]: -1, **({tgt: 1} if tgt is not None else {})},
             1.0 if tgt is not None else 0.0, "k_ht", "zeo"))
    # Z7b dehydrogenative aromatization (Lewis-acid Ga/Zn sites)
    for n in range(6, min(N, NA) + 1):
        add(("Z7b_dh", {I[f"nap{n}"]: 1}, {I[f"nap{n}"]: -1, I[f"aro{n}"]: 1, I["H2"]: 3},
             1.0, "k_dh", "zeo"))
    # Z8 aromatic methylation
    for n in range(6, NA):
        f = s_shape(n + 1, sh["n_capA"], sh["wA"])
        add(("Z8_ameth", {I[f"aro{n}"]: 1, I["MeOH"]: 1},
             {I[f"aro{n}"]: -1, I["MeOH"]: -1, I[f"aro{n+1}"]: 1, I["H2O"]: 1}, f, "k_ameth", "zeo"))
        add(("Z8_ameth", {I[f"aro{n}"]: 1, I["DME"]: 1},
             {I[f"aro{n}"]: -1, I["DME"]: -1, I[f"aro{n+1}"]: 1, I["MeOH"]: 1}, f, "k_ameth", "zeo"))
    # Z9 dealkylation (paring), polymethylbenzenes release C2=/C3=
    for n in range(8, NA + 1):
        base = math.exp(sh["alpha_dealk"] * (n - 10))
        for k, wk in ((2, 0.4), (3, 0.6)):
            if n - k >= 6:
                add(("Z9_dealk", {I[f"aro{n}"]: 1}, {I[f"aro{n}"]: -1, I[f"aro{n-k}"]: 1, Ox(k): 1},
                     base * wk, "k_dealk", "zeo"))
    _vectorize(net)
    return net


def _vectorize(net):
    nr, ns = len(net.rxns), len(net.names)
    net.S = np.zeros((ns, nr))
    net.r_i = np.full(nr, -1)
    net.r_j = np.full(nr, -1)
    net.fac = np.zeros(nr)
    net.kname = []
    net.site_ox = np.zeros(nr, bool)
    net.is_ht = np.zeros(nr, bool)
    for r, (fam, reac, st, f, kn, site) in enumerate(net.rxns):
        for i, nu in st.items():
            net.S[i, r] += nu
        items = []
        for i, o in reac.items():
            items += [i] * o
        net.r_i[r] = items[0]
        net.r_j[r] = items[1] if len(items) > 1 else -1
        net.fac[r] = f
        net.kname.append(kn)
        net.site_ox[r] = site == "ox"
        net.is_ht[r] = fam == "Z7_HT"
    net.families = [x[0] for x in net.rxns]
    net.knames = sorted(set(net.kname))
    net.kidx = np.array([net.knames.index(k) for k in net.kname])
    O = [net.idx[f"ole{n}"] for n in range(2, net.N + 1)]
    P = [net.idx["CH4" if n == 1 else f"par{n}"] for n in range(2, net.N + 1)]
    net.O_idx, net.P_idx = np.array(O), np.array(P)


# --------------------------------------------------------------------------- model
class Model:
    def __init__(self, params=None, shape=None, N=20, NA=16):
        self.p = {k: tuple(v) for k, v in DEFAULT.items()}
        if params:
            for k, v in params.items():
                self.p[k] = v if isinstance(v, tuple) else (v, self.p[k][1])
        self.net = build_network(N, NA, shape)

    def kT(self, name, T):
        k0, ea = self.p[name]
        return k0 * math.exp(-ea / R * (1.0 / T - 1.0 / TREF))

    def rhs(self, W, F, c, T, Ptot):
        net = self.net
        I = net.idx
        Ft = F.sum()
        P = np.maximum(F, 0.0) / Ft * Ptot
        # ---------------- oxide reversible steps
        pco2, ph2, pco, pw, pm, pd = (P[I[x]] for x in ("CO2", "H2", "CO", "H2O", "MeOH", "DME"))
        inh = 1.0 + self.kT("K_w", T) * pw
        ph2s = max(ph2, 1e-12)
        q1 = pm * pw / (max(pco2, 1e-12) * ph2s ** 3) / k_meoh(T)
        q2 = pco * pw / (max(pco2, 1e-12) * ph2s) / k_rwgs(T)
        r1 = self.kT("k_meoh", T) * pco2 * ph2 * (1.0 - q1) / inh
        r2 = self.kT("k_rwgs", T) * pco2 * math.sqrt(ph2s) * (1.0 - q2) / inh
        rd = self.kT("k_dme", T) * (pm * pm - pd * pw / k_dme(T))
        dF = np.zeros_like(F)
        phi = c.phi_ox
        for nm, nu in (("CO2", -1), ("H2", -3), ("MeOH", 1), ("H2O", 1)):
            dF[I[nm]] += phi * nu * r1
        for nm, nu in (("CO2", -1), ("H2", -1), ("CO", 1), ("H2O", 1)):
            dF[I[nm]] += phi * nu * r2
        for nm, nu in (("MeOH", -2), ("DME", 1), ("H2O", 1)):
            dF[I[nm]] += (1 - phi) * nu * rd
        # ---------------- lumped irreversible network
        kv = np.array([self.kT(k, T) for k in net.knames])[net.kidx]
        pi = P[net.r_i]
        pj = np.where(net.r_j >= 0, P[np.maximum(net.r_j, 0)], 1.0)
        rate = kv * net.fac * pi * pj
        # hydride transfer: r = k P_N * sum(P_olefin); consumes 3 olefins -> 3 paraffins
        po = P[net.O_idx]
        spo = po.sum()
        rate[net.is_ht] *= spo
        wsite = np.where(net.site_ox, phi, 1 - phi)
        rate = rate * wsite
        dF += net.S @ rate
        r_ht = rate[net.is_ht].sum()
        if spo > 0:
            frac = po / spo
            dF[net.O_idx] -= 3 * r_ht * frac
            dF[net.P_idx] += 3 * r_ht * frac
        self._last = (r1 * phi, r2 * phi)
        return dF

    def feed(self, c):
        I = self.net.idx
        F = np.zeros(len(self.net.names))
        tot = c.GHSV * c.W / 22414.0          # mol/h of feed gas
        x_co2 = 1.0 / (1.0 + c.H2_CO2 + c.inert)
        F[I["CO2"]] = tot * x_co2
        F[I["H2"]] = tot * x_co2 * c.H2_CO2
        F[I["inert"]] = tot * x_co2 * c.inert + 1e-12
        F[I["ole2"]] = 1e-12 * tot       # trace olefin seed (autocatalysis)
        return F

    def run(self, c: Conditions, n_out=1):
        T, Ptot = c.T_C + 273.15, c.P_MPa * 10.0
        F0 = self.feed(c)
        sol = solve_ivp(lambda w, f: self.rhs(w, f, c, T, Ptot), (0, c.W), F0, method="BDF",
                        rtol=1e-6, atol=1e-14 * F0.sum(), dense_output=False)
        if not sol.success:
            raise RuntimeError(sol.message)
        return Result(self.net, F0, sol.y[:, -1], c)


class Result:
    def __init__(self, net, F0, F, c):
        self.net, self.F0, self.F, self.c = net, F0, np.maximum(F, 0), c
        I = net.idx
        self.conv = 1 - self.F[I["CO2"]] / F0[I["CO2"]]
        dco2 = F0[I["CO2"]] - self.F[I["CO2"]]
        self.S_CO = self.F[I["CO"]] / dco2
        self.S_oxy = (self.F[I["MeOH"]] + 2 * self.F[I["DME"]]) / dco2
        hc = [i for i, k in enumerate(net.cls) if k in ("O", "P", "N", "A")]
        self.hc = hc
        cflow = self.F * net.nC
        self.C_hc = cflow[hc].sum()
        self.S_HC = self.C_hc / dco2

    def dist(self):
        """Carbon-based selectivity within hydrocarbons, by class and carbon number."""
        net = self.net
        out = {}
        for i in self.hc:
            n, k = net.nC[i], net.cls[i]
            out.setdefault(k, {})
            out[k][n] = self.F[i] * n / self.C_hc
        return out

    def summary(self, iso_frac=None):
        d = self.dist()
        def tot(kset, lo=1, hi=99):
            return sum(v for k in kset for n, v in d.get(k, {}).items() if lo <= n <= hi)
        s = {
            "X_CO2": self.conv, "S_CO": self.S_CO, "S_MeOH+DME": self.S_oxy, "S_HC": self.S_HC,
            "CH4": tot("P", 1, 1), "C2-C4=": tot("O", 2, 4), "C2-C4o": tot("P", 2, 4),
            "C5+": tot("OPNA", 5), "aromatics": tot("A"), "olefins_C5+": tot("O", 5),
            "paraffins_C5+": tot("P", 5), "naphthenes": tot("N"),
            "C8-C16": tot("OPNA", 8, 16), "C5-C11": tot("OPNA", 5, 11), "C17+": tot("OPNA", 17),
        }
        cut = {k: tot(k, 8, 16) for k in "OPNA"}
        tc = sum(cut.values()) or 1.0
        s["jet_cut_comp"] = {k: v / tc for k, v in cut.items()}
        # SAF-usable fraction: C8-C16 with aromatics capped at 25 % of the cut (ASTM D7566 blend limit)
        nonar = cut["O"] + cut["P"] + cut["N"]
        s["SAF_usable"] = nonar + min(cut["A"], nonar * 0.25 / 0.75)
        s["Y_C8-C16"] = s["C8-C16"] * s["S_HC"] * s["X_CO2"]
        ar = d.get("A", {})
        car = sum(ar.values())
        s["nA_mean"] = (sum(ar.values()) / sum(v / n for n, v in ar.items())) if car > 0 else 0.0
        return s
