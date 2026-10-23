"""
Reaction network defined as data, not code.

The original model hand-wrote 34 rate lines (``r[0] = k[0]*p[0]*theta[0]*...``)
and 15 ODE lines (``dt[0] = -(r[0]-r[1]) - ...``).  Both are reconstructed here
from a few matrices so that adding/altering a step means editing a table, and
the rate evaluation becomes a couple of vectorised NumPy operations.

Layout
------
* 15 surface species (``SPECIES``)
* 5 gas partial pressures (``PRESSURES``)
* 17 reversible elementary reactions -> 34 directional rates
  (forward = even index 2n, reverse = odd index 2n+1)

Matrices
--------
``COV_EXP`` : (34, 15) integer powers of each coverage in each directional rate
``P_EXP``   : (34, 5)  integer powers of each partial pressure
``COORD_KIND`` / ``SITE_KIND`` : the geometric prefactor of each directional
              rate, expressed symbolically so the numeric value is built once
              per geometry by :func:`build_geom`.
``STOICH``  : (15, 17) net stoichiometric matrix; dtheta/dt = STOICH @ net_rate
"""
from __future__ import annotations

import numpy as np

from .sites import SiteCounts

# --------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------
SPECIES = [
    "*", "CH4*", "CH3*", "CH2*", "CH*", "C*", "H*", "O*", "CO*", "CHO*", "C2*",
    "CO2v", "CO3v", "OHv", "Ov", "v",
]
N_METAL = 11            # first 11 species live on Ni sites, the rest on CeO2
PRESSURES = ["CH4", "CO2", "H2", "CO", "H2O"]
REACTIONS = [f"R{i + 1}" for i in range(20)]

N_SPECIES = len(SPECIES)        # 15
N_PRESSURE = len(PRESSURES)     # 5
N_REACTIONS = len(REACTIONS)    # 17
N_RATES = 2 * N_REACTIONS       # 34

_S = {name: i for i, name in enumerate(SPECIES)}     # species index lookup

# --------------------------------------------------------------------------
# Per-directional-rate definition.
# Each entry: (coverage powers as {species: power},
#              {pressure_index: power} or {},
#              coordination kind, site kind)
# The order is exactly r[0]..r[33] of the original `rates()` function.
#An empty dict {} means no gas products are released in that direction. Example: {0: 1} releases one H₂ molecule.
#"one" is unimolecular
# --------------------------------------------------------------------------
_RATE_DEFS = [
    # R1  CH4(g) + * -> CH4*
    ({"*": 1}, {0: 1}, "one", "NSm"),          # r0  f
    ({"CH4*": 1}, {}, "one", "NSm"),           # r1  r
    # R2  CH4* + * -> CH3* + H*
    ({"CH4*": 1, "*": 1}, {}, "Zmm", "NSm"),   # r2
    ({"CH3*": 1, "H*": 1}, {}, "Zmm", "NSm"),  # r3
    # R3  CH3* + * -> CH2* + H*
    ({"CH3*": 1, "*": 1}, {}, "Zmm", "NSm"),   # r4
    ({"CH2*": 1, "H*": 1}, {}, "Zmm", "NSm"),  # r5
    # R4  CH2* + * -> CH* + H*
    ({"CH2*": 1, "*": 1}, {}, "Zmm", "NSm"),   # r6
    ({"CH*": 1, "H*": 1}, {}, "Zmm", "NSm"),   # r7
    # R5  CH* + * -> C* + H*
    ({"CH*": 1, "*": 1}, {}, "Zmm", "NSm"),    # r8
    ({"C*": 1, "H*": 1}, {}, "Zmm", "NSm"),    # r9
    # R6  2H* -> H2(g) + 2*
    ({"H*": 2}, {}, "Zmm", "NSm"),             # r10
    ({"*": 2}, {2: 1}, "Zmm", "NSm"),          # r11
    # R7  CH* + O* -> CHO*
    ({"CH*": 1, "O*": 1}, {}, "Zmm", "NSm"),   # r12
    ({"CHO*": 1, "*": 1}, {}, "Zmm", "NSm"),   # r13
    # R8  CHO* + * -> CO* + H*
    ({"CHO*": 1, "*": 1}, {}, "Zmm", "NSm"),   # r14
    ({"CO*": 1, "H*": 1}, {}, "Zmm", "NSm"),   # r15
    # R9  C* + O* -> CO* + *
    ({"C*": 1, "O*": 1}, {}, "Zmm", "NSm"),    # r16
    ({"CO*": 1, "*": 1}, {}, "Zmm", "NSm"),    # r17
    # R10 CO* -> CO + *
    ({"CO*": 1}, {}, "one", "NSm"),            # r18
    ({"*": 1}, {3: 1}, "one", "NSm"),          # r19
    # R11 Ov + * -> O* + v
    ({"*": 1, "Ov": 1}, {}, "Zms", "NIs"),     # r20
    ({"v": 1, "O*": 1}, {}, "Zms", "NIs"),     # r21
    # R12 CO2(g) + v -> CO2v
    ({"v": 1}, {1: 1}, "one", "NSs"),          # r22
    ({"CO2v": 1}, {}, "one", "NSs"),           # r23
    # R13 CO2v -> CO + Ov
    ({"CO2v": 1}, {}, "one", "NSs"),           # r24
    ({"Ov": 1}, {3: 1}, "one", "NSs"),         # r25
    # R14 CO2(g) + Ov -> CO3v
    ({"Ov": 1}, {1: 1}, "one", "NSs"),         # r26
    ({"CO3v": 1}, {}, "one", "NSs"),           # r27
    # R15 CO3v + * -> CO2v + O*
    ({"*": 1, "CO3v": 1}, {}, "Zms", "NIs"),   # r28
    ({"CO2v": 1, "O*": 1}, {}, "Zms", "NIs"),  # r29
    # R16 H* + Ov -> OHv + *
    ({"Ov": 1, "H*": 1}, {}, "Zsm", "NIm"),    # r30
    ({"*": 1, "OHv": 1}, {}, "Zsm", "NIm"),    # r31
    # R17 OHv + H* -> H2O(g) + v + *
    ({"OHv": 1, "H*": 1}, {}, "Zsm", "NIm"),   # r32
    ({"*": 1, "v": 1}, {4: 1}, "Zsm", "NIm"),  # r33
    # ---------------- coking extension (metal sites) ----------------
    # R18 C* + C* -> C2* + *        (C-C coupling, coke nucleation)
    ({"C*": 2}, {}, "Zmm", "NSm"),             # r34 f
    ({"C2*": 1, "*": 1}, {}, "Zmm", "NSm"),    # r35 r
    # R19 C2* + O* -> CO* + C*      (oxidative scission of the dimer)
    ({"C2*": 1, "O*": 1}, {}, "Zmm", "NSm"),   # r36 f
    ({"CO*": 1, "C*": 1}, {}, "Zmm", "NSm"),   # r37 r
    # R20 C2* -> C(dep) + *         (irreversible growth into graphitic /
    #                                 filamentous coke; carbon leaves the lattice)
    ({"C2*": 1}, {}, "one", "NSm"),            # r38 f
    ({"*": 1}, {}, "one", "NSm"),              # r39 r  (k set to 0)
]
# Blocking variant (scenario B): deposited carbon stays on Ni sites as Cg*
_RATE_DEFS_BLOCKING = _RATE_DEFS[:38] + [
    # R20b C2* + * -> 2 Cg*          (graphitic carbon occupying metal sites)
    ({"C2*": 1, "*": 1}, {}, "Zmm", "NSm"),    # r38 f
    ({"Cg*": 2}, {}, "Zmm", "NSm"),            # r39 r (k = 0)
    # R21  Cg* + O* -> CO* + *       (slow gasification of graphitic carbon)
    ({"Cg*": 1, "O*": 1}, {}, "Zmm", "NSm"),   # r40 f
    ({"CO*": 1, "*": 1}, {}, "Zmm", "NSm"),    # r41 r
]
assert len(_RATE_DEFS) == N_RATES # Ensure the number of rate definitions matches the expected number

VARIANT = "flux"

def _build_exponents(defs, species):
    S = {name: i for i, name in enumerate(species)}
    cov_exp = np.zeros((len(defs), len(species)), dtype=int)
    p_exp = np.zeros((len(defs), N_PRESSURE), dtype=int)
    coord, site = [], []
    for j, (cov, pres, c, st) in enumerate(defs):
        for sp, power in cov.items():
            cov_exp[j, S[sp]] = power
        for pidx, power in pres.items():
            p_exp[j, pidx] = power
        coord.append(c); site.append(st)
    return cov_exp, p_exp, coord, site

COV_EXP, P_EXP, COORD_KIND, SITE_KIND = _build_exponents(_RATE_DEFS, SPECIES)

# --------------------------------------------------------------------------
# Net stoichiometric matrix STOICH (15 species x 17 reactions).
# Column n is the net change in each species per turnover of reaction R(n+1),
# transcribed directly from the original dt_ode().  dtheta/dt = STOICH @ Rnet.
# --------------------------------------------------------------------------
def _build_stoich(variant="flux") -> np.ndarray:
    species = SPECIES_BLOCKING if variant == "blocking" else SPECIES
    nrxn = 21 if variant == "blocking" else N_REACTIONS
    Sidx = {name: i for i, name in enumerate(species)}
    S = np.zeros((len(species), nrxn))
    rows = {
        "*":    {1: -1, 2: -1, 3: -1, 4: -1, 5: -1, 6: 2, 7: 1, 8: -1,
                 9: 1, 10: 1, 11: -1, 15: -1, 16: 1, 17: 1, 18: 1, 20: 1},
        "CH4*": {1: 1, 2: -1},
        "CH3*": {2: 1, 3: -1},
        "CH2*": {3: 1, 4: -1},
        "CH*":  {4: 1, 5: -1, 7: -1},
        "C*":   {5: 1, 9: -1, 18: -2, 19: 1},
        "H*":   {2: 1, 3: 1, 4: 1, 5: 1, 6: -2, 8: 1, 16: -1, 17: -1},
        "O*":   {7: -1, 9: -1, 11: 1, 15: 1, 19: -1},
        "CO*":  {8: 1, 9: 1, 10: -1, 19: 1},
        "C2*":  {18: 1, 19: -1, 20: -1},
        "CHO*": {7: 1, 8: -1},
        "CO2v": {12: 1, 13: -1, 15: 1},
        "CO3v": {14: 1, 15: -1},
        "OHv":  {16: 1, 17: -1},
        "Ov":   {11: -1, 13: 1, 14: -1, 16: -1},
        "v":    {11: 1, 12: -1, 17: 1},
    }
    if variant == "blocking":
        rows["*"] = {**rows["*"], 20: -1, 21: 1}          # R20b consumes a site, R21 frees one
        rows["Cg*"] = {20: 2, 21: -1}
        rows["O*"] = {**rows["O*"], 21: -1}
        rows["CO*"] = {**rows["CO*"], 21: 1}
    for sp, contrib in rows.items():
        for reaction_number, coeff in contrib.items():   # reaction_number is 1-based
            S[Sidx[sp], reaction_number - 1] = coeff
    return S


SPECIES_BLOCKING = SPECIES[:11] + ["Cg*"] + SPECIES[11:]
STOICH = _build_stoich()


def use_variant(name: str):
    """Switch the module-level network between 'flux' (C2* -> off-lattice
    carbon) and 'blocking' (carbon stays on Ni sites as Cg*, removed by O*)."""
    global VARIANT, SPECIES, REACTIONS, N_SPECIES, N_REACTIONS, N_RATES, N_METAL
    global COV_EXP, P_EXP, COORD_KIND, SITE_KIND, STOICH
    if name == "flux":
        SPECIES = _SPECIES_FLUX; defs = _RATE_DEFS; N_METAL = 11; N_REACTIONS = 20
    elif name == "blocking":
        SPECIES = SPECIES_BLOCKING; defs = _RATE_DEFS_BLOCKING; N_METAL = 12; N_REACTIONS = 21
    else:
        raise ValueError(name)
    VARIANT = name
    REACTIONS = [f"R{i + 1}" for i in range(N_REACTIONS)]
    N_SPECIES = len(SPECIES); N_RATES = 2 * N_REACTIONS
    COV_EXP, P_EXP, COORD_KIND, SITE_KIND = _build_exponents(defs, SPECIES)
    STOICH = _build_stoich(name)

_SPECIES_FLUX = list(SPECIES)


def build_geom(sites: SiteCounts, Zsm: float, Zms: float,
               na: float = 6.022e23) -> np.ndarray:
    """Geometric prefactor for each of the 34 directional rates.

    ``rate_j = k_j * geom_j * prod(theta**COV_EXP_j) * prod(p**P_EXP_j)``.
    ``geom_j = coordination_j * site_count_j / na`` and depends only on the
    catalyst geometry, so it is computed once per simulation.
    """
    coord_val = {"one": 1.0, "Zmm": sites.Zmm, "Zms": Zms, "Zsm": Zsm}
    site_val = {"NSm": sites.NSm, "NSs": sites.NSs,
                "NIm": sites.NIm, "NIs": sites.NIs}
    geom = np.array([coord_val[c] for c in COORD_KIND]) * \
        np.array([site_val[s] for s in SITE_KIND]) / na
    return geom
