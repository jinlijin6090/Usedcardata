"""
Coking extension: parameters for R18-R20 and the function that appends them
to the 34 published rate constants.

    R18  C* + C*  <->  C2* + *       C-C coupling (coke nucleation)
    R19  C2* + O* <->  CO* + C*      oxidative scission of the dimer
    R20  C2*       ->  C(dep) + *    irreversible growth into graphitic /
                                     filamentous carbon (off-lattice sink)

Thermodynamic consistency: R18 + R19 == R9 (C* + O* -> CO* + *), so
Keq19 = Keq9 / Keq18 is enforced rather than chosen.

Parameter provenance (see report): DFT literature for C+C coupling on Ni(111)
spans roughly 1.0-1.6 eV for the barrier with reaction energies between about
-0.5 and +0.5 eV depending on facet and coverage; graphitic growth / carbon
dissolution barriers reported ~1.0-1.5 eV. These are therefore treated as
scan variables, with the nominal values below as the base case.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .config import DEFAULT_T, DEFAULT_WORKBOOK, R_EV
from .parameters import RateConstants, arrhenius_from_excel

KBT_H = 2.0268e13  # kB*T/h at 973.15 K [1/s], same value as the workbook


@dataclass
class CokeParams:
    Ea_cc: float = 1.30      # R18 forward barrier [eV]
    dH_cc: float = -0.30     # R18 reaction enthalpy [eV]  (C2* vs 2 C*)
    dS_cc: float = 0.0       # R18 reaction entropy [eV/K] (ln(Af/Ar)*kB)
    Ea_ox: float = 1.40      # R19 forward barrier [eV]
    Ea_dep: float = 1.20     # R20 barrier [eV]
    A_cc: float = KBT_H      # prefactors [1/s]
    A_ox: float = KBT_H
    A_dep: float = KBT_H
    # blocking variant only
    dH_g: float = -0.50      # C* -> Cg* enthalpy per C atom [eV] (graphitic C more stable)
    Ea_gas: float = 1.80     # R21 Cg* + O* -> CO* + * barrier [eV]
    A_gas: float = KBT_H

    def with_(self, **kw):
        return replace(self, **kw)


def extend_rate_constants(rc: RateConstants, cp: CokeParams,
                          T: float = DEFAULT_T) -> RateConstants:
    """Append k18..k20 (interleaved f/r) to the published 34-vector."""
    kT = R_EV * T
    # R18
    kf18 = cp.A_cc * np.exp(-cp.Ea_cc / kT)
    Keq18 = np.exp(-(cp.dH_cc - T * cp.dS_cc) / kT)
    kr18 = kf18 / Keq18
    # R19, tied to the R9 / R18 cycle
    Keq19 = rc.Keq[8] / Keq18
    kf19 = cp.A_ox * np.exp(-cp.Ea_ox / kT)
    kr19 = kf19 / Keq19
    # R20, irreversible sink
    kf20 = cp.A_dep * np.exp(-cp.Ea_dep / kT)
    kr20 = 0.0

    k = np.concatenate([rc.k, [kf18, kr18, kf19, kr19, kf20, kr20]])
    Keq = np.concatenate([rc.Keq, [Keq18, Keq19, np.inf]])
    Ea = np.concatenate([rc.Ea, [cp.Ea_cc, cp.Ea_cc - cp.dH_cc,
                                 cp.Ea_ox, cp.Ea_ox + kT * np.log(Keq19),
                                 cp.Ea_dep, np.inf]])
    k0 = np.concatenate([rc.k0, [cp.A_cc, cp.A_cc, cp.A_ox, cp.A_ox,
                                 cp.A_dep, 0.0]])
    from . import network as net
    if net.VARIANT == "blocking":
        # R21 tied to the cycle  R21 = R9 - (C* -> Cg*)
        Kg = np.exp(-cp.dH_g / kT)
        Keq21 = rc.Keq[8] / Kg
        kf21 = cp.A_gas * np.exp(-cp.Ea_gas / kT)
        kr21 = kf21 / Keq21
        k = np.concatenate([k, [kf21, kr21]]); Keq = np.concatenate([Keq, [Keq21]])
        Ea = np.concatenate([Ea, [cp.Ea_gas, cp.Ea_gas + kT * np.log(Keq21)]])
        k0 = np.concatenate([k0, [cp.A_gas, cp.A_gas]])
    return RateConstants(k=k, Keq=Keq, Ea=Ea, k0=k0)


def load_model(cp: CokeParams | None = None, T: float = DEFAULT_T,
               workbook=DEFAULT_WORKBOOK, coke: bool = True) -> RateConstants:
    """Published parameters (+ coking steps unless coke=False, in which case
    the three extra steps are present but switched off with k=0)."""
    rc = arrhenius_from_excel(workbook, "Refined params", T)
    cp = cp or CokeParams()
    ext = extend_rate_constants(rc, cp, T)
    if not coke:
        ext.k[34:] = 0.0  # also switches off R21 in the blocking variant
    return ext
