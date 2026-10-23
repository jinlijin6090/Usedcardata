"""
Energy-space perturbations of the kinetic parameters (Campbell-type sensitivity).

* species n raised by d [eV]  -> every reaction with n as reactant: k_f *= exp(|nu| d/kT)
                                 every reaction with n as product : k_r *= exp(|nu| d/kT)
  (transition-state energies fixed, so barriers shift; Keq shifts consistently)
* transition state i raised by d -> k_f and k_r of reaction i *= exp(-d/kT)
Base-case rate constants are the workbook values, so the base case is unchanged.
"""
import numpy as np
from mkmc import network as net
from mkmc.parameters import RateConstants
from mkmc.config import R_EV

GAS_AND_REF = {"*", "v"}

def perturb(rc: RateConstants, dG_species: dict, dG_ts: dict, T: float) -> RateConstants:
    kT = R_EV * T
    k = rc.k.copy()
    for sp, d in dG_species.items():
        n = net.SPECIES.index(sp)
        for i in range(net.N_REACTIONS):
            s = net.STOICH[n, i]
            if s < 0:
                k[2 * i] *= np.exp(-s * d / kT)
            elif s > 0:
                k[2 * i + 1] *= np.exp(s * d / kT)
    for rx, d in dG_ts.items():
        i = int(rx[1:]) - 1
        k[2 * i] *= np.exp(-d / kT); k[2 * i + 1] *= np.exp(-d / kT)
    with np.errstate(divide="ignore", invalid="ignore"):
        Keq = np.where(k[1::2] > 0, k[0::2] / k[1::2], np.inf)
    return RateConstants(k=k, Keq=Keq, Ea=rc.Ea, k0=rc.k0)

def species_to_vary():
    return [s for s in net.SPECIES if s not in GAS_AND_REF]

def ts_to_vary():
    return list(net.REACTIONS)
