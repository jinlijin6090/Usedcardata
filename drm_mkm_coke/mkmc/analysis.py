"""
Reusable analysis primitives shared by the entry-point scripts:

* ``overall_rate``            - scalar objective (default: net CH4 dehydrogenation)
* ``degree_of_rate_control``  - Campbell's DRC for every reaction
* ``evaluate_perturbed``      - solve with a multiplicatively perturbed k-vector
                                (building block for Sobol / parameter sweeps)
"""
from __future__ import annotations

import numpy as np

from .config import ProcessParams
from .network import REACTIONS
from .solver import steady_state
from .rates import rate_vector, net_rates

# Index (0-based) of the reaction whose net rate is the control objective.
# 1 -> R2 (CH4* + * -> CH3* + H*), the first C-H scission = overall DRM flux.
DEFAULT_TARGET = 1


def overall_rate(k: np.ndarray, pp: ProcessParams, target=DEFAULT_TARGET,
                 theta0=None) -> float:
    """Steady-state net rate of the target reaction."""
    t, theta, geom = steady_state(k, pp, theta0=theta0)
    r = rate_vector(theta[-1], k, pp.p, geom)
    return float(net_rates(r)[target])


def degree_of_rate_control(k: np.ndarray, Keq: np.ndarray, pp: ProcessParams,
                           target=DEFAULT_TARGET, dk: float = 0.01):
    """Campbell degree of rate control X_rc,i for each reaction.

    Each forward constant k_i is perturbed by ``dk`` (fractional) while the
    reverse is updated to hold Keq fixed, the steady state is re-solved, and

        X_rc,i = (k_i / r) * dr/dk_i

    Returns
    -------
    labels : list[str]   reaction labels
    Xrc    : np.ndarray  rate-control coefficient per reaction
    """
    base = overall_rate(k, pp, target)
    Xrc = np.zeros(len(REACTIONS))
    for i in range(len(REACTIONS)):
        kp = k.copy()
        f = 2 * i                      # forward index in the interleaved vector
        delta = dk * kp[f]
        kp[f] += delta
        kp[f + 1] = kp[f] / Keq[i]     # hold equilibrium constant
        r_pert = overall_rate(kp, pp, target)
        dr_dk = (r_pert - base) / delta
        Xrc[i] = (k[f] / base) * dr_dk
    return list(REACTIONS), Xrc


def evaluate_perturbed(log10_factors: np.ndarray, k0: np.ndarray,
                       Keq: np.ndarray, pp: ProcessParams,
                       target=DEFAULT_TARGET) -> float:
    """Overall rate when each forward k_i is scaled by 10**factor_i.

    Reverse constants are rescaled to preserve Keq, so only kinetics (not
    thermodynamics) are perturbed.  Used as the model function for global
    (Sobol) sensitivity analysis.
    """
    k = k0.copy()
    scale = 10.0 ** np.asarray(log10_factors, float)
    k[0::2] *= scale
    k[1::2] = k[0::2] / Keq
    try:
        return overall_rate(k, pp, target)
    except Exception:
        return np.nan
