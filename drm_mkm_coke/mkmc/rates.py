"""
Vectorised rate laws and ODE right-hand side.

Replaces the 34 hand-written ``r[i] = ...`` lines and the 15 hand-written
``dt[i] = ...`` lines with three array operations, driven by the matrices in
:mod:`mkm.network`.  Results are identical to the original code (see
``tests/test_equivalence.py``).
"""
from __future__ import annotations

import numpy as np

from . import network as net
from .config import NA, ProcessParams
from .sites import calculate_sites


def rate_vector(theta: np.ndarray, k: np.ndarray, p: np.ndarray,
                geom: np.ndarray) -> np.ndarray:
    """All 34 directional rates as a NumPy array.

    Parameters
    ----------
    theta : (15,) surface coverages
    k     : (34,) interleaved forward/reverse rate constants
    p     : (5,)  partial pressures
    geom  : (34,) geometric prefactors from :func:`mkm.network.build_geom`
    """
    theta = np.asarray(theta, dtype=float)
    p = np.asarray(p, dtype=float)
    cov = np.prod(theta ** net.COV_EXP, axis=1)   # 0**0 == 1, integer powers
    pres = np.prod(p ** net.P_EXP, axis=1)
    return k * geom * cov * pres


def net_rates(r: np.ndarray) -> np.ndarray:
    """Forward-minus-reverse rate of each of the 17 reactions."""
    return r[0::2] - r[1::2]


def dtheta_dt(theta: np.ndarray, k: np.ndarray, p: np.ndarray,
              geom: np.ndarray) -> np.ndarray:
    """d(theta)/dt = STOICH @ net_rates."""
    r = rate_vector(theta, k, p, geom)
    return net.STOICH @ net_rates(r)


# --------------------------------------------------------------------------
# Backward-compatible wrappers (same call signature as the original code),
# so the analysis scripts and the equivalence test can call either form.
# --------------------------------------------------------------------------
def rates(theta, k, p, cat_wt, m_load, rm, ssa, rsm, Zsm, Zms):
    """Original-signature wrapper returning the 34 directional rates (list)."""
    sites = calculate_sites(cat_wt, m_load, rm, ssa, rsm, Zsm, Zms)
    geom = net.build_geom(sites, Zsm, Zms, na=NA)
    return list(rate_vector(theta, k, p, geom))


def dt_ode(t, theta, k, p, cat_wt, m_load, rm, ssa, rsm, Zsm, Zms):
    """Original-signature ODE RHS (ignores t; system is autonomous)."""
    sites = calculate_sites(cat_wt, m_load, rm, ssa, rsm, Zsm, Zms)
    geom = net.build_geom(sites, Zsm, Zms, na=NA)
    return list(dtheta_dt(theta, k, p, geom))


def make_geom(pp: ProcessParams) -> np.ndarray:
    """Geometric prefactor vector for a :class:`ProcessParams`."""
    sites = calculate_sites(*pp.as_site_args())
    return net.build_geom(sites, pp.Zsm, pp.Zms, na=NA)
