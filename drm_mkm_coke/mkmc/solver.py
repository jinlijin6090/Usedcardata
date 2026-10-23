"""
Steady-state solver for the surface-coverage ODE system.

Integrates with a stiff BDF method on a logarithmic time grid and keeps
extending the horizon until the coverages stop changing.  Same algorithm as
the original ``solve_ode`` but built on the vectorised RHS and a
:class:`ProcessParams`.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .config import NA, DEFAULT_THETA0, ProcessParams
from . import network as net
from .rates import rate_vector, net_rates
from .sites import calculate_sites


def steady_state(k: np.ndarray, pp: ProcessParams,
                 theta0: np.ndarray | None = None,
                 t_final: float = 1e7, n_points: int = 400,
                 tol: float = 1e-7, max_attempts: int = 1000):
    """Integrate to steady state.

    Returns
    -------
    t      : (n_points,) time grid of the final (converged) integration
    theta  : (n_points, 15) coverage trajectory
    geom   : (34,) geometric prefactors (reused by callers for rate analysis)
    """
    if theta0 is None:
        theta0 = np.zeros(len(net.SPECIES))
        theta0[net.SPECIES.index("*")] = 1.0
        theta0[net.SPECIES.index("Ov")] = 0.7
        theta0[net.SPECIES.index("v")] = 0.3
    sites = calculate_sites(*pp.as_site_args())
    geom = net.build_geom(sites, pp.Zsm, pp.Zms, na=NA)

    def rhs(_t, y):
        return net.STOICH @ net_rates(rate_vector(y, k, pp.p, geom))

    t_eval = np.logspace(0, np.log10(t_final), n_points)
    for _ in range(max_attempts):
        res = solve_ivp(rhs, (t_eval[0], t_eval[-1]), theta0,
                        method="BDF", t_eval=t_eval, rtol=1e-6, atol=1e-9)
        if res.y.shape[1] >= 2:
            max_diff = np.max(np.abs(res.y[:, -1] - res.y[:, -2]))
        else:
            max_diff = np.inf
        if max_diff < tol:
            return res.t, res.y.T, geom
        theta0 = np.maximum(res.y[:, -1], 0)
        t_final *= 10
        t_eval = np.logspace(0, np.log10(t_final), n_points)

    raise RuntimeError(f"Coverage did not converge after {max_attempts} attempts.")


def steady_rates(k: np.ndarray, pp: ProcessParams, **kw):
    """Convenience: steady-state coverages, all 34 rates, and 17 net rates."""
    t, theta, geom = steady_state(k, pp, **kw)
    ss = theta[-1]
    r = rate_vector(ss, k, pp.p, geom)
    return ss, r, net_rates(r)
