"""
Geometric site quantification for a Ni / CeO2(111) catalyst.

Unchanged in physics from the original ``calculate_sites`` -- only cleaned up
and given a structured return value so callers read by name instead of by
position.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import ProcessParams


@dataclass
class SiteCounts:
    """Total (catalyst-wide) site counts and coordination numbers."""
    Zmm: float      # nearest Ni neighbours of a Ni site
    Zss: float      # nearest Ov neighbours of an Ov site
    NSm: float      # total surface metal sites
    NSs: float      # total support (Ov) sites
    NIm: float      # total interfacial metal sites
    NIs: float      # total interfacial support sites


def calculate_sites(cat_wt, m_load, rm, ssa, rsm, Zsm, Zms) -> SiteCounts:
    """Quantify metal and support sites for a hemispherical Ni particle.

    Signature kept identical to the original so it is a drop-in replacement;
    prefer :func:`site_counts` when working with a :class:`ProcessParams`.
    """
    na = 6.022e23

    # ---- Ni (metal) site quantification ----
    mass_ni = m_load * cat_wt
    moles_ni = mass_ni / 58.693
    m = moles_ni * na                                   # total Ni atoms
    vol_ni = (2 / 3) * np.pi * (rm ** 3)                # hemispherical particle
    r_fcc = 0.3506 / (2 * (2 ** 0.5))                   # Ni fcc atomic radius
    vol_1ni_fcc = (4 / 3) * np.pi * (r_fcc ** 3)
    n_ni = vol_ni / vol_1ni_fcc                         # atoms per particle
    csa = 2 * np.pi * (rm ** 2)                         # curved surface area
    vol_shell = csa * (2 * r_fcc)
    nem = vol_shell / vol_1ni_fcc                       # surface atoms / particle
    nim = (np.pi * 2 * rm) / (2 * r_fcc)               # interfacial atoms / particle
    zmm = ((nim * 4) + (nem - nim) * 6) / nem           # Ni coordination number
    n_particles = m / n_ni
    nem_net = nem * n_particles                         # NSm
    nim_net = nim * n_particles                         # NIm
    tsa_ni_projected = n_particles * np.pi * ((rm * 1e-9) ** 2)

    # ---- CeO2(111) support site quantification ----
    zss = 3
    areal_density_o = 15.45
    s = (areal_density_o / 1e-18) * ssa                 # NSs (total Ov sites)
    rs = rsm + rm
    nis = np.pi * ((rs ** 2) - (rm ** 2)) * areal_density_o
    nis_net = nis * n_particles                         # NIs

    return SiteCounts(Zmm=zmm, Zss=zss, NSm=nem_net, NSs=s,
                      NIm=nim_net, NIs=nis_net)


def site_counts(pp: ProcessParams) -> SiteCounts:
    """Convenience wrapper that takes a :class:`ProcessParams`."""
    return calculate_sites(*pp.as_site_args())
