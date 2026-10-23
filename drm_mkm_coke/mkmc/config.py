"""
Central configuration: file locations, physical constants, and the
catalyst / process parameter container used across the whole package.

Everything that used to be a magic number scattered through the scripts
lives here, so a study is reproduced by editing one place.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# Paths (resolved relative to the project root, so scripts work from anywhere)
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT
OUTPUT_DIR = ROOT / "results"
FIGURE_DIR = ROOT / "figures"
DEFAULT_WORKBOOK = ROOT / "MKM_parameters.xlsx"

for _d in (OUTPUT_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Physical constants
# --------------------------------------------------------------------------
NA = 6.022e23                  # Avogadro number [1/mol]
EV_TO_KJ = 96.49               # 1 eV  ->  kJ/mol
R_EV = 8.314 / (1000 * EV_TO_KJ)   # gas constant in eV/K  (== 8.613e-5)
R_J = R_EV * EV_TO_KJ * 1000   # gas constant in J/mol.K (== 8.311)
DEFAULT_T = 973.15             # operating temperature [K]


@dataclass
class ProcessParams:
    """Catalyst geometry, site coordination and operating conditions.

    Replaces the loose ``cat_wt, m_load, rm, ssa, rsm, Zsm, Zms`` argument
    list that every old function had to thread through by hand.
    """
    p: np.ndarray = field(                 # [pCH4, pCO2, pH2, pCO, pH2O] (bar)
        default_factory=lambda: np.array([0.5, 0.5, 0.0, 0.0, 0.0]))
    cat_wt: float = 1.0                    # total catalyst weight [g]
    m_load: float = 0.04                   # metal (Ni) loading [g per g cat]
    rm: float = 4.0                        # Ni particle radius [nm]
    ssa: float = 70.0                      # specific surface area [m^2/g]
    rsm: float = 0.38                      # interfacial radius [nm]
    Zsm: float = 3.0                       # support->metal coordination
    Zms: float = 2.0                       # metal->support coordination
    T: float = DEFAULT_T                   # temperature [K]

    def as_site_args(self):
        """Tuple in the order expected by :func:`mkm.sites.calculate_sites`."""
        return (self.cat_wt, self.m_load, self.rm, self.ssa,
                self.rsm, self.Zsm, self.Zms)


# Default initial coverage guess: bare metal (*), and a clean support split
# between O-vacancy-free lattice (Ov, theta[13]) and empty support (v, theta[14]).
DEFAULT_THETA0 = np.array(
    [1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0, 0, 0.0, 0.7, 0.3]
)  # 16 species: C2* inserted after CHO* (metal block = first 11)
