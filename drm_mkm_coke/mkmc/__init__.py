"""Ni/CeO2 DRM microkinetic model (Naveen et al., ACS Catal. 2026) extended
with explicit coking steps (C-C coupling, dimer oxidation, carbon deposition)."""
from .config import ProcessParams, DEFAULT_T
from .parameters import arrhenius_from_excel, RateConstants
from .coke import CokeParams, extend_rate_constants, load_model
from .solver import steady_state, steady_rates
from .rates import rate_vector, net_rates
from .network import SPECIES, REACTIONS, N_METAL
