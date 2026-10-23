"""
Read kinetic parameters from the Excel workbook and build Arrhenius rate
constants.

Same logic as the original ``arrhenius_from_excel`` but parameterised on the
workbook path and returning a small dataclass instead of four loose arrays.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .config import DEFAULT_T, DEFAULT_WORKBOOK, R_EV


@dataclass
class RateConstants:
    """Interleaved forward/reverse kinetic data (length 34) plus Keq (17)."""
    k: np.ndarray       # rate constants kf,kr,kf,kr,...
    Keq: np.ndarray     # equilibrium constant per reaction (kf/kr)
    Ea: np.ndarray      # activation enthalpies (interleaved)
    k0: np.ndarray      # pre-exponential factors (interleaved)


def arrhenius_from_excel(excel_path: str | Path = DEFAULT_WORKBOOK,
                         sheet: str = "Refined params",
                         T: float = DEFAULT_T) -> RateConstants:
    """Load forward/reverse barriers and prefactors and compute k and Keq.

    The header row is detected by searching for the ``rate index`` cell, so the
    constants block at the top of the sheet is skipped automatically.
    """
    raw = pd.read_excel(excel_path, sheet_name=sheet, header=None)
    header_row = raw[raw.apply(
        lambda r: r.astype(str).str.contains("rate index", case=False).any(),
        axis=1)].index[0]
    df = pd.read_excel(excel_path, sheet_name=sheet, header=header_row)
    df = df[df["rate index"].notna()].reset_index(drop=True)

    Ea_f = pd.to_numeric(df["∆Hf,Ri (eV)"], errors="coerce").to_numpy(float)
    Ea_r = pd.to_numeric(df["∆Hr,Ri (eV)"], errors="coerce").to_numpy(float)
    k0_f = pd.to_numeric(df["Af,Ri (s-1)"], errors="coerce").to_numpy(float)
    k0_r = pd.to_numeric(df["Ar,Ri (s-1)"], errors="coerce").to_numpy(float)

    def interleave(a, b):
        out = np.empty(a.size + b.size)
        out[0::2], out[1::2] = a, b
        return out

    kf = k0_f * np.exp(-Ea_f / (R_EV * T))
    kr = k0_r * np.exp(-Ea_r / (R_EV * T))

    return RateConstants(
        k=interleave(kf, kr),
        Keq=kf / kr,
        Ea=interleave(Ea_f, Ea_r),
        k0=interleave(k0_f, k0_r),
    )
