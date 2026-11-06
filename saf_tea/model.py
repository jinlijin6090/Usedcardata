"""Net-production-cost (NPC) model reproducing the TEA of
Sendlhofer et al., Fuel 421 (2026) 139021, extended with a China region
and an EtJ-2 (thermocatalytic CO2-to-ethanol) pathway.

All results in USD (2025). Equations follow the paper:
  ACC = (FCI+WC)*CRF(IR,PL) - WC/PL                      (Eq. 12)
  SR  = floor(PL*h/SL - 1);  ARC = SR*CAPEX_EL*0.77/PL    (Eq. 13-14)
  h_labor = 2.13*PC^0.242*n*h/24                         (Eq. 15)
  NPC = (ACC + ARC + fixOpex + varOpex - BPR)/m_SAF      (Eq. 16; BPR inside varOpex)
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).parent
EUR_PER_USD = 0.8449          # paper: USD = 0.8449 EUR
USD_PER_EUR = 1 / EUR_PER_USD
USD_PER_SAR = 1 / 3.75        # SAR is pegged to USD

# ---------------------------------------------------------------- plant basis
HOURS = 7884.0                # h/a, implied by SI varOpex tables (CO2 cost / CO2 feed)
SAF_KGH = 13700.0             # kg/h SAF in every route (SI stream tables)
M_SAF = SAF_KGH * HOURS       # kg/a  (= 108 kt/a; the paper quotes 120 kt nameplate)
PL = 20                       # plant lifetime, a
STACK_LIFE = 60000.0          # h
EL_KWH_PER_KG = 55.0          # kWh_el per kg H2 (AEL)
EL_EUR_PER_KW = 800.0         # AEL system cost, EUR/kW (paper)

# H2 from electrolysis, kg/h (SI stream tables, streams with origin "AEL").
# EtJ hydrogenation make-up is not given net in the SI; it is calibrated so that
# the electrolyser CAPEX matches Fig. 7a (~880 MUSD).
H2_KGH = {
    "FTS": 11754.25 + 297.44,
    "MtJ-1": 10290.28 + 211.80,
    "MtJ-2": 9231.46 + 217.95,
    "EtJ": 9690.72 + 267.0,
}

# Process steps n (SI Table 34) - count of non-zero entries incl. AEL and BoP
N_STEPS = {"FTS": 6, "MtJ-1": 8, "MtJ-2": 7, "EtJ": 7, "EtJ-2": 7}

# Surcharge factors (SI Tables 31-33)
DIRECT = 0.47 + 0.36 + 0.68 + 0.11 + 0.24 + 0.10 + 0.55
INDIRECT = 0.33 + 0.41 + 0.04
FCI_FACTOR = (1 + DIRECT + INDIRECT) * (1 + 0.05 + 0.10)            # = 4.93 x EC
FCI_EL_FACTOR = (1 + 0.47 + 0.10 + 0.04) + 0.15 * (0.47 + 0.10 + 0.04)  # = 1.70 x EC_h
ETJ_SAFETY = 1.20             # +20 % on all EtJ CAPEX (paper Sect. 2.2.1)


# ---------------------------------------------------------------- regions
@dataclass
class Region:
    name: str
    ir: float                       # interest rate
    labor_usd_h: float
    elec: float                     # USD/MWh, already incl. carbon-intensity surcharge
    water: float                    # USD/m3 purified feed water
    co2: float                      # USD/t CO2 feedstock
    mp_steam: float                 # USD/t (credit)
    lp_steam: float
    lights: float                   # USD/kg (credit)
    naphtha: float
    diesel: float
    wax: float
    cooling: float                  # USD per 1000 m3
    waste: float                    # USD/m3
    capex_factor: float = 1.0       # multiplier on the EU equipment list (US/ME use own lists)
    el_usd_per_kw: float = EL_EUR_PER_KW * USD_PER_EUR
    notes: dict = field(default_factory=dict)


def eur(x: float) -> float:
    return x * USD_PER_EUR


REGIONS = {
    "EU": Region("EU", 0.0392, eur(40.0), eur(67.70), eur(9.56), eur(110), eur(25.7), eur(26.3),
                 eur(1.07), eur(0.59), eur(0.98), eur(1.88), eur(0.18), eur(3.52)),
    "US": Region("US", 0.0432, 38.20, 32.81, 9.94, 127.24, 30.96, 30.25,
                 0.65, 0.54, 0.69, 1.45, 0.21, 3.67),
    "ME": Region("ME", 0.0422, 27.33, 98.01, 22.17, 127.34, 30.96, 30.25,
                 0.66, 0.54, 0.90, 1.22, 0.409, 5.0),
}


# ---------------------------------------------------------------- data loading
def load_equipment() -> dict[str, dict[str, float]]:
    """Total equipment cost per route and region, USD."""
    ec: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for r in csv.DictReader(open(HERE / "data/equipment_costs.csv")):
        ec[r["route"]]["EU"] += float(r["EU_EUR"]) * USD_PER_EUR
        ec[r["route"]]["US"] += float(r["US_USD"])
        ec[r["route"]]["ME"] += float(r["ME_SAR"]) * USD_PER_SAR
    return ec


def load_varopex_si() -> dict[tuple[str, str], dict[str, float]]:
    """SI varOpex line items, converted to USD (EU tables are in EUR)."""
    out: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    for r in csv.DictReader(open(HERE / "data/varopex_si.csv")):
        k = (r["route"], r["region"])
        v = float(r["value"]) * (USD_PER_EUR if r["region"] == "EU" else 1.0)
        cat = canonical(r["category"])
        out[k][cat] = out[k].get(cat, 0.0) + v
    return out


def canonical(cat: str) -> str:
    c = cat.lower()
    for key, words in [
        ("TOTAL", ["variable opex"]), ("electricity", ["electric"]), ("water", ["water"]),
        ("co2", ["co2"]), ("lights", ["light"]), ("naphtha", ["naph"]), ("diesel", ["diesel"]),
        ("wax", ["wax"]), ("lp_steam", ["lp"]), ("mp_steam", ["mp"]), ("oxygen", ["oxygen"]),
    ]:
        if any(w in c for w in words):
            if key == "water" and ("waste" in c or "cool" in c):
                continue
            return key
    if "cool" in c:
        return "cooling"
    if "waste" in c:
        return "waste"
    return c


def quantities_from_eu(vo_si) -> dict[str, dict[str, float]]:
    """Back out physical annual quantities from the EU varOpex table and EU prices."""
    eu = REGIONS["EU"]
    price = {"electricity": eu.elec, "water": eu.water, "co2": eu.co2, "lights": eu.lights,
             "naphtha": eu.naphtha, "diesel": eu.diesel, "wax": eu.wax, "lp_steam": eu.lp_steam,
             "mp_steam": eu.mp_steam, "cooling": eu.cooling, "waste": eu.waste}
    q: dict[str, dict[str, float]] = {}
    for route in ["FTS", "MtJ-1", "MtJ-2", "EtJ"]:
        items = vo_si[(route, "EU")]
        q[route] = {k: abs(items.get(k, 0.0)) / p for k, p in price.items()}
    return q


def varopex_from_quantities(q: dict[str, float], reg: Region) -> dict[str, float]:
    costs = {
        "electricity": q["electricity"] * reg.elec,
        "water": q["water"] * reg.water,
        "co2": q["co2"] * reg.co2,
        "waste": q["waste"] * reg.waste,
        "cooling": q["cooling"] * reg.cooling,
        "lights": -q["lights"] * reg.lights,
        "naphtha": -q["naphtha"] * reg.naphtha,
        "diesel": -q["diesel"] * reg.diesel,
        "wax": -q["wax"] * reg.wax,
        "lp_steam": -q["lp_steam"] * reg.lp_steam,
        "mp_steam": -q["mp_steam"] * reg.mp_steam,
    }
    costs["TOTAL"] = sum(costs.values())
    return costs


# ---------------------------------------------------------------- core economics
def crf(ir: float, n: int) -> float:
    return ir * (1 + ir) ** n / ((1 + ir) ** n - 1)


def npc(route: str, reg: Region, ec_usd: float, h2_kgh: float, varopex_usd: float,
        fixopex_lab_basis: str = "OLC") -> dict[str, float]:
    """Return NPC breakdown in USD/kg SAF."""
    safety = ETJ_SAFETY if route.startswith("EtJ") else 1.0
    fci_jet = ec_usd * FCI_FACTOR * safety
    p_el_kw = h2_kgh * EL_KWH_PER_KG
    ec_el = p_el_kw * reg.el_usd_per_kw
    fci_el = ec_el * FCI_EL_FACTOR
    fci = fci_jet + fci_el
    wc = fci / 0.9 * 0.1                      # WC = 10 % of TCI
    acc = (fci + wc) * crf(reg.ir, PL) - wc / PL
    sr = math.floor(PL * HOURS / STACK_LIFE - 1)
    arc = sr * fci_el * 0.77 / PL
    # fixed OPEX (SI Table 33), applied to the jet-process FCI
    h_lab = 2.13 * SAF_KGH ** 0.242 * N_STEPS[route] * HOURS / 24
    olc = h_lab * reg.labor_usd_h
    os_ = 0.15 * olc
    mt = 0.02 * fci_jet
    mat = 0.02 * fci_jet
    sup = 0.15 * olc
    lab = 0.2 * (ec_usd * safety if fixopex_lab_basis == "EC" else olc)
    ins = 0.02 * fci_jet
    po = 0.6 * (olc + os_ + mt)
    adm = 0.25 * po
    fix = olc + os_ + mt + mat + sup + lab + ins + po + adm
    res = {
        "FCI_jet_MUSD": fci_jet / 1e6, "FCI_el_MUSD": fci_el / 1e6, "FCI_MUSD": fci / 1e6,
        "ACC": acc / M_SAF, "ARC": arc / M_SAF, "fixOpex": fix / M_SAF,
        "varOpex": varopex_usd / M_SAF,
    }
    res["NPC"] = res["ACC"] + res["ARC"] + res["fixOpex"] + res["varOpex"]
    return res
