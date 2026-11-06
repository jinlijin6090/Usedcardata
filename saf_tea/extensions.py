"""Extensions to the replicated TEA:
  1. a China (CN) region,
  2. an EtJ-2 pathway: one-step thermocatalytic CO2 hydrogenation to ethanol
     (2 CO2 + 6 H2 -> C2H5OH + 3 H2O), followed by the same dehydration /
     oligomerisation / hydrogenation / distillation train as EtJ.

Everything that is an assumption (not taken from the paper/SI) is collected in
CN_ASSUMPTIONS and ETJ2_ASSUMPTIONS so it can be changed in one place.
"""
from __future__ import annotations

import csv
import math
from dataclasses import replace

import model as m

# ============================================================== China region
CNY_PER_USD = 7.1128          # PBOC central parity, 19 Sep 2025 (same date as the paper's EUR rate)


def cny(x: float) -> float:
    return x / CNY_PER_USD


CN_ASSUMPTIONS = {
    # electricity: direct green-power purchase price for a hydrogen project in Bayannur,
    # Inner Mongolia (SMM, 2025) + carbon-intensity surcharge as in the paper
    "elec_base_cny_kwh": 0.25,
    "grid_ci_kg_per_kwh": 0.5306,     # 2023 national average CO2 emission factor (MEE/NBS)
    "co2_price_cny_t": 62.36,         # 2025 average CEA price (LSEG)
    "ir": 0.035,                      # 5-year LPR, Sep 2025
    "labor_cny_h": 81.0,              # ~120k CNY/a wage + ~35 % social charges, /2000 h (estimate)
    "location_factor": 0.61,          # China (indigenous) vs US Gulf Coast, Towler & Sinnott
    "el_usd_per_kw": 450.0,           # AEL system, equipment basis -> ~765 USD/kW installed (IEA 2025: 500-1100)
    "co2_feed_cny_t": 200.0,          # high-concentration source (coal chemical / ammonia), estimate
    "water_cny_m3": 15.0,             # demineralised water, estimate
    "waste_cny_m3": 5.0,              # wastewater treatment, estimate
    "cooling_usd_1000m3": 0.21,       # same as US (paper uses US value as default)
    "mp_steam_cny_t": 180.0,
    "lp_steam_cny_t": 160.0,
    "lights_usd_kg": 0.65,            # LPG/propane, ex-tax
    "naphtha_usd_kg": 0.62,
    "diesel_usd_kg": 0.75,
    "wax_usd_kg": 1.20,
}


def china_region(a: dict = CN_ASSUMPTIONS) -> m.Region:
    surcharge_usd_mwh = cny(a["grid_ci_kg_per_kwh"] * a["co2_price_cny_t"])   # kg/kWh * CNY/t = CNY/MWh
    elec = cny(a["elec_base_cny_kwh"] * 1000) + surcharge_usd_mwh
    return m.Region(
        "CN", a["ir"], cny(a["labor_cny_h"]), elec, cny(a["water_cny_m3"]), cny(a["co2_feed_cny_t"]),
        cny(a["mp_steam_cny_t"]), cny(a["lp_steam_cny_t"]), a["lights_usd_kg"], a["naphtha_usd_kg"],
        a["diesel_usd_kg"], a["wax_usd_kg"], a["cooling_usd_1000m3"], cny(a["waste_cny_m3"]),
        capex_factor=a["location_factor"], el_usd_per_kw=a["el_usd_per_kw"],
    )


# ============================================================== EtJ-2
ETJ2_ASSUMPTIONS = {
    # aligned with MtJ-2 (SI Table 9): methanol C / fresh CO2 C = 1514.5/1524.8 = 99.33 %,
    # H2 to synthesis / stoichiometric H2 = 9231.46/9159.7 = 1.0078; no CH4 formation
    "carbon_yield": 1.00,        # share of converted CO2 carbon ending up in ethanol (rest -> CH4, purged)
    "single_pass_conv": 0.20,    # CO2 single-pass conversion of the synthesis loop
    "mtj2_loop_conv": 0.50,      # assumed once-through conversion of the MtJ-2 3-stage loop (scaling ref.)
    "purge_loss": 1 - 1514.5 / 1524.8,     # CO2 lost with the loop purge (as MtJ-2)
    "h2_purge_loss": 1 - 9159.7 / 9231.46, # H2 lost with the loop purge (as MtJ-2)
    "p_syn_bar": 50.0,           # synthesis pressure
    "steam_recovery": 0.5,       # share of reaction heat exported as MP steam
    "capex_safety": 1.20,        # low-TRL contingency (same 20 % the paper puts on EtJ)
    "loop_exp": 0.6,             # scaling exponent, loop equipment
    "feed_comp_exp": 0.8,        # scaling exponent, feed compressors
}

ETOH_KGH = 34948.86          # EtJ ethanol to dehydration (SI Table 21), same SAF output
M_ETOH, M_CO2, M_H2, M_H2O = 46.07, 44.01, 2.016, 18.015
DH_ETOH = 173.6e3            # kJ/kmol, 2CO2+6H2->EtOH(g)+3H2O(g), exothermic
DH_CH4 = 165.0e3             # kJ/kmol, CO2 methanation
EtJ_DOWNSTREAM_H2 = 267.0    # kg/h hydrogenation make-up (calibrated, see model.H2_KGH)
GASFERM_WASTE_M3A = 1288922.96 / 3.52   # SI Table 38 (EU): "Waste - Gasferment" / EU waste price


def etj2_balance(a: dict = ETJ2_ASSUMPTIONS) -> dict:
    n_etoh = ETOH_KGH / M_ETOH                     # kmol/h
    c_etoh = 2 * n_etoh
    c_conv = c_etoh / a["carbon_yield"]
    n_ch4 = c_conv - c_etoh
    co2 = c_conv * M_CO2 / (1 - a["purge_loss"])   # kg/h fresh CO2
    h2_syn = (6 * n_etoh + 4 * n_ch4) * M_H2 / (1 - a.get("h2_purge_loss", a["purge_loss"]))
    h2 = h2_syn + EtJ_DOWNSTREAM_H2
    water_made = (3 * n_etoh + 2 * n_ch4) * M_H2O  # kg/h, goes with crude ethanol to dehydration
    heat_mw = (n_etoh * DH_ETOH + n_ch4 * DH_CH4) / 3600 / 1000
    # compression (isothermal / 0.75, T = 313 K)
    rt = 8.314 * 313
    w_co2 = co2 / M_CO2 / 3600 * rt * math.log(a["p_syn_bar"] / 1.0) / 0.75 / 1000      # MW
    w_h2 = h2_syn / M_H2 / 3600 * rt * math.log(a["p_syn_bar"] / 30.0) / 0.75 / 1000
    fresh_kmol_s = (co2 / M_CO2 + h2_syn / M_H2) / 3600
    circ = fresh_kmol_s / a["single_pass_conv"]
    w_rec = circ * rt * math.log((a["p_syn_bar"] + 5) / a["p_syn_bar"]) / 0.75 / 1000
    return {
        "EtOH_kgh": ETOH_KGH, "CO2_kgh": co2, "H2_kgh": h2, "H2_syn_kgh": h2_syn,
        "water_made_kgh": water_made, "reaction_heat_MW": heat_mw,
        "W_comp_MW": w_co2 + w_h2 + w_rec, "fresh_kgh": co2 + h2_syn,
    }


def load_items():
    items = []
    for r in csv.DictReader(open(m.HERE / "data/equipment_costs.csv")):
        items.append(r)
    return items


def etj2_equipment(a: dict = ETJ2_ASSUMPTIONS) -> dict[str, float]:
    """EtJ-2 equipment cost (USD) per region: EtJ downstream + scaled MtJ-2 synthesis loop."""
    bal = etj2_balance(a)
    mtj2_fresh = 67105.0 + 9231.46                 # kg/h CO2 + H2 to MtJ-2 methanol synthesis
    f_feed = bal["fresh_kgh"] / mtj2_fresh
    f_loop = (bal["fresh_kgh"] / a["single_pass_conv"]) / (mtj2_fresh / a["mtj2_loop_conv"])
    out = {"EU": 0.0, "US": 0.0, "ME": 0.0}
    conv = {"EU": m.USD_PER_EUR, "US": 1.0, "ME": m.USD_PER_SAR}
    col = {"EU": "EU_EUR", "US": "US_USD", "ME": "ME_SAR"}
    for r in load_items():
        item = r["item"]
        if r["route"] == "EtJ" and not item.startswith("GAS-FER"):
            scale = 1.0                                    # downstream identical to EtJ
        elif r["route"] == "MtJ-2" and item.startswith("MEOH-SYN"):
            if "DIST" in item or "COND-DIS" in item or "COOL-DIS" in item or "ME-REF" in item:
                continue                                   # no alcohol distillation: crude EtOH/H2O goes to dehydration
            if item.endswith(("COMP-CO2", "COMP-F")):
                scale = f_feed ** a["feed_comp_exp"]
            else:
                scale = f_loop ** a["loop_exp"]
        else:
            continue
        for reg in out:
            out[reg] += float(r[col[reg]]) * conv[reg] * scale
    return out


def etj2_quantities(q_etj: dict[str, float], a: dict = ETJ2_ASSUMPTIONS,
                    etj_process_mw: float = 0.0) -> dict[str, float]:
    """Annual physical quantities for EtJ-2, derived from EtJ (same downstream)."""
    bal = etj2_balance(a)
    h = m.HOURS
    q = dict(q_etj)
    # electricity: electrolysis + EtJ downstream process power + synthesis compression
    q["electricity"] = (bal["H2_kgh"] * m.EL_KWH_PER_KG / 1000 + etj_process_mw + bal["W_comp_MW"]) * h
    # purified water: electrolysis water + dehydration co-feed minus reaction water already in the crude
    elec_water = 8.93e-3 * bal["H2_kgh"]           # m3/h (EtJ: 8.93 L per kg H2)
    deh_cofeed = max(104846.59 - bal["water_made_kgh"], 0) / 1000
    q["water"] = (elec_water + deh_cofeed) * h
    q["co2"] = bal["CO2_kgh"] / 1000 * h
    # no fermentation broth to treat
    q["waste"] = max(q["waste"] - GASFERM_WASTE_M3A, 0.0)
    q["mp_steam"] = a["steam_recovery"] * bal["reaction_heat_MW"] * 3600 / 2.0 * h / 1000  # t/a, ~2 GJ/t
    return q


def etj_process_mw(vo_si) -> float:
    """Non-electrolysis power of EtJ, MW (from EU electricity bill)."""
    e_mwh = vo_si[("EtJ", "EU")]["electricity"] / m.REGIONS["EU"].elec
    return e_mwh / m.HOURS - m.H2_KGH["EtJ"] * m.EL_KWH_PER_KG / 1000
