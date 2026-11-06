"""China-only routes and the China sensitivity study.

New routes (same SAF output, 13.7 t/h):
  HEFA   - hydroprocessing of waste lipids (UCO) in max-jet mode, H2 from on-site AEL.
  DRM-FT - dry reforming of methane (CH4 + CO2 -> 2 CO + 2 H2) replacing the rWGS of
           the paper's FTS route; FT / hydrocracking / distillation unchanged;
           electrolytic H2 only tops up the H2/CO ratio.

All route assumptions sit in HEFA_ASSUMPTIONS / DRM_ASSUMPTIONS; China market data in
extensions.CN_ASSUMPTIONS (plus the keys added below).
"""
from __future__ import annotations

import copy
import csv

import matplotlib.pyplot as plt
from openpyxl import load_workbook

import extensions as x
import model as m
import run_all as ra   # sets up fonts, equipment & varOpex data

OUT = ra.OUT
cny = x.cny

# extra China market data used by the new routes
x.CN_ASSUMPTIONS.update({
    "uco_cny_t": 7500.0,        # waste cooking oil, ~1,050 USD/t (2024 FOB 910-970 + loss of 13 % export rebate)
    "ch4_cny_nm3": 3.0,         # industrial pipeline gas / biomethane, ~2.5-3.5 CNY/Nm3 in 2025
    "hours": m.HOURS,           # operating hours (sensitivity only)
})

HEFA_ASSUMPTIONS = {
    # max-jet mode, per t lipid feed (typical literature ranges, e.g. Pearlson et al. 2013)
    "saf_yield": 0.50,          # t SAF / t feed  (0.45-0.55)
    "rd_yield": 0.17,           # renewable diesel, credited at diesel price
    "naphtha_yield": 0.07,
    "lpg_yield": 0.06,          # propane from glycerol backbone + light ends
    "h2_per_t": 0.040,          # t H2 / t feed (0.035-0.045)
    "elec_kwh_per_t": 60.0,     # process power, kWh / t feed
    "fuel_gj_per_t": 1.2,       # process heat (natural gas), GJ / t feed
    "water_m3_per_t": 0.5,
    "waste_m3_per_t": 0.3,
    "fci_usd_per_tpa": 600.0,   # China FCI per t/a of total liquid product (~4,300 CNY/t; 50 万 t/a ~ 21 亿元)
}

DRM_ASSUMPTIONS = {
    "ch4_to_co2": 1.0,          # molar CH4 : CO2 in fresh feed (1 = stoichiometric DRM)
    "carbon_eff": (13700 / 14.027) / (84414.5 / 44.01),  # SAF C / feed C, as in the paper's FTS route
    "reformer_ec_usd": 15e6,    # extra equipment for the electrically heated DRM reformer (US basis)
    "capex_safety": 1.20,       # DRM catalyst coking / TRL 5-6 contingency
    "dh_drm_kj_mol": 247.0,     # CH4 + CO2 -> 2CO + 2H2
}

NM3_PER_KMOL = 22.414
FTS_CO2_KGH = 84414.50                  # SI Table 1, fresh CO2 of the FTS route
SAF_C_KMOLH = m.SAF_KGH / 14.027        # SAF approximated as (CH2)n


# ---------------------------------------------------------------- HEFA
def hefa(cn: dict, a: dict = HEFA_ASSUMPTIONS, hours: float | None = None) -> dict:
    hours = hours or cn["hours"]
    reg = x.china_region(cn)
    feed_tph = m.SAF_KGH / 1000 / a["saf_yield"]
    h2_kgh = feed_tph * a["h2_per_t"] * 1000
    liquids_tpa = feed_tph * (a["saf_yield"] + a["rd_yield"] + a["naphtha_yield"] + a["lpg_yield"]) * m.HOURS
    fci = a["fci_usd_per_tpa"] * liquids_tpa
    gas_nm3h = feed_tph * a["fuel_gj_per_t"] / 0.0359      # 35.9 MJ/Nm3 LHV
    h = m.HOURS
    vo = {
        "feed": feed_tph * h * cny(cn["uco_cny_t"]),
        "electricity": (h2_kgh * m.EL_KWH_PER_KG / 1000 + feed_tph * a["elec_kwh_per_t"] / 1000) * h * reg.elec,
        "fuel": gas_nm3h * h * cny(cn["ch4_cny_nm3"]),
        "water": (feed_tph * a["water_m3_per_t"] + 8.93e-3 * h2_kgh) * h * reg.water,
        "waste": feed_tph * a["waste_m3_per_t"] * h * reg.waste,
        "co2": 0.0,
        "diesel": -feed_tph * a["rd_yield"] * 1000 * h * reg.diesel,
        "naphtha": -feed_tph * a["naphtha_yield"] * 1000 * h * reg.naphtha,
        "lights": -feed_tph * a["lpg_yield"] * 1000 * h * reg.lights,
    }
    vo["TOTAL"] = sum(vo.values())
    res = m.npc("HEFA", reg, 0.0, h2_kgh, vo["TOTAL"], hours=hours, fci_jet_direct=fci)
    res.update(route="HEFA", H2_kgh=h2_kgh, power_cost=vo["electricity"] / m.M_SAF,
               feed_cost=vo["feed"] / m.M_SAF, co2_cost=0.0, feed_tph=feed_tph)
    return res


# ---------------------------------------------------------------- DRM-FT
def drm_balance(a: dict = DRM_ASSUMPTIONS) -> dict:
    c_feed = SAF_C_KMOLH / a["carbon_eff"]                 # kmol C/h fresh
    r = a["ch4_to_co2"]
    n_ch4 = c_feed * r / (1 + r)
    n_co2 = c_feed - n_ch4
    # H2 balance relative to the paper's FTS route (CO2 only): per C need 1 H2 (CH2) + 2 H2 per CO2 (O removal),
    # CH4 brings 2 H2. Extra consumption (purge, paraffin end groups) taken from the FTS route.
    fts_c = FTS_CO2_KGH / 44.01
    extra = m.H2_KGH["FTS"] / m_h2() - 3 * fts_c
    h2_kmol = c_feed + 2 * n_co2 - 2 * n_ch4 + extra
    heat_mw = (n_ch4 * a["dh_drm_kj_mol"] - fts_c * 41.0) / 3600
    return {"CH4_kmolh": n_ch4, "CO2_kmolh": n_co2, "CO2_kgh": n_co2 * 44.01,
            "CH4_nm3h": n_ch4 * NM3_PER_KMOL, "H2_kgh": max(h2_kmol, 0) * m_h2(),
            "extra_heat_MW": heat_mw}


def m_h2() -> float:
    return 2.016


def drm_ft(cn: dict, a: dict = DRM_ASSUMPTIONS, hours: float | None = None) -> dict:
    hours = hours or cn["hours"]
    reg = x.china_region(cn)
    bal = drm_balance(a)
    q = dict(ra.Q["FTS"])
    fts_process_mw = q["electricity"] / m.HOURS - m.H2_KGH["FTS"] * m.EL_KWH_PER_KG / 1000
    q["electricity"] = (bal["H2_kgh"] * m.EL_KWH_PER_KG / 1000 + fts_process_mw + bal["extra_heat_MW"]) * m.HOURS
    q["water"] = q["water"] - 8.93e-3 * (m.H2_KGH["FTS"] - bal["H2_kgh"]) * m.HOURS
    q["co2"] = bal["CO2_kgh"] / 1000 * m.HOURS
    vo = m.varopex_from_quantities(q, reg)
    ch4_cost = bal["CH4_nm3h"] * m.HOURS * cny(cn["ch4_cny_nm3"])
    vo["ch4"] = ch4_cost
    vo["TOTAL"] += ch4_cost
    ec = (ra.EC["FTS"]["US"] + a["reformer_ec_usd"]) * reg.capex_factor
    res = m.npc("DRM-FT", reg, ec, bal["H2_kgh"], vo["TOTAL"], hours=hours, safety=a["capex_safety"])
    res.update(route="DRM-FT", H2_kgh=bal["H2_kgh"], power_cost=vo["electricity"] / m.M_SAF,
               feed_cost=ch4_cost / m.M_SAF, co2_cost=vo["co2"] / m.M_SAF)
    return res


# ---------------------------------------------------------------- unified China scenario
ROUTES_CN = ["FTS", "MtJ-1", "MtJ-2", "EtJ", "EtJ-2", "DRM-FT", "HEFA"]


def cn_scenario(route: str, cn: dict | None = None, e2: dict | None = None,
                hefa_a: dict | None = None, drm_a: dict | None = None) -> dict:
    cn = cn or x.CN_ASSUMPTIONS
    if route == "HEFA":
        return hefa(cn, hefa_a or HEFA_ASSUMPTIONS)
    if route == "DRM-FT":
        return drm_ft(cn, drm_a or DRM_ASSUMPTIONS)
    reg = x.china_region(cn)
    if route == "EtJ-2":
        e2 = e2 or x.ETJ2_ASSUMPTIONS
        ec = x.etj2_equipment(e2)["US"]
        h2 = x.etj2_balance(e2)["H2_kgh"]
        vo = m.varopex_from_quantities(x.etj2_quantities(ra.Q["EtJ"], e2, ra.ETJ_PMW), reg)
        safety = e2["capex_safety"]
    else:
        ec, h2 = ra.EC[route]["US"], m.H2_KGH[route]
        vo = m.varopex_from_quantities(ra.Q[route], reg)
        safety = None
    res = m.npc(route, reg, ec * reg.capex_factor, h2, vo["TOTAL"], hours=cn["hours"], safety=safety)
    res.update(route=route, H2_kgh=h2, power_cost=vo["electricity"] / m.M_SAF,
               co2_cost=vo["co2"] / m.M_SAF, feed_cost=0.0)
    return res


# ---------------------------------------------------------------- sensitivity
COMMON = [  # key, label, low, high
    ("elec_base_cny_kwh", "绿电价格 0.15–0.40 元/kWh", 0.15, 0.40),
    ("el_usd_per_kw", "电解槽 250–947 USD/kW", 250, 947),
    ("hours", "年运行小时 5000–8500 h", 5000, 8500),
    ("co2_feed_cny_t", "CO₂ 原料 100–905 元/t", 100, 905),
    ("location_factor", "设备地区系数 0.5–1.0", 0.5, 1.0),
    ("ir", "贷款利率 2.5–5.0%", 0.025, 0.05),
    ("co2_price_cny_t", "碳价 40–200 元/t", 40, 200),
    ("uco_cny_t", "地沟油 5,500–9,500 元/t", 5500, 9500),
    ("ch4_cny_nm3", "天然气/生物甲烷 2.0–5.0 元/Nm³", 2.0, 5.0),
]
SPECIFIC = {
    "HEFA": [("hefa", "fci_usd_per_tpa", "HEFA 投资 400–1,200 USD/(t/a)", 400, 1200),
             ("hefa", "saf_yield", "HEFA 航煤收率 0.45–0.55", 0.45, 0.55)],
    "DRM-FT": [("drm", "reformer_ec_usd", "重整炉设备 5–40 MUSD", 5e6, 40e6),
               ("drm", "ch4_to_co2", "CH₄:CO₂ 0.8–1.5", 0.8, 1.5)],
    "EtJ-2": [("e2", "carbon_yield", "乙醇碳选择性 0.65–0.95", 0.65, 0.95)],
}


def run_sensitivity():
    rows = []
    base = {r: cn_scenario(r)["NPC"] for r in ROUTES_CN}
    for route in ROUTES_CN:
        for key, label, lo, hi in COMMON:
            vals = [cn_scenario(route, cn=dict(x.CN_ASSUMPTIONS, **{key: v}))["NPC"] for v in (lo, hi)]
            rows.append({"route": route, "param": label, "NPC_base": base[route],
                         "NPC_low": vals[0], "NPC_high": vals[1], "swing": abs(vals[1] - vals[0])})
        for grp, key, label, lo, hi in SPECIFIC.get(route, []):
            vals = []
            for v in (lo, hi):
                kw = {"e2": None, "hefa_a": None, "drm_a": None}
                src = {"e2": x.ETJ2_ASSUMPTIONS, "hefa": HEFA_ASSUMPTIONS, "drm": DRM_ASSUMPTIONS}[grp]
                arg = {"e2": "e2", "hefa": "hefa_a", "drm": "drm_a"}[grp]
                kw[arg] = dict(src, **{key: v})
                vals.append(cn_scenario(route, **kw)["NPC"])
            rows.append({"route": route, "param": label, "NPC_base": base[route],
                         "NPC_low": vals[0], "NPC_high": vals[1], "swing": abs(vals[1] - vals[0])})
    return base, rows


def sweep(key, values):
    return {r: [cn_scenario(r, cn=dict(x.CN_ASSUMPTIONS, **{key: v}))["NPC"] for v in values] for r in ROUTES_CN}


def crossover(route_a, route_b, key, lo, hi):
    """Value of `key` at which route_a and route_b have equal NPC (bisection), or None."""
    f = lambda v: (cn_scenario(route_a, cn=dict(x.CN_ASSUMPTIONS, **{key: v}))["NPC"]
                   - cn_scenario(route_b, cn=dict(x.CN_ASSUMPTIONS, **{key: v}))["NPC"])
    if f(lo) * f(hi) > 0:
        return None
    for _ in range(50):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


# ---------------------------------------------------------------- plots
STYLE = {"FTS": ("#9AA5B1", "-"), "MtJ-1": ("#5B8DD6", "-"), "MtJ-2": ("#E8A33D", "-"),
         "EtJ": ("#2F4A6B", "-"), "EtJ-2": ("#2F4A6B", "--"), "DRM-FT": ("#4C9A2A", "-"),
         "HEFA": ("#C0392B", "-")}


def plot_bars(res):
    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=150)
    for i, r in enumerate(ROUTES_CN):
        d = res[r]
        parts = [("原料(地沟油/CH₄)", d.get("feed_cost", 0), "#8E6C8A"),
                 ("其他可变成本", d["varOpex"] - d.get("feed_cost", 0), ra.COLORS["varOpex"]),
                 ("fixOpex", d["fixOpex"], ra.COLORS["fixOpex"]), ("ARC", d["ARC"], ra.COLORS["ARC"]),
                 ("ACC", d["ACC"], ra.COLORS["ACC"])]
        bottom = 0
        for lab, v, c in parts:
            ax.bar(i, v, bottom=bottom, color=c, edgecolor="#333", lw=0.6, label=lab if i == 0 else None)
            bottom += v
        ax.text(i, bottom + 0.06, f"{bottom:.2f}", ha="center", fontsize=10)
    ax.axhline(1.036, ls="--", color="k", lw=1.2, label="化石 Jet A-1")
    ax.set_xticks(range(len(ROUTES_CN))); ax.set_xticklabels(ROUTES_CN)
    ax.set_ylabel("NPC (USD/kg SAF)"); ax.set_title("中国：7 条路线的净生产成本（基准情景）")
    ax.legend(fontsize=8, ncol=3, loc="upper left"); ax.grid(axis="y", ls=":", alpha=0.5)
    ax.set_ylim(0, max(res[r]["NPC"] for r in ROUTES_CN) * 1.3)
    fig.tight_layout(); fig.savefig(OUT / "cn_7routes_npc.png"); plt.close(fig)


def plot_tornado_grid(rows, base):
    fig, axes = plt.subplots(2, 4, figsize=(20, 9), dpi=130)
    for ax, route in zip(axes.flat, ROUTES_CN):
        rr = [r for r in rows if r["route"] == route and r["swing"] > 1e-3]
        rr = sorted(rr, key=lambda r: r["swing"])[-7:]
        b = base[route]
        for i, r in enumerate(rr):
            for v, c in ((r["NPC_low"], "#3B6EA8"), (r["NPC_high"], "#E8A33D")):
                ax.barh(i, v - b, left=b, color=c, edgecolor="#333", lw=0.4)
        ax.axvline(b, color="k", lw=1)
        ax.set_yticks(range(len(rr))); ax.set_yticklabels([r["param"] for r in rr], fontsize=8)
        ax.set_title(f"{route}（基准 {b:.2f}）", fontsize=11)
        ax.set_xlim(0.8, 4.2); ax.grid(axis="x", ls=":", alpha=0.5)
    axes.flat[-1].axis("off")
    axes.flat[-1].text(0.05, 0.6, "蓝 = 参数取区间下限\n橙 = 参数取区间上限\n横轴：NPC (USD/kg)\n每图列出影响最大的 7 项",
                       fontsize=12, va="top")
    fig.suptitle("中国 7 条路线的单因素敏感性", fontsize=14)
    fig.tight_layout(); fig.savefig(OUT / "cn_tornado_grid.png"); plt.close(fig)


def plot_sweep(key, values, xlabel, fname, title, scale=1.0):
    data = sweep(key, values)
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=150)
    for r in ROUTES_CN:
        c, ls = STYLE[r]
        ax.plot([v * scale for v in values], data[r], color=c, ls=ls, lw=2, label=r)
    ax.axhline(1.036, ls=":", color="k", lw=1)
    ax.text(values[0] * scale, 1.08, "化石 Jet A-1", fontsize=8)
    base_v = x.CN_ASSUMPTIONS[key] * scale
    ax.axvline(base_v, color="#999", lw=1, ls="--")
    ax.set_xlabel(xlabel); ax.set_ylabel("NPC (USD/kg SAF)"); ax.set_title(title)
    ax.legend(ncol=2, fontsize=9); ax.grid(ls=":", alpha=0.5)
    fig.tight_layout(); fig.savefig(OUT / fname); plt.close(fig)
    return data


def main():
    res = {r: cn_scenario(r) for r in ROUTES_CN}
    base, rows = run_sensitivity()
    plot_bars(res); plot_tornado_grid(rows, base)
    elec_vals = [round(0.10 + 0.025 * i, 3) for i in range(15)]
    hours_vals = list(range(3000, 8761, 250))
    plot_sweep("elec_base_cny_kwh", elec_vals, "绿电价格（元/kWh）", "cn_sweep_electricity.png",
               "中国：NPC 随绿电价格变化")
    plot_sweep("hours", hours_vals, "年运行小时（h）", "cn_sweep_hours.png", "中国：NPC 随年运行小时变化")
    plot_sweep("co2_feed_cny_t", list(range(0, 1001, 50)), "CO₂ 原料价（元/t）", "cn_sweep_co2.png",
               "中国：NPC 随 CO₂ 原料价变化")
    cross = {
        "MtJ-2 追平 HEFA 的绿电价 (元/kWh)": crossover("MtJ-2", "HEFA", "elec_base_cny_kwh", 0.0, 0.6),
        "DRM-FT 追平 MtJ-2 的天然气价 (元/Nm³)": crossover("DRM-FT", "MtJ-2", "ch4_cny_nm3", 0.0, 8.0),
        "HEFA 追平 MtJ-2 的地沟油价 (元/t)": crossover("HEFA", "MtJ-2", "uco_cny_t", 2000, 15000),
        "DRM-FT 追平 HEFA 的天然气价 (元/Nm³)": crossover("DRM-FT", "HEFA", "ch4_cny_nm3", 0.0, 8.0),
    }

    tab = [{"路线": r, "NPC": d["NPC"], "varOpex": d["varOpex"], "其中电费": d["power_cost"],
            "其中原料(地沟油/CH₄)": d.get("feed_cost", 0), "其中CO₂": d["co2_cost"], "fixOpex": d["fixOpex"],
            "ARC": d["ARC"], "ACC": d["ACC"], "FCI工艺_MUSD": d["FCI_jet_MUSD"], "FCI电解_MUSD": d["FCI_el_MUSD"],
            "H2_kg/h": d["H2_kgh"]} for r, d in res.items()]
    ra.write_csv(OUT / "cn_7routes.csv", tab)
    ra.write_csv(OUT / "cn_sensitivity.csv", rows)
    wb = load_workbook(OUT / "SAF_TEA_results.xlsx")
    for name in ["中国7路线", "中国敏感性", "中国交叉点", "HEFA假设", "DRM-FT假设"]:
        if name in wb.sheetnames:
            del wb[name]
    wb.save(OUT / "SAF_TEA_results.xlsx")
    from openpyxl import load_workbook as lw
    wb = lw(OUT / "SAF_TEA_results.xlsx")
    ra.sheet(wb, "中国7路线", tab, "中国基准情景，USD/kg SAF；HEFA 与 DRM-FT 仅对中国计算")
    ra.sheet(wb, "中国敏感性", rows, "单因素，其余参数取基准；不适用的参数（如 HEFA 的 CO₂ 价）摆幅为 0")
    ra.sheet(wb, "中国交叉点", [{"指标": k, "数值": v} for k, v in cross.items()])
    ra.sheet(wb, "HEFA假设", [{"参数": k, "取值": v} for k, v in HEFA_ASSUMPTIONS.items()])
    ra.sheet(wb, "DRM-FT假设", [{"参数": k, "取值": v} for k, v in DRM_ASSUMPTIONS.items()]
             + [{"参数": "物料平衡: " + k, "取值": v} for k, v in drm_balance().items()])
    wb.save(OUT / "SAF_TEA_results.xlsx")

    for t in tab:
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in t.items()})
    print("cross", cross)
    print("DRM bal", drm_balance())
    for r in ROUTES_CN:
        top = sorted([z for z in rows if z["route"] == r], key=lambda z: -z["swing"])[:3]
        print(r, [(z["param"], round(z["NPC_low"], 2), round(z["NPC_high"], 2)) for z in top])


if __name__ == "__main__":
    main()
