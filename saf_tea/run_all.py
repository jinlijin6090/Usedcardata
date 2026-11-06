"""Run: (1) replication of the 12 published scenarios, (2) 5 routes x 4 regions,
(3) sensitivities for China and EtJ-2. Writes results/*.csv, results/*.png and
results/SAF_TEA_results.xlsx.
"""
from __future__ import annotations

import copy
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

import extensions as x
import model as m

OUT = m.HERE / "results"
OUT.mkdir(exist_ok=True)
for f in font_manager.findSystemFonts():
    if "wqy-zenhei" in f:
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = "WenQuanYi Zen Hei"
plt.rcParams["axes.unicode_minus"] = False

ROUTES = ["FTS", "MtJ-1", "MtJ-2", "EtJ", "EtJ-2"]
REGS = ["EU", "US", "ME", "CN"]
REG_CN = {"EU": "欧洲", "US": "美国", "ME": "中东", "CN": "中国"}
COLORS = {"varOpex": "#F2C14E", "fixOpex": "#5B8DD6", "ARC": "#8CC084", "ACC": "#D9534F"}

EC = m.load_equipment()
VO = m.load_varopex_si()
Q = m.quantities_from_eu(VO)
ETJ_PMW = x.etj_process_mw(VO)


def scenario(route: str, reg_name: str, cn: dict | None = None, e2: dict | None = None,
             use_si_varopex: bool = True) -> dict:
    cn = cn or x.CN_ASSUMPTIONS
    e2 = e2 or x.ETJ2_ASSUMPTIONS
    reg = x.china_region(cn) if reg_name == "CN" else m.REGIONS[reg_name]
    list_reg = "US" if reg_name == "CN" else reg_name           # CN: US equipment list x location factor
    if route == "EtJ-2":
        ec = x.etj2_equipment(e2)[list_reg]
        h2 = x.etj2_balance(e2)["H2_kgh"]
        q = x.etj2_quantities(Q["EtJ"], e2, ETJ_PMW)
        vo = m.varopex_from_quantities(q, reg)
    else:
        ec = EC[route][list_reg]
        h2 = m.H2_KGH[route]
        if use_si_varopex and reg_name != "CN":
            vo = VO[(route, reg_name)]
        else:
            vo = m.varopex_from_quantities(Q[route], reg)
    ec *= reg.capex_factor
    safety_backup = m.ETJ_SAFETY
    if route == "EtJ-2":
        m.ETJ_SAFETY = e2["capex_safety"]
    res = m.npc(route, reg, ec, h2, vo["TOTAL"])
    m.ETJ_SAFETY = safety_backup
    res.update(route=route, region=reg_name, H2_kgh=h2, elec_usd_mwh=reg.elec,
               power_cost=vo.get("electricity", 0) / m.M_SAF, co2_cost=vo.get("co2", 0) / m.M_SAF)
    return res


# ------------------------------------------------------------------ 1. replication
def replication() -> list[dict]:
    fig = {(r["region"], r["route"]): r for r in csv.DictReader(open(m.HERE / "data/fig8_digitized.csv"))}
    rows = []
    for reg in ["EU", "US", "ME"]:
        for route in ROUTES[:4]:
            r = scenario(route, reg)
            f = fig[(reg, route)]
            rows.append({"region": reg, "route": route,
                         **{f"{k}_model": round(r[k], 3) for k in ["varOpex", "fixOpex", "ARC", "ACC", "NPC"]},
                         **{f"{k}_paper": float(f[k]) for k in ["varOpex", "fixOpex", "ARC", "ACC", "NPC"]},
                         "dev_%": round((r["NPC"] / float(f["NPC"]) - 1) * 100, 1)})
    return rows


# ------------------------------------------------------------------ 3. sensitivities
def tornado(base_fn, params: list[tuple[str, str, float, float]], store: dict):
    base = base_fn(store)["NPC"]
    out = []
    for key, label, lo, hi in params:
        vals = []
        for v in (lo, hi):
            s = copy.deepcopy(store)
            s[key] = v
            vals.append(base_fn(s)["NPC"])
        out.append({"param": label, "low_value": lo, "high_value": hi,
                    "NPC_low": round(vals[0], 3), "NPC_high": round(vals[1], 3), "base": round(base, 3)})
    return out


def cn_sens():
    params = [
        ("elec_base_cny_kwh", "绿电价格 元/kWh (0.15–0.40)", 0.15, 0.40),
        ("el_usd_per_kw", "电解槽设备价 USD/kW (250–947)", 250, 947),
        ("location_factor", "设备地区系数 (0.5–1.0)", 0.5, 1.0),
        ("co2_feed_cny_t", "CO₂原料 元/t (100–905)", 100, 905),
        ("ir", "贷款利率 (2.5%–5%)", 0.025, 0.05),
        ("co2_price_cny_t", "全国碳价 元/t (40–200)", 40, 200),
    ]
    return tornado(lambda s: scenario("MtJ-2", "CN", cn=s), params, x.CN_ASSUMPTIONS)


def etj2_sens(region="CN"):
    params = [
        ("carbon_yield", "乙醇碳选择性 (0.80–1.00)", 0.80, 1.00),
        ("single_pass_conv", "CO₂单程转化率 (0.10–0.35)", 0.10, 0.35),
        ("capex_safety", "低TRL投资裕量 (1.0–1.5)", 1.0, 1.5),
        ("steam_recovery", "反应热回收为蒸汽 (0–0.8)", 0.0, 0.8),
        ("p_syn_bar", "合成压力 bar (30–80)", 30, 80),
    ]
    return tornado(lambda s: scenario("EtJ-2", region, e2=s), params, x.ETJ2_ASSUMPTIONS)


def decomposition(route="MtJ-2"):
    """Walk from the US case to the China case one parameter group at a time."""
    from dataclasses import replace
    us, cn = m.REGIONS["US"], x.china_region()
    steps = [("美国基准", {}),
             ("电价 → 中国", {"elec": cn.elec}),
             ("电解槽单价 → 中国", {"el_usd_per_kw": cn.el_usd_per_kw}),
             ("设备地区系数 → 0.61", {"capex_factor": cn.capex_factor}),
             ("CO₂原料价 → 中国", {"co2": cn.co2}),
             ("利率 → 3.5%", {"ir": cn.ir}),
             ("人工、水、副产品等 → 中国", {f: getattr(cn, f) for f in ["labor_usd_h", "water", "mp_steam", "lp_steam",
                                                     "lights", "naphtha", "diesel", "wax", "cooling", "waste"]})]
    reg, rows, prev = us, [], None
    for label, ch in steps:
        reg = replace(reg, **ch)
        ec = (x.etj2_equipment()["US"] if route == "EtJ-2" else EC[route]["US"]) * reg.capex_factor
        q = x.etj2_quantities(Q["EtJ"], x.ETJ2_ASSUMPTIONS, ETJ_PMW) if route == "EtJ-2" else Q[route]
        h2 = x.etj2_balance()["H2_kgh"] if route == "EtJ-2" else m.H2_KGH[route]
        npc = m.npc(route, reg, ec, h2, m.varopex_from_quantities(q, reg)["TOTAL"])["NPC"]
        rows.append({"step": label, "NPC": round(npc, 3), "delta": round(npc - prev, 3) if prev else 0.0})
        prev = npc
    return rows


def etj2_breakeven(region="CN"):
    target = scenario("EtJ", region)["NPC"]
    lo, hi = 0.5, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2
        s = dict(x.ETJ2_ASSUMPTIONS, carbon_yield=mid)
        if scenario("EtJ-2", region, e2=s)["NPC"] > target:
            lo = mid
        else:
            hi = mid
    return hi, target


# ------------------------------------------------------------------ outputs
def plot_npc(results):
    fig, ax = plt.subplots(figsize=(15, 6.5), dpi=150)
    xs, labels = [], []
    pos = 0
    for reg in REGS:
        for route in ROUTES:
            r = next(rr for rr in results if rr["region"] == reg and rr["route"] == route)
            bottom = 0
            for k in ["varOpex", "fixOpex", "ARC", "ACC"]:
                ax.bar(pos, r[k], bottom=bottom, color=COLORS[k], edgecolor="#333", linewidth=0.6,
                       label=k if pos == 0 else None)
                bottom += r[k]
            ax.text(pos, bottom + 0.08, f"{bottom:.2f}", ha="center", fontsize=8.5)
            xs.append(pos); labels.append(route)
            pos += 1
        ax.text(pos - 3, -1.05, REG_CN[reg], ha="center", fontsize=12, fontweight="bold")
        pos += 0.8
    ax.axhline(1.036, ls="--", color="k", lw=1.2, label="化石 Jet A-1 (1.036)")
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("净生产成本 NPC (USD/kg SAF)")
    ax.set_ylim(0, 9.2)
    ax.legend(ncol=5, loc="upper left", fontsize=9)
    ax.set_title("5 条路线 × 4 个地区的 SAF 净生产成本（2025 年美元）")
    ax.grid(axis="y", ls=":", alpha=0.5)
    fig.tight_layout(); fig.subplots_adjust(bottom=0.13)
    fig.savefig(OUT / "npc_20_scenarios.png"); plt.close(fig)


def plot_parity(rep):
    fig, ax = plt.subplots(figsize=(5.5, 5.5), dpi=150)
    p = [r["NPC_paper"] for r in rep]; mm = [r["NPC_model"] for r in rep]
    ax.plot([2, 9], [2, 9], "k--", lw=1)
    for r in rep:
        ax.scatter(r["NPC_paper"], r["NPC_model"], s=40,
                   color={"EU": "#3B6EA8", "US": "#D9534F", "ME": "#4C9A2A"}[r["region"]])
        ax.annotate(f'{r["region"]}-{r["route"]}', (r["NPC_paper"], r["NPC_model"]), fontsize=7,
                    xytext=(4, -8), textcoords="offset points")
    ax.set_xlabel("论文 图 8 读数 (USD/kg)"); ax.set_ylabel("本模型 (USD/kg)")
    ax.set_title("复现精度：12 个场景")
    ax.grid(ls=":", alpha=0.5); fig.tight_layout()
    fig.savefig(OUT / "replication_parity.png"); plt.close(fig)


def plot_tornado(rows, title, fname):
    rows = sorted(rows, key=lambda r: abs(r["NPC_high"] - r["NPC_low"]))
    base = rows[0]["base"]
    fig, ax = plt.subplots(figsize=(9, 0.6 * len(rows) + 1.6), dpi=150)
    for i, r in enumerate(rows):
        for v, c in ((r["NPC_low"], "#3B6EA8"), (r["NPC_high"], "#E8A33D")):
            ax.barh(i, v - base, left=base, color=c, edgecolor="#333", lw=0.5)
            ax.text(v + (0.02 if v >= base else -0.02), i, f"{v:.2f}", va="center",
                    ha="left" if v >= base else "right", fontsize=8)
    ax.axvline(base, color="k", lw=1)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r["param"] for r in rows], fontsize=9)
    ax.set_xlabel(f"NPC (USD/kg)，基准 = {base:.2f}；蓝 = 区间下限取值，橙 = 上限取值")
    ax.set_title(title); ax.grid(axis="x", ls=":", alpha=0.5)
    fig.tight_layout(); fig.savefig(OUT / fname); plt.close(fig)


def write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def sheet(wb, title, rows, note=None):
    ws = wb.create_sheet(title)
    r0 = 1
    if note:
        ws.cell(1, 1, note).font = Font(italic=True, color="555555")
        r0 = 3
    keys = list(rows[0].keys())
    for j, k in enumerate(keys, 1):
        c = ws.cell(r0, j, k); c.font = Font(bold=True, color="FFFFFF")
        c.fill = PatternFill("solid", fgColor="14213D"); c.alignment = Alignment(wrap_text=True)
        ws.column_dimensions[get_column_letter(j)].width = max(12, min(40, len(str(k)) + 4))
    for i, row in enumerate(rows, r0 + 1):
        for j, k in enumerate(keys, 1):
            v = row[k]
            ws.cell(i, j, round(v, 4) if isinstance(v, float) else v)
    ws.freeze_panes = ws.cell(r0 + 1, 1)


def main():
    rep = replication()
    results = [scenario(rt, rg) for rg in REGS for rt in ROUTES]
    res_rows = [{"地区": REG_CN[r["region"]], "路线": r["route"], "NPC": r["NPC"], "varOpex": r["varOpex"],
                 "其中电费": r["power_cost"], "其中CO2原料": r["co2_cost"], "fixOpex": r["fixOpex"],
                 "ARC(电堆更换)": r["ARC"], "ACC(资本回收)": r["ACC"], "FCI工艺段_MUSD": r["FCI_jet_MUSD"],
                 "FCI电解_MUSD": r["FCI_el_MUSD"], "H2_kg/h": r["H2_kgh"], "电价_USD/MWh": r["elec_usd_mwh"]}
                for r in results]
    cns = cn_sens(); e2s = etj2_sens("CN"); e2s_eu = etj2_sens("EU")
    bal = x.etj2_balance()
    dec = decomposition("MtJ-2")
    be = {rg: etj2_breakeven(rg) for rg in REGS}
    write_csv(OUT / "decomposition_us_to_cn_mtj2.csv", dec)

    write_csv(OUT / "replication.csv", rep)
    write_csv(OUT / "npc_20_scenarios.csv", res_rows)
    write_csv(OUT / "sens_china_mtj2.csv", cns)
    write_csv(OUT / "sens_etj2_china.csv", e2s)
    plot_npc(results); plot_parity(rep)
    plot_tornado(cns, "中国 MtJ-2：参数敏感性", "tornado_china_mtj2.png")
    plot_tornado(e2s, "中国 EtJ-2：工艺假设敏感性", "tornado_etj2_china.png")

    regions = {k: m.REGIONS[k] for k in ["EU", "US", "ME"]} | {"CN": x.china_region()}
    reg_rows = [{"参数": f, **{REG_CN[k]: getattr(v, f) for k, v in regions.items()}}
                for f in ["ir", "labor_usd_h", "elec", "water", "co2", "mp_steam", "lp_steam", "lights",
                          "naphtha", "diesel", "wax", "cooling", "waste", "capex_factor", "el_usd_per_kw"]]
    wb = Workbook(); wb.remove(wb.active)
    sheet(wb, "20场景结果", res_rows, "单位：USD/kg SAF（2025年美元）；SAF 产量 13.7 t/h × 7884 h = 108 kt/a")
    sheet(wb, "复现对比", rep, "论文值为图 8 像素读数（±0.03 USD/kg）")
    sheet(wb, "地区假设(USD)", reg_rows, "EU/US/ME 取自论文表 1–3（EU 已按 0.8449 换成 USD）；CN 见 中国假设 表")
    sheet(wb, "中国假设", [{"参数": k, "取值": v} for k, v in x.CN_ASSUMPTIONS.items()],
          f"汇率 1 USD = {x.CNY_PER_USD} CNY（2025-09-19 中间价）")
    sheet(wb, "EtJ-2假设", [{"参数": k, "取值": v} for k, v in x.ETJ2_ASSUMPTIONS.items()]
          + [{"参数": "物料平衡: " + k, "取值": v} for k, v in bal.items()])
    sheet(wb, "美国→中国分解(MtJ-2)", dec, "从美国 MtJ-2 出发，逐项替换成中国参数（varOpex 按物理量×单价计算）")
    sheet(wb, "EtJ-2盈亏平衡", [{"地区": REG_CN[k], "EtJ_NPC": v[1], "EtJ-2追平EtJ所需乙醇碳选择性": v[0]}
                               for k, v in be.items()])
    sheet(wb, "敏感性_中国MtJ-2", cns)
    sheet(wb, "敏感性_中国EtJ-2", e2s)
    sheet(wb, "敏感性_欧洲EtJ-2", e2s_eu)
    wb.save(OUT / "SAF_TEA_results.xlsx")

    print("Replication |dev| max %.1f %%" % max(abs(r["dev_%"]) for r in rep))
    for r in res_rows:
        print(f'{r["地区"]} {r["路线"]:6s} NPC {r["NPC"]:.2f}  var {r["varOpex"]:.2f} (电 {r["其中电费"]:.2f}, CO2 {r["其中CO2原料"]:.2f}) '
              f'fix {r["fixOpex"]:.2f} ARC {r["ARC(电堆更换)"]:.2f} ACC {r["ACC(资本回收)"]:.2f}  FCIjet {r["FCI工艺段_MUSD"]:.0f} FCIel {r["FCI电解_MUSD"]:.0f}')
    print("decomp", dec); print("breakeven", be)
    print("CN sens", cns); print("EtJ2 CN sens", e2s); print("bal", bal)


if __name__ == "__main__":
    main()
