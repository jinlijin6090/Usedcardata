"""Entry point:  python -m ecspace.main --config config.ini [--set section.key=value ...] [--stages ...]

Stages (default: all, in order):  estimate  enumerate  derive  solvent  stats  viewer
Every stage is resumable: finished layers are skipped when the job is restarted.
The closed-shell / derived space does not depend on the solvent; `solvent` (only when [solvent] solvents is set)
adds the bimolecular steps with the solvent molecules into out/solvent/<TAG>/, so one base run serves every
solvent system:  --stages solvent,stats,viewer  with another preset and the same run.out.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import random
import shutil
import sys
import time
from pathlib import Path

from . import chem, config, pipeline
from .pipeline import log

TYPE_CN = {
    "closed_shell": "闭壳层中性分子", "radical": "自由基 R•", "anion": "负离子 R⁻", "cation": "碳正离子 R⁺", "protonated": "质子化 MH⁺",
    "radical_anion": "自由基阴离子 M•⁻", "radical_cation": "自由基阳离子 M•⁺", "distonic_radical_anion": "远端自由基阴离子 •R⁻",
    "distonic_radical_cation": "远端自由基阳离子 •R⁺", "dianion": "二价阴离子", "biradical": "双自由基", "zwitterion": "两性离子",
}
FAM_CN = {
    "h_loss": "失 H•（R–H ⇌ R• + H•）", "deprotonation": "去质子", "hydride_loss": "失负氢 H⁻", "protonation": "O 质子化",
    "electron_capture": "得电子（M + e⁻ ⇌ M•⁻）", "electron_loss": "失电子（M ⇌ M•⁺ + e⁻）", "site_electron_transfer": "位点电子转移（R• ⇌ R⁻ / R⁺）",
    "li_binding": "Li⁺ 配位", "homolysis": "均裂 / 自由基复合", "heterolysis": "异裂 / 离子复合", "beta_scission": "自由基 β-断裂 / 加成",
    "radical_anion_fragmentation": "自由基阴离子断键（含开环）", "radical_cation_fragmentation": "自由基阳离子断键", "pair_steps": "双位点物种的逐步形成",
    "solvent_h_abstraction": "从溶剂分子抽 H（X• + S ⇌ XH + S•）", "solvent_deprotonation": "夺溶剂质子（X⁻ + S ⇌ XH + S⁻）",
    "solvent_acyl_transfer": "亲核酰基取代：酯交换 / 开环 / C-酰化（X⁻ + S）", "solvent_alkyl_transfer": "SN2 烷基转移：甲基化 / 醚断裂（X⁻ + S）",
}


def type_cn(t):
    if t.startswith("li"):
        k, base = t[2], t.split("_", 1)[1]
        return f"{TYPE_CN.get(base, base)} · {k} Li⁺"
    return TYPE_CN.get(t, t)


# ------------------------------------------------------------------ estimate
def estimate(cfg):
    """Run the first layers for real in a scratch directory, then extrapolate counts, disk and CPU time."""
    tmp = cfg.out / "_estimate"
    if tmp.exists():
        shutil.rmtree(tmp)
    small = config.Config(**{**cfg.__dict__})
    small.out = tmp
    small.max_heavy = min(cfg.max_heavy, 7)
    small.derive_max_heavy = 0
    small.count_sample = 1.0
    log("estimate: enumerating the first layers exactly (scratch run)")
    pipeline.enumerate_closed(small)
    counts = [json.loads((tmp / "closed" / f"L{n}" / "DONE.json").read_text())["count"] for n in range(1, small.max_heavy + 1)]
    # per-parent derivation cost and yield from a sample of the last exact layer
    last = small.max_heavy
    pars = [s for f in pipeline.layer_files(tmp, last) for s in pipeline.read_lines(f)]
    random.seed(1)
    smp = random.sample(pars, min(300, len(pars)))
    D = chem.Deriver(cfg)
    t0 = time.time()
    nsp = nrx = 0
    for p in smp:
        sp, rx = D.derive(p)
        nsp += len(sp)
        nrx += len(rx)
    dt = (time.time() - t0) / len(smp)
    sp_per, rx_per = nsp / len(smp), nrx / len(smp)
    # growth of the closed-shell layers: the measured ratio keeps rising by ~0.8 per layer for C/H/O
    ratios = [counts[i] / counts[i - 1] for i in range(2, len(counts))]
    r, step = ratios[-1], ratios[-1] - ratios[-2]
    est = counts[:]
    for n in range(len(counts) + 1, cfg.max_heavy + 1):
        r = r + max(0.6, min(1.0, step))
        est.append(int(est[-1] * r * (0.5 if n > max(cfg.max_C, cfg.max_O) else 1.0)))
    # calibrated on test runs: per-parent yields/costs grow linearly with size; compressed bytes per row
    B_CLOSED, B_SP, B_RX, DUP = 4.5, 7.1, 8.4, 4.8
    rows = []
    tot_sp = tot_rx = cpu = disk = 0
    tmp_peak = 0
    for n in range(1, cfg.max_heavy + 1):
        c = est[n - 1]
        d_n = max(0, n - last)
        spf, rxf = sp_per + 19 * d_n, rx_per + 45 * d_n
        s, x = c * spf, c * rxf
        w = n <= cfg.derive_max_heavy
        gen_cpu = (est[n - 2] if n > 1 else 0) * 0.0115 * (1.45 ** (n - 8))
        der_cpu = c * (dt + 0.010 * d_n) * (1.0 if w else cfg.count_sample)
        cpu += gen_cpu + der_cpu
        dsk = c * B_CLOSED + ((s * B_SP + x * B_RX) if w else 0)
        disk += dsk
        tmp_peak = max(tmp_peak, c * DUP * 6)
        tot_sp += s
        tot_rx += x
        rows.append((n, c, s, x, w, dsk, (gen_cpu + der_cpu) / 3600))
    W = max(1, cfg.workers)
    lines = ["", "规模估算（重原子 ≤ %d 为本机精确枚举，其余按实测增长倍数外推；派生物种、反应、CPU 和磁盘按测试运行标定，误差约 ±50%%）" % last, "",
             "| 重原子 | 闭壳层分子 | 全部物种（含自由基 / 离子 / Li） | 基元反应 | 逐条写出 | 磁盘（压缩后） | CPU 核时 |",
             "|---|---|---|---|---|---|---|"]
    for n, c, s, x, w, d, h in rows:
        lines.append(f"| {n} | {c:,.0f}{'' if n <= last else ' *'} | {s:,.0f} | {x:,.0f} | {'是' if w else '抽样计数'} | {d / 1e9:,.2f} GB | {h:,.1f} |")
    lines += ["", f"合计：物种约 {tot_sp:.2e}，基元反应约 {tot_rx:.2e}；写出数据约 {disk / 1e9:,.1f} GB，"
              f"枚举时临时文件峰值约 {tmp_peak / 1e9:,.1f} GB；CPU 约 {cpu / 3600:,.0f} 核时，用 {W} 核约 {cpu / 3600 / W:,.1f} 小时。",
              "内存：每个进程约 0.3–1 GB，去重时每个分片约 0.5 GB；建议按每核 2 GB 申请。（* 为外推值）", ""]
    text = "\n".join(lines)
    print(text, flush=True)
    (cfg.out / "ESTIMATE.md").write_text(text, encoding="utf-8")
    shutil.rmtree(tmp)


# ------------------------------------------------------------------ stats
def stats(cfg):
    out = cfg.out
    layers = []
    for n in range(1, cfg.max_heavy + 1):
        cj = out / "closed" / f"L{n}" / "DONE.json"
        dj = out / "derive" / f"L{n}.json"
        if not cj.exists():
            break
        c = json.loads(cj.read_text())
        d = json.loads(dj.read_text()) if dj.exists() else {}
        layers.append({"n": n, "closed_shell": c["count"], "exact": d.get("exact", True),
                       "species": d.get("species", {}), "reactions": d.get("reactions", {})})
    types = sorted({t for L in layers for t in L["species"] if not t.startswith("_")}, key=lambda t: (t.startswith("li"), t))
    fams = sorted({f for L in layers for f in L["reactions"]})
    tot_sp = sum(sum(v for k, v in L["species"].items() if not k.startswith("_")) for L in layers)
    tot_rx = sum(sum(L["reactions"].values()) for L in layers)
    res = {"config": cfg.source, "max_heavy": cfg.max_heavy, "derive_max_heavy": cfg.derive_max_heavy,
           "closed_shell_total": sum(L["closed_shell"] for L in layers), "species_total": tot_sp, "reactions_total": tot_rx,
           "layers": layers}
    from . import solvent as solv
    tag = solv.tag_of(cfg)
    sol = []
    if tag:
        for n in range(1, cfg.max_heavy + 1):
            sj = out / "solvent" / tag / f"L{n}.json"
            if sj.exists():
                sol.append(json.loads(sj.read_text()))
        res["solvent"] = {"tag": tag, "solvents": [k for _, k in solv.resolve(cfg.solvents)], "layers": sol,
                          "reactions_total": sum(v for L in sol for k, v in L["reactions"].items() if not k.startswith("_"))}
    (out / "stats.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    md = [f"# EC 化学空间枚举结果（重原子 ≤ {cfg.max_heavy}）", "",
          f"- 闭壳层中性分子：**{res['closed_shell_total']:,}**（精确）",
          f"- 全部物种（含自由基、离子、Li⁺ 配位物）：**{tot_sp:,}**",
          f"- 基元反应：**{tot_rx:,}**",
          f"- 重原子 ≤ {cfg.derive_max_heavy} 的物种和反应已逐条写入 species/ 与 reactions/；更大的层为抽样估算（标 ≈）", "",
          "## 按重原子数", "", "| 重原子 | 闭壳层分子 | 全部物种 | 基元反应 |", "|---|---|---|---|"]
    for L in layers:
        a = "" if L["exact"] else "≈ "
        md.append(f"| {L['n']} | {L['closed_shell']:,} | {a}{sum(v for k, v in L['species'].items() if not k.startswith('_')):,} | {a}{sum(L['reactions'].values()):,} |")
    md += ["", "## 按物种类别", "", "| 类别 | 数量 |", "|---|---|"]
    for t in types:
        md.append(f"| {type_cn(t)} | {sum(L['species'].get(t, 0) for L in layers):,} |")
    md += ["", "## 按反应类型", "", "| 反应类型 | 数量 |", "|---|---|"]
    for f in fams:
        md.append(f"| {FAM_CN.get(f, f)} | {sum(L['reactions'].get(f, 0) for L in layers):,} |")
    if tag and sol:
        sfam = sorted({f for L in sol for f in L["reactions"] if not f.startswith("_")})
        names = "、".join(f"{lab}（{k}）" for lab, k in solv.resolve(cfg.solvents))
        md += ["", f"## 与溶剂分子的双分子反应：{tag}", "", f"溶剂：{names}。每个自由基 / 负离子（定域、不含 Li）与每种溶剂分子之间的反应；"
               f"结果在 solvent/{tag}/。", "", "| 重原子（反应物 X） | 反应数 |", "|---|---|"]
        for L in sol:
            md.append(f"| {L['n']} | {'' if L['exact'] else '≈ '}{sum(v for k, v in L['reactions'].items() if not k.startswith('_')):,} |")
        md += ["", "| 反应类型 | 数量 |", "|---|---|"]
        for f in sfam:
            md.append(f"| {FAM_CN.get(f, f)} | {sum(L['reactions'].get(f, 0) for L in sol):,} |")
        md.append(f"| 合计 | {res['solvent']['reactions_total']:,} |")
    name = f"SUMMARY_{tag}.md" if tag else "SUMMARY.md"
    (out / name).write_text("\n".join(md) + "\n", encoding="utf-8")
    if tag:
        (out / f"stats_{tag}.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print("\n".join(md), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description="EC chemical-space brute-force enumeration")
    ap.add_argument("--config", default=None)
    ap.add_argument("--set", action="append", default=[], help="override, e.g. space.max_heavy=10")
    ap.add_argument("--stages", default="estimate,enumerate,derive,solvent,stats,viewer")
    a = ap.parse_args(argv)
    cfg = config.load(a.config, a.set)
    cfg.out.mkdir(parents=True, exist_ok=True)
    stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    log(f"chemical space  solvents={'+'.join(cfg.solvents) or '-'}  start={','.join(cfg.start)}  max_heavy={cfg.max_heavy}  derive_max_heavy={cfg.derive_max_heavy}  workers={cfg.workers}  out={cfg.out}")
    (cfg.out / "config_used.json").write_text(json.dumps({k: (str(v) if isinstance(v, Path) else v) for k, v in cfg.__dict__.items()},
                                                         indent=1, ensure_ascii=False))
    if "estimate" in stages or cfg.dry_run:
        estimate(cfg)
        if cfg.dry_run:
            log("dry_run = yes: stopping after the estimate")
            return
    if "enumerate" in stages:
        pipeline.enumerate_closed(cfg)
    if "derive" in stages:
        pipeline.write_headers(cfg)
        pipeline.derive_all(cfg)
    if "solvent" in stages:
        from . import solvent
        solvent.run(cfg)
    if "stats" in stages:
        stats(cfg)
    if "viewer" in stages:
        from . import viewer
        viewer.build(cfg)
    log("all stages finished")


if __name__ == "__main__":
    main()
