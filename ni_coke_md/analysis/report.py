#!/usr/bin/env python
"""Campaign-level analysis: analyse every case, make comparison figures, write REPORT.md.

    python analysis/report.py [task ...] [--jobs 4] [--gif] [--out analysis_out]

Works on partially finished campaigns too (uses whatever frames exist), so you can run it
while jobs are still going.
"""
import argparse
import json
import os
import sys
import warnings
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from analyze_case import analyze  # noqa: E402

TASKS = ["task1_flux", "task2_anneal", "task3_codep"]
warnings.filterwarnings("ignore", category=RuntimeWarning)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK, INK2, GRID, BG = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": BG, "axes.facecolor": BG})


# ------------------------------------------------------------------ loading
def case_dirs(task):
    d = os.path.join(ROOT, "runs", task)
    if not os.path.isdir(d):
        return []
    return sorted(os.path.join(d, c) for c in os.listdir(d)
                  if os.path.exists(os.path.join(d, c, "config.json")))


def _analyze(c):
    try:
        return c, analyze(c)
    except Exception as e:  # keep going if one trajectory is broken
        return c, f"ERROR {e}"


def load(task, jobs):
    dirs = case_dirs(task)
    with Pool(jobs) as p:
        res = p.map(_analyze, dirs)
    out = []
    for c, m in res:
        if not m or str(m).startswith("ERROR"):
            if m:
                print("  !", c, m)
            continue
        d = np.genfromtxt(m, delimiter=",", names=True)
        if d.size < 2:
            continue
        cfg = json.load(open(os.path.join(c, "config.json")))
        st = os.path.join(c, "state.json")
        done = os.path.exists(st) and json.load(open(st)).get("phase") == "done"
        out.append(dict(dir=c, cfg=cfg, m=d, done=done))
    return out


# ------------------------------------------------------------------ reductions
def by_deposited(m, key, nmax):
    """metric value at the end of each deposition interval (last frame with C_ins == n)."""
    x = np.arange(1, nmax + 1)
    y = np.full(len(x), np.nan)
    for i, n in enumerate(x):
        sel = np.where(m["C_ins"] == n)[0]
        if len(sel):
            y[i] = m[key][sel[-1]]
    return x, y


def metric(m, key):
    if key == "C_bound":
        return m["nC"] - m["C_detached"]
    if key == "f6":
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(m["rings"] > 0, m["r6"] / np.maximum(m["rings"], 1), np.nan)
    if key == "sp2_frac":
        return m["sp2"] / np.maximum(m["nC"], 1)
    return m[key]


def with_derived(m):
    names = list(m.dtype.names) + ["C_bound", "f6", "sp2_frac"]
    out = {n: m[n] for n in m.dtype.names}
    for k in ("C_bound", "f6", "sp2_frac"):
        out[k] = metric(m, k)
    return out


def band(ys):
    ys = np.array(ys, float)
    with np.errstate(all="ignore"):
        return np.nanmean(ys, 0), np.nanmin(ys, 0), np.nanmax(ys, 0)


# ------------------------------------------------------------------ plotting
def panel(ax, groups, title, ylabel, xlabel, integer=True):
    for (lab, col, x, mu, lo, hi) in groups:
        ok = ~np.isnan(mu)
        if not ok.any():
            continue
        ax.fill_between(x[ok], lo[ok], hi[ok], color=col, alpha=0.15, lw=0)
        ax.plot(x[ok], mu[ok], color=col, lw=2, label=lab)
        ax.annotate(lab, (x[ok][-1], mu[ok][-1]), xytext=(4, 0), textcoords="offset points",
                    color=INK2, fontsize=8, va="center")
    ax.set_title(title, loc="left", color=INK, fontsize=10)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    if integer:
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_ylim(bottom=0, top=max(ax.get_ylim()[1], 1))
    else:
        ax.set_ylim(bottom=0, top=max(ax.get_ylim()[1], 0.05))
    xmax = max((g[2][~np.isnan(g[3])][-1] for g in groups if (~np.isnan(g[3])).any()), default=1)
    ax.set_xlim(0, xmax * 1.16)
    if len(groups) > 1:
        ax.legend(frameon=False, loc="upper left", labelcolor=INK2, fontsize=8)


def final_table(groups, keys, labels):
    """groups: list of (label, [metric dicts]); value = last frame of each case; mean ± sd."""
    head = "| group | n | " + " | ".join(labels) + " |\n|---|---|" + "---|" * len(labels) + "\n"
    rows = ""
    stats = {}
    for lab, ms in groups:
        vals = {k: np.array([m[k][-1] for m in ms], float) for k in keys}
        stats[lab] = vals
        cells = [f"{np.nanmean(v):.2f} ± {np.nanstd(v, ddof=1) if len(v) > 1 else 0:.2f}" for v in vals.values()]
        rows += f"| {lab} | {len(ms)} | " + " | ".join(cells) + " |\n"
    return head + rows, stats


# ------------------------------------------------------------------ tasks
KEYS_C = [("C_interior", "C inside Ni"), ("C_surface", "C on surface"), ("C6p", "C in ≥C6 clusters"),
          ("largest", "largest C cluster"), ("rings", "C rings (all sizes)"), ("sp2", "3-fold C (sp2-like)")]


def task1(cases, out, md):
    if not cases:
        return
    groups = {}
    for c in cases:
        groups.setdefault(c["cfg"]["c_interval_ps"], []).append(c)
    ivs = sorted(groups)
    nC = max(c["cfg"]["deposit"]["C"] for c in cases)
    fig, axs = plt.subplots(2, 3, figsize=(14, 8), dpi=140)
    for ax, (key, lab) in zip(axs.flat, KEYS_C):
        gs = []
        for i, iv in enumerate(ivs):
            ys = []
            for c in groups[iv]:
                x, y = by_deposited(with_derived(c["m"]), key, nC)
                ys.append(y)
            gs.append((f"{iv:g} ps/C", SERIES[i], x, *band(ys)))
        panel(ax, gs, lab, "number", "C atoms deposited")
    fig.suptitle("Task 1 · flux convergence — mean of seeds (line) and min–max (band), "
                 "values at the end of each insertion interval", x=0.01, ha="left", color=INK, fontsize=11)
    fig.tight_layout()
    f = os.path.join(out, "task1_flux.png")
    fig.savefig(f)
    plt.close(fig)
    md.append("## Task 1 — flux convergence (C only)\n\n![task1](task1_flux.png)\n\n")
    tab, stats = final_table([(f"{iv:g} ps/C", [with_derived(c["m"]) for c in groups[iv]]) for iv in ivs],
                             [k for k, _ in KEYS_C], [l for _, l in KEYS_C])
    md.append("Final frame of each case (after the post-deposition relaxation):\n\n" + tab + "\n")
    if len(ivs) >= 2:
        a, b = f"{ivs[-2]:g} ps/C", f"{ivs[-1]:g} ps/C"
        md.append(f"**Convergence check** — difference between the two slowest fluxes ({a} vs {b}) "
                  "in units of the seed-to-seed standard deviation (|Δ|/σ < 1 ≈ converged):\n\n")
        for k, lab in KEYS_C:
            va, vb = stats[a][k], stats[b][k]
            sd = np.sqrt((np.nanvar(va, ddof=1) + np.nanvar(vb, ddof=1)) / 2) if min(len(va), len(vb)) > 1 else np.nan
            d = np.nanmean(vb) - np.nanmean(va)
            ratio = abs(d) / sd if sd and sd > 0 else np.nan
            flag = "✓" if ratio < 1 else ("?" if np.isnan(ratio) else "✗")
            md.append(f"- {lab}: Δ = {d:+.2f}, |Δ|/σ = {ratio:.2f} {flag}\n")
        md.append("\n")
    for iv in ivs:
        n_done = sum(c["done"] for c in groups[iv])
        md.append(f"- {iv:g} ps/C: {n_done}/{len(groups[iv])} cases finished\n")
    md.append("\n")


def task2(cases, out, md):
    if not cases:
        return
    groups = {}
    for c in cases:
        groups.setdefault(c["cfg"]["T"], []).append(c)
    fig, axs = plt.subplots(2, 3, figsize=(14, 8), dpi=140)
    specs = [("r5", "5-rings"), ("r6", "6-rings"), ("r7", "7-rings"), ("f6", "fraction of rings that are 6-rings"),
             ("sp2_frac", "fraction of C that is 3-fold (sp2-like)"), ("largest", "largest C cluster")]
    for ax, (key, lab) in zip(axs.flat, specs):
        gs = []
        for i, T in enumerate(sorted(groups)):
            ys, tmax = [], 0
            for c in groups[T]:
                m = with_derived(c["m"])
                tmax = max(tmax, m["t_ps"][-1] - m["t_ps"][0])
            grid = np.linspace(0, tmax, 200)
            for c in groups[T]:
                m = with_derived(c["m"])
                t = m["t_ps"] - m["t_ps"][0]
                y = np.interp(grid, t, m[key], right=np.nan)
                ys.append(y)
            gs.append((f"{T:g} K", SERIES[i], grid, *band(ys)))
        panel(ax, gs, lab, "fraction" if "frac" in key or key == "f6" else "number",
              "anneal time (ps)", integer=not ("frac" in key or key == "f6"))
    fig.suptitle("Task 2 · annealing of the deposited carbon — mean of seeds and min–max band",
                 x=0.01, ha="left", color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "task2_anneal.png"))
    plt.close(fig)
    md.append("## Task 2 — annealing\n\n![task2](task2_anneal.png)\n\n")
    keys = [("r5", "5-rings"), ("r6", "6-rings"), ("r7", "7-rings"), ("f6", "6-ring fraction"),
            ("sp2_frac", "sp2 fraction"), ("largest", "largest cluster"), ("C_interior", "C inside Ni")]
    first = [(f"{T:g} K start", [{k: with_derived(c["m"])[k][:1] for k, _ in keys} for c in groups[T]])
             for T in sorted(groups)]
    last = [(f"{T:g} K end", [with_derived(c["m"]) for c in groups[T]]) for T in sorted(groups)]
    tab, _ = final_table(first + last, [k for k, _ in keys], [l for _, l in keys])
    md.append("Start vs end of the anneal:\n\n" + tab + "\n")


def task3(cases, out, md):
    if not cases:
        return
    groups = {}
    for c in cases:
        groups.setdefault(c["cfg"].get("O_per_C", 0), []).append(c)
    rs = sorted(groups)
    nC = max(c["cfg"]["deposit"]["C"] for c in cases)
    specs = [("C_bound", "C retained on the particle"), ("C_rem", "C removed (desorbed fragments)"),
             ("C_rem_CO", "C removed as CO"), ("C6p", "C in ≥C6 clusters"),
             ("rings", "C rings"), ("O_surface", "O on the surface")]
    fig, axs = plt.subplots(2, 3, figsize=(14, 8), dpi=140)
    for ax, (key, lab) in zip(axs.flat, specs):
        gs = []
        for i, r in enumerate(rs):
            ys = []
            for c in groups[r]:
                x, y = by_deposited(with_derived(c["m"]), key, nC)
                ys.append(y)
            gs.append((f"O/C = {r:g}", SERIES[i], x, *band(ys)))
        panel(ax, gs, lab, "number", "C atoms deposited")
    fig.suptitle("Task 3 · C + O co-deposition (O = oxygen from CO2 dissociation) — mean of seeds and min–max band",
                 x=0.01, ha="left", color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "task3_codep.png"))
    plt.close(fig)
    md.append("## Task 3 — C + O co-deposition\n\n![task3](task3_codep.png)\n\n")
    keys = [("C_bound", "C retained"), ("C_rem", "C removed"), ("C_rem_CO", "as CO"), ("C_rem_CO2", "as CO2"),
            ("C6p", "C in ≥C6"), ("rings", "rings"), ("O_interior", "O inside Ni"), ("O_surface", "O on surface"),
            ("O_rem", "O removed")]
    tab, _ = final_table([(f"O/C = {r:g}", [with_derived(c["m"]) for c in groups[r]]) for r in rs],
                         [k for k, _ in keys], [l for _, l in keys])
    md.append("Final frame of each case:\n\n" + tab + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tasks", nargs="*", default=[])
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--gif", action="store_true", help="also render a GIF for seed 1 of every group")
    ap.add_argument("--out", default=os.path.join(ROOT, "analysis_out"))
    a = ap.parse_args()
    tasks = a.tasks or TASKS
    os.makedirs(a.out, exist_ok=True)
    md = ["# Ni carbon-deposition campaign — report\n\n",
          "Generated by `analysis/report.py`. Lines = mean over seeds, bands = min–max over seeds.\n\n"]
    fns = dict(task1_flux=task1, task2_anneal=task2, task3_codep=task3)
    for t in tasks:
        print("analysing", t, flush=True)
        cs = load(t, a.jobs)
        print(f"  {len(cs)} case(s) with data")
        fns[t](cs, a.out, md)
        if a.gif:
            seen = set()
            for c in cs:
                key = (c["cfg"].get("c_interval_ps"), c["cfg"].get("O_per_C"), c["cfg"].get("T"))
                if key in seen:
                    continue
                seen.add(key)
                name = f"{t}_{os.path.basename(c['dir'])}.gif"
                os.system(f'"{sys.executable}" "{os.path.join(HERE, "render_gif.py")}" "{c["dir"]}" '
                          f'--out "{os.path.join(a.out, name)}"')
                md.append(f"![{name}]({name})\n\n")
    with open(os.path.join(a.out, "REPORT.md"), "w") as f:
        f.writelines(md)
    print("wrote", os.path.join(a.out, "REPORT.md"))


if __name__ == "__main__":
    main()
