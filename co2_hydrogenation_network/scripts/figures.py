"""Figures: species and reaction counts vs carbon number, and the reduction ladder."""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
FIG = os.path.join(HERE, "..", "figures")

S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"   # categorical slots 1-3
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "figure.facecolor": SURF,
    "axes.facecolor": SURF, "savefig.facecolor": SURF, "legend.frameon": False,
})


def fig_species(est):
    by = est["by_carbon"]
    cs = list(range(1, 13))
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for key, col, lab in (("mol", S1, "closed-shell molecules"), ("rad", S2, "monoradicals"),
                          ("mol_p", S3, "molecules, plausible tier")):
        y = [sum(by[f"{c},{o}"][key] for o in range(3)) for c in cs]
        ex = [all(by[f"{c},{o}"]["exact"] for o in range(3)) for c in cs]
        last = max(i for i, e in enumerate(ex) if e)
        ax.plot(cs[:last + 1], y[:last + 1], "-o", color=col, lw=2, ms=5, label=lab)
        ax.plot(cs[last:], y[last:], "--o", color=col, lw=2, ms=5, mfc=SURF)
        ax.annotate(f"{y[-1]:.1e}", (cs[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=8, color=INK2)
    ax.set_yscale("log")
    ax.set_xlabel("carbon number (O = 0-2 summed)")
    ax.set_ylabel("number of species (constitutional)")
    ax.set_xticks(cs)
    ax.set_xlim(0.5, 13.6)
    ax.legend(loc="upper left")
    ax.set_title("C/H/O species per carbon number\nfilled = exact enumeration, open/dashed = extrapolated",
                 fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_species_vs_carbon.png"), dpi=180)


def fig_reactions(est):
    bc = est["by_cmax"]
    cs = list(range(1, 13))
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for key, col, lab in (("quadratic_total", S2, "H/atom-transfer families (∝ N_rad²)"),
                          ("linear_total", S1, "all other families (∝ N_species)"),
                          ("species", S3, "species (molecules + radicals)")):
        y = [bc[str(c)]["full"][key] for c in cs]
        ax.plot(cs[:8], y[:8], "-o", color=col, lw=2, ms=5, label=lab)
        ax.plot(cs[7:], y[7:], "--o", color=col, lw=2, ms=5, mfc=SURF)
        ax.annotate(f"{y[-1]:.1e}", (cs[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=8, color=INK2)
    ax.set_yscale("log")
    ax.set_xlabel("carbon cap C_max of the network (O ≤ 2)")
    ax.set_ylabel("count")
    ax.set_xticks(cs)
    ax.set_xlim(0.5, 13.6)
    ax.legend(loc="upper left")
    ax.set_title("Full CO₂-hydrogenation network size vs carbon cap\n"
                 "open/dashed = extrapolated (some C·O classes not enumerated)",
                 fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_reactions_vs_cmax.png"), dpi=180)


def fig_ladder(est, red):
    full = est["central"]
    rows = [("Full space C≤12, O≤2", full["full"]["reactions_total"]),
            ("Plausible tier", full["plausible"]["reactions_total"])]
    names = {"A_methanation": "A methanation", "B_methanol": "B methanol", "C_CO2-FTS": "C CO₂-FTS",
             "D_CO2-to-aromatics": "D aromatics"}
    for tag in ("D_CO2-to-aromatics_O", "D_CO2-to-aromatics", "C_CO2-FTS_O", "C_CO2-FTS",
                "B_methanol_O", "B_methanol", "A_methanation_O", "A_methanation"):
        if tag not in red:
            continue
        base = tag[:-2] if tag.endswith("_O") else tag
        lab = names[base] + (" (+O, small-radical pool)" if tag.endswith("_O") else " (skeleton, pool)")
        rows.append((lab, red[tag]["reactions_pool"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    y = list(range(len(rows)))[::-1]
    vals = [r[1] for r in rows]
    cols = [S2 if i < 2 else S1 for i in range(len(rows))]
    ax.barh(y, vals, color=cols, height=0.62)
    for yi, v in zip(y, vals):
        ax.text(v * 1.6, yi, f"{v:.1e}", va="center", fontsize=8, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(10, max(vals) * 1e3)
    ax.set_xlabel("number of reversible elementary reactions (log)")
    ax.grid(axis="y", visible=False)
    ax.set_title("How much a known product slate shrinks the network", fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig3_reduction_ladder.png"), dpi=180)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    est = json.load(open(os.path.join(DATA, "estimate_C12.json")))
    fig_species(est)
    fig_reactions(est)
    p = os.path.join(DATA, "reduced_networks.json")
    if os.path.exists(p):
        fig_ladder(est, json.load(open(p)))
