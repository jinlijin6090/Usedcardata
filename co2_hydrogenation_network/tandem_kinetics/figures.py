"""Figures for the OX-ZEO kinetic study."""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from model import Conditions, Model  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "font.size": 9.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "savefig.facecolor": SURF, "legend.frameon": False})
CLS = [("P", "paraffins", S[0]), ("O", "olefins", S[1]), ("N", "cyclic (naphthenic)", S[2]),
       ("A", "aromatics", S[3])]
LABEL = {"G1_GaZrOx_GaZSM5": "G1 GaZrOx/Ga-ZSM-5 (Guo, 2026)",
         "G2_ZnZrO_ZSM5": "G2 ZnO/ZrO2-ZSM-5 (Guo/Song, 2019)",
         "B1_In2O3_HZSM5": "B1 In2O3/HZSM-5 (benchmark, Gao 2017)"}


def fig_distributions(cal):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for ax, (name, r) in zip(axes, cal.items()):
        p = {k: tuple(v) for k, v in r["params"].items()}
        d = Model(p).run(Conditions(**r["cond"])).dist()
        ns = np.arange(1, 21)
        bottom = np.zeros(len(ns))
        for k, lab, col in CLS:
            y = np.array([d.get(k, {}).get(n, 0) * 100 for n in ns])
            ax.bar(ns, y, bottom=bottom, color=col, width=0.8, label=lab, edgecolor=SURF, linewidth=0.8)
            bottom += y
        ax.axvspan(7.5, 16.5, color=INK2, alpha=0.06, lw=0)
        ax.text(12, 0.97, "C8-C16 (jet)", ha="center", va="top", fontsize=8, color=INK2,
                transform=ax.get_xaxis_transform())
        ax.set_title(LABEL[name], fontsize=9, loc="left", color=INK)
        ax.set_xlabel("carbon number")
        ax.set_xticks([1, 4, 8, 12, 16, 20])
    axes[0].set_ylabel("carbon selectivity in hydrocarbons (%)")
    axes[-1].legend(loc="center right", fontsize=8)
    fig.suptitle("Calibrated model: predicted hydrocarbon distribution (only X, S_CO, C5+/aromatics, CH4 were fitted)",
                 fontsize=9.5, x=0.01, ha="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "f1_distributions.png"), dpi=180)


def fig_identifiability(ens):
    good = [e for e in ens if e["ok"]]
    keys = [("X_CO2", "CO2 conv.", True), ("S_CO", "CO sel.", True), ("C5+", "C5+", True),
            ("CH4", "CH4", True), ("aromatics", "aromatics", False), ("C2-C4o", "C2-C4 paraffins", False),
            ("C8-C16", "C8-C16", False), ("SAF_usable", "SAF-usable C8-C16", False)]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    rng = np.random.default_rng(0)
    for i, (k, lab, fitted) in enumerate(keys):
        v = np.array([e["pred"][k] for e in good]) * 100
        x = i + rng.uniform(-0.15, 0.15, len(v))
        ax.scatter(x, v, s=18, color=S[0] if fitted else S[1], alpha=0.8, edgecolor=SURF, linewidth=0.5)
        ax.text(i, v.max() + 3, f"{v.min():.0f}-{v.max():.0f}", ha="center", fontsize=8, color=INK2)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([k[1] for k in keys], rotation=20, ha="right")
    ax.set_ylabel("%")
    ax.scatter([], [], color=S[0], label="fitted to G1 data")
    ax.scatter([], [], color=S[1], label="not measured -> predicted")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title(f"{len(good)} parameter sets that all reproduce the G1 data: what is and is not pinned down",
                 fontsize=9.5, loc="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "f2_identifiability.png"), dpi=180)


def fig_sweeps(sw):
    vals, mem = sw["values"], sw["members"]
    xl = {"T_C": "T (°C)", "P_MPa": "P (MPa)", "GHSV": "GHSV (mL g⁻¹ h⁻¹)", "H2_CO2": "H2/CO2",
          "phi_ox": "oxide mass fraction"}
    rows = [("X_CO2", "CO2 conversion (%)"), ("S_CO", "CO selectivity (%)"),
            ("C8-C16", "C8-C16 in HC (%)")]
    fig, axes = plt.subplots(3, 5, figsize=(12, 6.6), sharey="row")
    for j, var in enumerate(vals):
        x = vals[var]
        for i, (k, lab) in enumerate(rows):
            ax = axes[i, j]
            Y = np.array([[r[k] * 100 if r else np.nan for r in m[var]] for m in mem])
            ax.fill_between(x, np.nanmin(Y[1:], 0), np.nanmax(Y[1:], 0), color=S[0], alpha=0.15, lw=0)
            ax.plot(x, Y[0], "-o", color=S[0], lw=2, ms=4)
            if var == "GHSV":
                ax.set_xscale("log")
            if i == 2:
                ax.set_xlabel(xl[var])
            if j == 0:
                ax.set_ylabel(lab)
    fig.suptitle("Operating sweeps around G1 (line = central fit, band = range over ensemble members)",
                 fontsize=9.5, x=0.01, ha="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "f3_sweeps.png"), dpi=180)


def fig_pore(caps, caps_e):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    ax = axes[0]
    for nca, col in ((10, S[0]), (12, S[1])):
        rows = [c for c in caps if c["n_capA"] == nca and c["pred"]]
        x = [c["n_cap"] for c in rows]
        ax.plot(x, [c["pred"]["C8-C16"] * 100 for c in rows], "-o", color=col, lw=2, ms=4,
                label=f"C8-C16, aromatic cap C{nca}")
        ax.plot(x, [c["pred"]["SAF_usable"] * 100 for c in rows], "--o", color=col, lw=2, ms=4,
                mfc=SURF, label=f"SAF-usable, aromatic cap C{nca}")
    ax.set_xlabel("chain-growth cap n_cap (pore size proxy: TON~9, MFI~11, BEA~13-14, FAU~16)")
    ax.set_ylabel("% of hydrocarbons")
    ax.legend(fontsize=7.5)
    ax.set_title("Pore-size proxy vs jet cut (G1 central fit)", fontsize=9.5, loc="left")
    ax = axes[1]
    rows = [c for c in caps if c["n_capA"] == 12 and c["pred"]]
    x = np.array([c["n_cap"] for c in rows])
    bottom = np.zeros(len(x))
    for k, lab, col in CLS:
        y = np.array([c["pred"]["jet_cut_comp"][k] * 100 for c in rows])
        ax.bar(x, y, bottom=bottom, color=col, label=lab, width=0.8, edgecolor=SURF)
        bottom += y
    ax.axhline(25, color=INK2, lw=1, ls=":")
    ax.text(x[0] - 0.4, 27, "25 % aromatics limit (D7566 blend)", fontsize=7.5, color=INK2)
    ax.set_xlabel("n_cap")
    ax.set_ylabel("composition of C8-C16 cut (%)")
    ax.legend(fontsize=7.5, loc="lower right")
    ax.set_title("What the jet cut is made of", fontsize=9.5, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "f4_pore_cap.png"), dpi=180)


def fig_drc(d):
    ks = list(d["drc"])
    outs = ["X_CO2", "S_CO", "C5+", "aromatics", "C8-C16", "Y_C8-C16"]
    M = np.array([[d["drc"][k][o] for o in outs] for k in ks])
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    lim = np.nanmax(np.abs(M))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(outs)))
    ax.set_xticklabels(outs, rotation=25, ha="right")
    ax.set_yticks(range(len(ks)))
    ax.set_yticklabels(ks)
    for i in range(len(ks)):
        for j in range(len(outs)):
            if abs(M[i, j]) >= 0.05:
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(M[i, j]) > 0.6 * lim else INK)
    ax.grid(False)
    fig.colorbar(im, ax=ax, label="d ln(output) / d ln(k)")
    ax.set_title("Sensitivity (degree of rate control), G1 central fit", fontsize=9.5, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "f5_drc.png"), dpi=180)


def fig_ceiling(c):
    fig, ax = plt.subplots(figsize=(6, 3.6))
    for p, row, col in zip(c["nP"], c["S_A_max"], S + [INK2]):
        ax.plot(c["nA"], np.array(row) * 100, "-", color=col, lw=2, label=f"mean paraffin C{p:g}")
    ax.axhline(70, color=INK2, ls=":", lw=1)
    ax.text(7.1, 72, "70 % aromatics reported for G2", fontsize=8, color=INK2)
    ax.set_xlabel("mean aromatic carbon number (MFI: ~C8-C10)")
    ax.set_ylabel("max aromatics, % C in HC")
    ax.legend(fontsize=7.5, ncol=2)
    ax.set_title("Hydride-transfer-only ceiling: S_A <= n_A / (n_A + 3 n_P)", fontsize=9.5, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "f6_ht_ceiling.png"), dpi=180)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    cal = json.load(open(os.path.join(HERE, "calibration.json")))
    fig_distributions(cal)
    fig_ceiling(json.load(open(os.path.join(RES, "ht_ceiling.json"))))
    if os.path.exists(os.path.join(RES, "ensemble_G1.json")):
        fig_identifiability(json.load(open(os.path.join(RES, "ensemble_G1.json"))))
    if os.path.exists(os.path.join(RES, "sweeps_G1.json")):
        fig_sweeps(json.load(open(os.path.join(RES, "sweeps_G1.json"))))
    if os.path.exists(os.path.join(RES, "pore_cap_G1.json")):
        fig_pore(json.load(open(os.path.join(RES, "pore_cap_G1.json"))),
                 json.load(open(os.path.join(RES, "pore_cap_ensemble.json"))))
    if os.path.exists(os.path.join(RES, "drc_G1.json")):
        fig_drc(json.load(open(os.path.join(RES, "drc_G1.json"))))
