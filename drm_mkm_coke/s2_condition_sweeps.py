"""Operating-condition and geometry sweeps of the coking model (flux variant)."""
from coke_tools import *
EA_SET = [1.5, 2.2]
def run(name, xs, pp_of, xlabel, logx=False):
    frames = []
    for Ea in EA_SET + [None]:
        if Ea is None:
            rc = load_model(coke=False); lab = "no coke (paper)"
        else:
            rc = load_model(CokeParams(Ea_cc=Ea)); lab = f"Ea_cc={Ea} eV"
        df = sweep(xs, pp_of, lambda v: rc); df["case"] = lab; frames.append(df)
    df = pd.concat(frames); df.to_csv(f"{RES}/s2_{name}.csv", index=False)
    fig, ax = plt.subplots(2, 3, figsize=(13, 7))
    ax = ax.ravel()
    for lab, g in df.groupby("case", sort=False):
        ax[0].plot(g.x, g.r_DRM, marker=".", label=lab)
        if "no coke" not in lab:
            ax[1].semilogy(g.x, g.r_coke.clip(1e-14), marker=".", label=lab)
            ax[2].semilogy(g.x, g.S_coke.clip(1e-14), marker=".", label=lab)
            ax[5].semilogy(g.x, g["th_C2*"].clip(1e-14), marker=".", label=lab)
        ax[3].semilogy(g.x, g["th_C*"].clip(1e-14), marker=".", label=lab)
        ax[4].semilogy(g.x, g["th_O*"].clip(1e-14), marker=".", label=lab)
    ax[0].set(ylabel="DRM rate r$_{H_2}$ [mol g$^{-1}$ s$^{-1}$]")
    ax[1].set(ylabel="coke rate [mol C g$^{-1}$ s$^{-1}$]")
    ax[2].set(ylabel="coke selectivity 2r$_{20}$/r$_{CH_4}$")
    ax[3].set(ylabel="θ(C*)"); ax[4].set(ylabel="θ(O*)"); ax[5].set(ylabel="θ(C$_2$*)")
    for a in ax:
        a.set_xlabel(xlabel)
        if logx: a.set_xscale("log")
    ax[0].legend(fontsize=8)
    fig.suptitle(f"Sweep: {name}", fontsize=10); fig.tight_layout()
    fig.savefig(f"{FIG}/s2_{name}.png"); plt.close(fig)
    return df

base_T = 973.15
d1 = run("pCH4", np.round(np.arange(0.1, 2.01, 0.1), 2), lambda v: pp_with(pCH4=v), "p$_{CH_4}$ [bar]  (p$_{CO_2}$=0.5)")
d2 = run("pCO2", np.concatenate([[0.01, 0.02, 0.05], np.round(np.arange(0.1, 2.01, 0.1), 2)]), lambda v: pp_with(pCO2=v), "p$_{CO_2}$ [bar]  (p$_{CH_4}$=0.5)", logx=True)
d3 = run("rm", np.arange(2, 10.01, 0.5), lambda v: pp_with(rm=v), "Ni radius r$_m$ [nm]")
d4 = run("rint", np.round(np.arange(0.38, 1.91, 0.19), 2), lambda v: pp_with(rsm=v), "interfacial radius r$_{int}$ [nm]")
d5 = run("pH2_cofeed", np.array([0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]), lambda v: pp_with(pH2=v), "p$_{H_2}$ cofeed [bar]")
# temperature: rate constants recomputed at each T (prefactors fixed, as in the workbook)
frames = []
for Ea in EA_SET + [None]:
    rows = []
    for T in np.arange(773.15, 1123.16, 25):
        rc = load_model(coke=False, T=T) if Ea is None else load_model(CokeParams(Ea_cc=Ea), T=T)
        o = solve(rc, pp_with(T=T)); o = {k: v for k, v in o.items() if not k.startswith("_")}; o["x"] = T; rows.append(o)
    df = pd.DataFrame(rows); df["case"] = "no coke (paper)" if Ea is None else f"Ea_cc={Ea} eV"; frames.append(df)
dT = pd.concat(frames); dT.to_csv(f"{RES}/s2_T.csv", index=False)
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
for lab, g in dT.groupby("case", sort=False):
    ax[0].semilogy(g.x, g.r_DRM, marker=".", label=lab)
    if "no coke" not in lab:
        ax[1].semilogy(g.x, g.r_coke.clip(1e-14), marker=".", label=lab); ax[2].semilogy(g.x, g.S_coke.clip(1e-14), marker=".", label=lab)
ax[0].set(xlabel="T [K]", ylabel="DRM rate"); ax[1].set(xlabel="T [K]", ylabel="coke rate [mol C g$^{-1}$ s$^{-1}$]"); ax[2].set(xlabel="T [K]", ylabel="coke selectivity")
ax[0].legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{FIG}/s2_T.png")
pd.set_option("display.width", 220)
for nm, d in [("pCH4", d1), ("pCO2", d2), ("rm", d3), ("rint", d4), ("pH2", d5), ("T", dT)]:
    print("=====", nm)
    print(d[d.case != "no coke (paper)"][["case", "x", "r_DRM", "r_coke", "S_coke", "th_C*", "th_O*", "th_C2*", "ratio_Cform_Otrans"]].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
