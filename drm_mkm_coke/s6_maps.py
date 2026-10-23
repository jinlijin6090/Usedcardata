"""2-D operating maps of coke selectivity (flux model)."""
from coke_tools import *
from matplotlib.colors import LogNorm
for Ea in [2.2, 1.8]:
    rc = load_model(CokeParams(Ea_cc=Ea))
    pch4 = np.round(np.linspace(0.1, 2.0, 16), 3); pco2 = np.round(np.logspace(-2, np.log10(2), 16), 4); rms = np.linspace(2, 10, 12)
    S1 = np.full((len(pco2), len(pch4)), np.nan); R1 = S1.copy(); S2 = np.full((len(rms), len(pch4)), np.nan)
    for i, b in enumerate(pco2):
        for j, a in enumerate(pch4):
            o = solve(rc, pp_with(pCH4=a, pCO2=b)); S1[i, j] = o["S_coke"]; R1[i, j] = o["r_coke"]
    for i, r in enumerate(rms):
        for j, a in enumerate(pch4):
            S2[i, j] = solve(rc, pp_with(pCH4=a, rm=r))["S_coke"]
    np.savez(f"{RES}/s6_maps_Ea{Ea}.npz", pch4=pch4, pco2=pco2, rms=rms, S_pco2=S1, r_pco2=R1, S_rm=S2)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    m = ax[0].pcolormesh(pch4, pco2, np.clip(S1, 1e-6, 1), norm=LogNorm(1e-6, 1), cmap="magma_r", shading="nearest")
    ax[0].set_yscale("log"); ax[0].set(xlabel="p$_{CH_4}$ [bar]", ylabel="p$_{CO_2}$ [bar]", title="coke selectivity (r$_m$=4 nm)")
    cs = ax[0].contour(pch4, pco2, S1, levels=[1e-3, 1e-2, 1e-1], colors="w", linewidths=0.8); ax[0].clabel(cs, fmt="%g", fontsize=7)
    ax[0].plot([0.1, 2], [0.1, 2], "c--", lw=0.8, label="CH$_4$:CO$_2$=1"); ax[0].legend(fontsize=7, loc="lower right")
    fig.colorbar(m, ax=ax[0])
    m2 = ax[1].pcolormesh(pch4, rms, np.clip(S2, 1e-6, 1), norm=LogNorm(1e-6, 1), cmap="magma_r", shading="nearest")
    ax[1].set(xlabel="p$_{CH_4}$ [bar]", ylabel="r$_m$ [nm]", title="coke selectivity (p$_{CO_2}$=0.5 bar)")
    cs = ax[1].contour(pch4, rms, S2, levels=[1e-3, 1e-2, 1e-1], colors="w", linewidths=0.8); ax[1].clabel(cs, fmt="%g", fontsize=7)
    fig.colorbar(m2, ax=ax[1]); fig.suptitle(f"Operating maps, Ea_cc = {Ea} eV, 973 K", fontsize=10); fig.tight_layout()
    fig.savefig(f"{FIG}/s6_maps_Ea{Ea}.png"); plt.close(fig)
    print(f"Ea={Ea}: S_coke at (0.5,0.5)={S1[np.argmin(abs(pco2-0.5)), np.argmin(abs(pch4-0.5))]:.3g}")
