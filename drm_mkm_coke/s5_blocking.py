"""Scenario B: deposited carbon stays on Ni sites (Cg*), slowly gasified by O*.
Steady-state blocking and transient deactivation curves."""
from coke_tools import *
from mkmc import network as net
from mkmc.sites import calculate_sites
from scipy.integrate import solve_ivp
net.use_variant("blocking")
EA_CC = 2.2
rows = []
# steady-state blocking vs conditions for two gasification barriers
for Ea_gas in [1.6, 1.8, 2.0]:
    cp = CokeParams(Ea_cc=EA_CC, Ea_gas=Ea_gas); rc = load_model(cp); rc_nc = load_model(coke=False)
    for var, xs, mk in [("pCH4", np.round(np.arange(0.1, 2.01, 0.1), 2), lambda v: pp_with(pCH4=v)),
                        ("rm", np.arange(2, 10.01, 1.0), lambda v: pp_with(rm=v)),
                        ("rint", np.round(np.arange(0.38, 1.91, 0.38), 2), lambda v: pp_with(rsm=v)),
                        ("pCO2", np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]), lambda v: pp_with(pCO2=v))]:
        for v in xs:
            pp = mk(v); o = solve(rc, pp); o0 = solve(rc_nc, pp)
            rows.append(dict(Ea_gas=Ea_gas, var=var, x=v, r_DRM=o["r_DRM"], r_DRM_nocoke=o0["r_DRM"],
                             activity=o["r_DRM"] / o0["r_DRM"], th_Cg=o["th_Cg*"], th_C=o["th_C*"], th_free=o["th_*"],
                             r_dep=2 * o["_rn"][19], r_gas=o["_rn"][20]))
df = pd.DataFrame(rows); df.to_csv(f"{RES}/s5_blocking_steady.csv", index=False)
fig, ax = plt.subplots(2, 4, figsize=(15, 6.5))
for j, var in enumerate(["pCH4", "pCO2", "rm", "rint"]):
    for Ea_gas, g in df[df["var"] == var].groupby("Ea_gas"):
        ax[0, j].plot(g.x, g.th_Cg, marker=".", label=f"Ea_gas={Ea_gas} eV")
        ax[1, j].plot(g.x, g.activity, marker=".", label=f"Ea_gas={Ea_gas} eV")
    ax[0, j].set(xlabel=var, ylabel="steady-state θ(Cg*)", ylim=(0, 1)); ax[1, j].set(xlabel=var, ylabel="activity r/r(no coke)", ylim=(0, 1.1))
    if var == "pCO2": ax[0, j].set_xscale("log"); ax[1, j].set_xscale("log")
ax[0, 0].legend(fontsize=8); fig.suptitle(f"Scenario B (site-blocking coke), Ea_cc={EA_CC} eV: steady-state blocking", fontsize=10); fig.tight_layout()
fig.savefig(f"{FIG}/s5_blocking_steady.png"); plt.close(fig)
pd.set_option("display.width", 200)
print(df[df.Ea_gas == 1.8].to_string(index=False, float_format=lambda v: f"{v:.3g}"))

# transient: integrate from clean surface, Ea_gas = 1.8
cp = CokeParams(Ea_cc=EA_CC, Ea_gas=1.8); rc = load_model(cp)
def transient(pp, t_end=1e3):
    sites = calculate_sites(*pp.as_site_args()); geom = net.build_geom(sites, pp.Zsm, pp.Zms)
    y0 = np.zeros(len(net.SPECIES)); y0[0] = 1; y0[net.SPECIES.index("Ov")] = 0.7; y0[net.SPECIES.index("v")] = 0.3
    norm = np.array([sites.NSm / 6.022e23] * net.N_METAL + [sites.NSs / 6.022e23] * (len(net.SPECIES) - net.N_METAL))
    rhs = lambda t, y: (net.STOICH @ net_rates(rate_vector(y, rc.k, pp.p, geom))) / norm   # site-normalised ODE (true time)
    t_eval = np.logspace(-5, np.log10(t_end), 400)
    res = solve_ivp(rhs, (0, t_end), y0, method="BDF", t_eval=t_eval, rtol=1e-7, atol=1e-11)
    rdrm = np.array([net_rates(rate_vector(y, rc.k, pp.p, geom))[5] for y in res.y.T])
    return res.t, res.y.T, rdrm, sites.NSm / 6.022e23
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8)); trows = []
for lab, pp in [("base", ProcessParams()), ("pCH4=1.0", pp_with(pCH4=1.0)), ("pCH4=1.5", pp_with(pCH4=1.5)), ("rm=8nm", pp_with(rm=8.0)), ("pCO2=0.1", pp_with(pCO2=0.1)), ("rint=1.14", pp_with(rsm=1.14))]:
    t, y, r, nsm = transient(pp)
    thCg = y[:, net.SPECIES.index("Cg*")]
    m = t > 1e-4                                   # skip the fast surface-filling transient
    r0 = r[np.argmax(m)]                           # clean-surface activity
    tof = r0 / nsm                                 # turnovers per surface Ni site per second
    ax[0].semilogx(t[m], thCg[m], label=lab); ax[1].semilogx(t[m], r[m] / r0, label=lab); ax[2].semilogx(t[m] * tof, r[m] / r0, label=lab)
    rel = r / r0; late = m & (rel < 0.5)
    half = t[np.argmax(late)] if late.any() else np.nan
    trows.append(dict(case=lab, r_clean=r0, r_final=r[-1], activity_final=r[-1] / r0, th_Cg_final=thCg[-1], t_half_s=half,
                      turnovers_to_half=half * tof if np.isfinite(half) else np.nan, coke_mgC_per_gcat_final=thCg[-1] * nsm * 12.011e3))
ax[0].set(xlabel="time [s]", ylabel="θ(Cg*)"); ax[1].set(xlabel="time [s]", ylabel="relative DRM activity"); ax[2].set(xlabel="cumulative turnovers per surface Ni site", ylabel="relative DRM activity")
ax[0].legend(fontsize=8); fig.suptitle("Scenario B transient from a clean surface (Ea_cc=2.2, Ea_gas=1.8 eV)", fontsize=10); fig.tight_layout()
fig.savefig(f"{FIG}/s5_blocking_transient.png")
tdf = pd.DataFrame(trows); tdf.to_csv(f"{RES}/s5_blocking_transient.csv", index=False); print(tdf.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
