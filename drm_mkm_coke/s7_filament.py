"""
Scenario C: growth-type carbon (whisker / filament vs surface carbon islands).

Part 1  full surface+bulk ODE (filament.py) for the two pure routes, to get the
        surface state of a "clean" particle and of a "whisker-growing" particle.
Part 2  particle-population model: every Ni particle is clean (c), whisker-growing
        (w, still active) or encapsulated (e, dead).  Transitions are nucleation
        events with effective barriers:
            c -> w   J_fil = A exp(-E_nuc/kT) (S-1)^2          (filament nucleus at the rear)
            c -> e   J_isl = A exp(-E_isl/kT) Zmm thC*^2 N_surf (graphene island on the surface,
            w -> e   same with the whisker-state coverage        followed by supply-limited growth
                                                                that covers the particle)
        Whisker particles have a much lower C* coverage (carbon drains through the
        bulk), so they nucleate islands much more slowly: the filament protects the
        surface.  Activity(t) = P_c r_c + P_w r_w, filament carbon = P_w r_fil.
"""
from coke_tools import *
from filament import FilamentModel, FilamentParams, CokeParams
from mkmc import network as net
from mkmc.config import R_EV
from mkmc.sites import calculate_sites
from mkmc.solver import steady_rates
KBT_H = 2.0268e13

def states(pp, fp=None):
    fp = fp or FilamentParams()
    kT = R_EV * pp.T
    net.use_variant("blocking")
    rc0 = load_model(coke=False)
    ss, r, rn = steady_rates(rc0.k, pp)
    i = {s: j for j, s in enumerate(net.SPECIES)}
    s = calculate_sites(*pp.as_site_args())
    # Ni atoms per particle from the geometry helper
    r_fcc = 0.3506 / (2 * 2 ** 0.5); n_ni = (2 / 3 * np.pi * pp.rm ** 3) / (4 / 3 * np.pi * r_fcc ** 3)
    n_part = (pp.m_load * pp.cat_wt / 58.693 * 6.022e23) / n_ni
    N_surf = s.NSm / n_part
    m = FilamentModel(pp, fp=fp.with_(islands=False))
    Keq = np.exp(-fp.dH_sol / kT)
    thC, thf = ss[i["C*"]], ss[i["*"]]
    x_c = Keq * thC / thf / (1 + Keq * thC / thf)
    clean = dict(thC=thC, thf=thf, r_DRM=rn[5], S=x_c / m.x_eq, Zmm=s.Zmm, N_surf=N_surf)
    df = m.run(t_end=2e3, n=150); e = df.iloc[-1]
    whisk = dict(thC=e.th_C, thf=e.th_free, r_DRM=e.r_DRM, S=e.S, r_fil=e.r_fil, f=e.f_nuc, x=e.x_bulk)
    return clean, whisk, m

def rates(clean, whisk, T, E_isl, E_nuc):
    kT = R_EV * T
    J_ic = KBT_H * np.exp(-E_isl / kT) * clean["Zmm"] * clean["thC"] ** 2 * clean["N_surf"]
    J_iw = KBT_H * np.exp(-E_isl / kT) * clean["Zmm"] * whisk["thC"] ** 2 * clean["N_surf"]
    J_f = KBT_H * np.exp(-E_nuc / kT) * max(clean["S"] - 1, 0) ** 2
    return J_ic, J_iw, J_f

def population(t, J_ic, J_iw, J_f):
    a = J_f + J_ic
    P_c = np.exp(-a * t)
    P_w = J_f / (a - J_iw) * (np.exp(-J_iw * t) - np.exp(-a * t)) if abs(a - J_iw) > 1e-30 else J_f * t * np.exp(-a * t)
    return P_c, P_w, 1 - P_c - P_w

def metrics(clean, whisk, T, E_isl, E_nuc, t_max=1e7):
    J_ic, J_iw, J_f = rates(clean, whisk, T, E_isl, E_nuc)
    t = np.logspace(-3, np.log10(t_max), 2000)
    P_c, P_w, P_e = population(t, J_ic, J_iw, J_f)
    act = (P_c * clean["r_DRM"] + P_w * whisk["r_DRM"]) / clean["r_DRM"]
    half = t[np.argmax(act < 0.5)] if (act < 0.5).any() else np.inf
    C_fil = np.concatenate([[0], np.cumsum(0.5 * (P_w[1:] + P_w[:-1]) * whisk["r_fil"] * np.diff(t))]) * 12.011e3  # mg C/g
    return dict(J_isl_clean=J_ic, J_isl_whisker=J_iw, J_fil=J_f, whisker_fraction=J_f / (J_f + J_ic),
                t_half_activity=half, t_whisker_nucleation=1 / J_f if J_f > 0 else np.inf, t_encap_clean=1 / J_ic,
                t_encap_whisker=1 / J_iw, fil_rate_w=whisk["r_fil"], act_w=whisk["r_DRM"] / clean["r_DRM"],
                C_fil_at_half=np.interp(half, t, C_fil) if np.isfinite(half) else C_fil[-1]), (t, P_c, P_w, P_e, act, C_fil)

E_ISL, E_NUC = 3.6, 3.4
pp0 = ProcessParams(); c0, w0, m0 = states(pp0)
print("clean state :", {k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in c0.items()})
print("whisker state:", {k: f"{v:.4g}" for k, v in w0.items()})
print(f"x_eq={m0.x_eq:.3e}  k_diff={m0.k_diff:.3e}  kf_sol={m0.kf22:.3e}")

# ---- Part 1 figure: the two pure routes, full ODE
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
for lab, fp, cp in [("surface islands only (encapsulation)", FilamentParams(whisker=False), CokeParams(Ea_cc=2.6, Ea_gas=1.8, dH_g=-0.4)),
                    ("whisker only (bulk route)", FilamentParams(islands=False), None)]:
    m = FilamentModel(pp0, cp=cp, fp=fp); d = m.run(t_end=1e3, n=300)
    ax[0].semilogx(d.t, d.r_DRM / 0.1243, label=lab); ax[1].semilogx(d.t, d.th_C, label=lab)
    if fp.whisker: ax[2].semilogx(d.t, d.C_fil_mgC_g, "--", label="filament carbon (whisker route)")
    else: ax[2].semilogx(d.t, d.C_island_mgC_g, label="island carbon on Ni (encapsulation route)")
    d.to_csv(f"{RES}/s7_ode_{'islands' if fp.whisker is False else 'whisker'}.csv", index=False)
ax[0].set(xlabel="time [s]", ylabel="DRM activity / paper steady state", ylim=(0, 1.3)); ax[0].legend(fontsize=7)
ax[1].set(xlabel="time [s]", ylabel="θ(C*)", yscale="log"); ax[2].set(xlabel="time [s]", ylabel="carbon [mg C / g cat]", yscale="log"); ax[2].legend(fontsize=7)
fig.suptitle("Scenario C, full ODE: encapsulation kills activity; the whisker route keeps the surface active while carbon flows through the particle", fontsize=9)
fig.tight_layout(); fig.savefig(f"{FIG}/s7_pure_routes.png"); plt.close(fig)

# ---- Part 2a: regime map over the two effective nucleation barriers
Ei = np.arange(3.0, 4.21, 0.1); En = np.arange(2.6, 4.01, 0.1)
WF = np.zeros((len(En), len(Ei))); TH = np.zeros_like(WF)
for a, en in enumerate(En):
    for b, ei in enumerate(Ei):
        mt, _ = metrics(c0, w0, pp0.T, ei, en); WF[a, b] = mt["whisker_fraction"]; TH[a, b] = mt["t_half_activity"]
np.savez(f"{RES}/s7_regime_map.npz", E_isl=Ei, E_nuc=En, whisker_fraction=WF, t_half=TH)
from matplotlib.colors import LogNorm
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
m1 = ax[0].pcolormesh(Ei, En, WF, cmap="viridis", vmin=0, vmax=1, shading="nearest"); fig.colorbar(m1, ax=ax[0], label="fraction of particles that grow whiskers")
cs = ax[0].contour(Ei, En, WF, levels=[0.1, 0.5, 0.9], colors="w", linewidths=0.8); ax[0].clabel(cs, fmt="%g", fontsize=7)
ax[0].plot(E_ISL, E_NUC, "r*", ms=10, label="reference pair"); ax[0].legend(fontsize=7, loc="upper left")
ax[0].set(xlabel="E_isl: effective island-nucleation barrier [eV]", ylabel="E_nuc: filament-nucleation barrier [eV]")
m2 = ax[1].pcolormesh(Ei, En, np.clip(TH, 1e-2, 1e7), norm=LogNorm(1e-2, 1e7), cmap="magma", shading="nearest"); fig.colorbar(m2, ax=ax[1], label="t½ of activity [s, intrinsic time]")
ax[1].plot(E_ISL, E_NUC, "c*", ms=10); ax[1].set(xlabel="E_isl [eV]", ylabel="E_nuc [eV]")
fig.suptitle("Scenario C regime map at the base condition (973 K, 0.5/0.5 bar, r_m = 4 nm)", fontsize=10); fig.tight_layout()
fig.savefig(f"{FIG}/s7_regime_map.png"); plt.close(fig)

# ---- Part 2b: example population curves
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
for ei, en, lab in [(3.2, 3.4, "encapsulation-dominated (E_isl=3.2)"), (E_ISL, E_NUC, f"reference (E_isl={E_ISL}, E_nuc={E_NUC})"), (4.0, 3.0, "whisker-dominated (E_isl=4.0, E_nuc=3.0)")]:
    mt, (t, Pc, Pw, Pe, act, Cf) = metrics(c0, w0, pp0.T, ei, en)
    ax[0].semilogx(t, act, label=lab); ax[1].semilogx(t, Pw, label=lab); ax[2].loglog(t, np.clip(Cf, 1e-3, None), label=lab)
ax[0].set(xlabel="time [s]", ylabel="relative DRM activity"); ax[0].legend(fontsize=7)
ax[1].set(xlabel="time [s]", ylabel="fraction of whisker-growing particles"); ax[2].set(xlabel="time [s]", ylabel="filament carbon [mg C / g cat]")
fig.suptitle("Scenario C population model: activity, whisker population and filament carbon", fontsize=10); fig.tight_layout()
fig.savefig(f"{FIG}/s7_population_examples.png"); plt.close(fig)

# ---- Part 2c: operating conditions at the reference pair
rows = []
def add(var, x, pp, fp=None):
    c, w, m = states(pp, fp); mt, _ = metrics(c, w, pp.T, E_ISL, E_NUC)
    rows.append(dict(var=var, x=x, thC_clean=c["thC"], thC_whisker=w["thC"], S_clean=c["S"], x_eq=m.x_eq, r_DRM_clean=c["r_DRM"], **mt))
for v in [0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]: add("pCH4", v, pp_with(pCH4=v))
for v in [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]: add("pCO2", v, pp_with(pCO2=v))
for v in [873.15, 923.15, 973.15, 1023.15, 1073.15]: add("T", v, pp_with(T=v))
for v in [2, 3, 4, 6, 8, 10]: add("rm", v, pp_with(rm=float(v)))
for v in [0.0, 0.005, 0.02, 0.05]: add("pH2", v, pp_with(pH2=v))
for v in [1.0, 1.2, 1.4, 1.6]: add("Ea_sol", v, pp0, FilamentParams(Ea_sol=v))
for v in [0.3, 0.4, 0.5, 0.6, 0.7]: add("dG_gr", v, pp0, FilamentParams(dG_gr=v))
df = pd.DataFrame(rows); df.to_csv(f"{RES}/s7_conditions.csv", index=False)
pd.set_option("display.width", 250)
print(df[["var", "x", "thC_clean", "thC_whisker", "S_clean", "J_isl_clean", "J_fil", "whisker_fraction", "t_half_activity", "t_encap_whisker", "fil_rate_w", "act_w"]].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
fig, ax = plt.subplots(2, 4, figsize=(15, 6.5))
for j, (var, xl, logx) in enumerate([("pCH4", "p$_{CH_4}$ [bar]", False), ("pCO2", "p$_{CO_2}$ [bar]", True), ("T", "T [K]", False), ("rm", "r$_m$ [nm]", False)]):
    g = df[df["var"] == var]
    ax[0, j].plot(g.x, g.whisker_fraction, "o-", label="whisker fraction"); ax[0, j].set(xlabel=xl, ylabel="fraction of particles → whisker", ylim=(0, 1.05))
    ax[1, j].semilogy(g.x, g.t_half_activity.clip(1e-2, 1e7), "o-", label="t½ activity"); ax[1, j].semilogy(g.x, g.t_whisker_nucleation.clip(1e-2, 1e7), "s--", label="1/J_fil")
    ax[1, j].semilogy(g.x, g.t_encap_clean.clip(1e-2, 1e7), "^:", label="1/J_isl (clean)"); ax[1, j].set(xlabel=xl, ylabel="time [s, intrinsic]")
    if logx: ax[0, j].set_xscale("log"); ax[1, j].set_xscale("log")
ax[1, 0].legend(fontsize=7); fig.suptitle(f"Scenario C vs operating conditions (E_isl={E_ISL}, E_nuc={E_NUC} eV)", fontsize=10); fig.tight_layout()
fig.savefig(f"{FIG}/s7_conditions.png"); plt.close(fig)
