"""Energy sensitivity of the whisker-vs-encapsulation outcome (scenario C population model)."""
from coke_tools import *
from energy_perturb import perturb, species_to_vary, ts_to_vary
from filament import FilamentModel, FilamentParams, CokeParams
from mkmc import network as net
from mkmc.config import R_EV
from mkmc.sites import calculate_sites
from mkmc.solver import steady_rates
from s7_filament import KBT_H, metrics
net.use_variant("blocking")
E_ISL, E_NUC = 3.6, 3.4
pp = ProcessParams(); kT = R_EV * pp.T
cp = CokeParams(Ea_cc=2.2, Ea_gas=1.8, dH_g=-0.40)
rc_base = load_model(cp)
SP = [s for s in species_to_vary() if s != "Cg*"]; TS = [r for r in ts_to_vary() if r != "R21"]
r_fcc = 0.3506 / (2 * 2 ** 0.5); n_ni = (2 / 3 * np.pi * pp.rm ** 3) / (4 / 3 * np.pi * r_fcc ** 3)
s = calculate_sites(*pp.as_site_args()); n_part = (pp.m_load / 58.693 * 6.022e23) / n_ni; N_surf = s.NSm / n_part
i = {k: j for j, k in enumerate(net.SPECIES)}

def whisker_fraction(ds, dt, fp):
    rc = perturb(rc_base, ds, dt, pp.T)
    k_off = rc.k.copy(); k_off[34:] = 0.0
    ss, r, rn = steady_rates(k_off, pp)
    m = FilamentModel(pp, cp=cp, fp=fp.with_(islands=False)); m.rc.k[:] = rc.k; m.rc.k[40:] = 0.0
    Keq = np.exp(-fp.dH_sol / kT); thC, thf = ss[i["C*"]], ss[i["*"]]
    x_c = Keq * thC / thf / (1 + Keq * thC / thf)
    clean = dict(thC=thC, thf=thf, r_DRM=rn[5], S=x_c / m.x_eq, Zmm=s.Zmm, N_surf=N_surf)
    e = m.run(t_end=2e3, n=120).iloc[-1]
    whisk = dict(thC=e.th_C, thf=e.th_free, r_DRM=e.r_DRM, S=e.S, r_fil=e.r_fil, f=e.f_nuc, x=e.x_bulk)
    mt, _ = metrics(clean, whisk, pp.T, E_ISL, E_NUC)
    return mt["whisker_fraction"], mt["t_half_activity"], clean["S"]

d = 0.01; rows = []
base = whisker_fraction({}, {}, FilamentParams())
print("base: whisker fraction %.3f  t_half %.3g  S %.2f" % base)
for kind, names in [("species", SP), ("ts", TS), ("whisker", ["dH_sol(C_b)", "dG_gr(C_b vs graphite)", "Ea_sol"])]:
    for nm in names:
        vals = []
        for sgn in (+1, -1):
            if kind == "species": v = whisker_fraction({nm: sgn * d}, {}, FilamentParams())
            elif kind == "ts": v = whisker_fraction({}, {nm: sgn * d}, FilamentParams())
            else:
                fp = FilamentParams()
                if nm.startswith("dH_sol"): fp = fp.with_(dH_sol=fp.dH_sol + sgn * d)
                elif nm.startswith("dG_gr"): fp = fp.with_(dG_gr=fp.dG_gr + sgn * d)
                else: fp = fp.with_(Ea_sol=fp.Ea_sol + sgn * d)
                v = whisker_fraction({}, {}, fp)
            vals.append(v)
        (wp, tp, _), (wm, tm, _) = vals
        # logit sensitivity of the whisker fraction and log sensitivity of t_half, per +0.1 eV
        lg = lambda w: np.log(w / (1 - w)) if 0 < w < 1 else np.nan
        rows.append(dict(kind=kind, name=nm, dlogit_wf_per_0p1eV=(lg(wp) - lg(wm)) / (2 * d) * 0.1,
                         dlog10_thalf_per_0p1eV=(np.log10(tp) - np.log10(tm)) / (2 * d) * 0.1))
df = pd.DataFrame(rows); df.to_csv(f"{RES}/s8b_energy_whisker.csv", index=False)
df["a"] = df.dlogit_wf_per_0p1eV.abs()
pd.set_option("display.width", 200)
print(df.sort_values("a", ascending=False).head(15)[["kind", "name", "dlogit_wf_per_0p1eV", "dlog10_thalf_per_0p1eV"]].to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
fig, ax = plt.subplots(1, 2, figsize=(12, 5.5))
for a, col, ttl in zip(ax, ["dlogit_wf_per_0p1eV", "dlog10_thalf_per_0p1eV"], ["Δ logit(whisker fraction) per +0.1 eV", "Δ log10(activity t½) per +0.1 eV"]):
    g = df.copy(); g["a"] = g[col].abs(); g = g.sort_values("a", ascending=False).head(14).iloc[::-1]
    colors = {"ts": "#c44e52", "species": "#4c72b0", "whisker": "#55a868"}
    a.barh(range(len(g)), g[col], color=[colors[k] for k in g.kind]); a.set_yticks(range(len(g)))
    a.set_yticklabels([("TS " if k == "ts" else "E(") + n + ("" if k == "ts" else ")") for k, n in zip(g.kind, g.name)], fontsize=8)
    a.axvline(0, color="k", lw=0.6); a.set_xlabel(ttl)
fig.suptitle("Scenario C: which energies decide whisker vs encapsulation (red TS, blue species, green bulk-carbon energies)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{FIG}/s8b_tornado_whisker.png")
