"""Which DFT energies control coking?  Local (Campbell X_RC / X_TRC) + global (Sobol) sensitivity
in energy space, flux model with Ea_cc = 2.2 eV."""
from coke_tools import *
from energy_perturb import perturb, species_to_vary, ts_to_vary
from multiprocessing import Pool
import sys
EA_CC = 2.2
rc0 = load_model(CokeParams(Ea_cc=EA_CC))
SP = species_to_vary(); TS = ts_to_vary()
OBJ = ["S_coke", "r_coke", "r_DRM"]
CONDS = {"base": ProcessParams(), "pCH4=1.0": pp_with(pCH4=1.0), "pCO2=0.1": pp_with(pCO2=0.1)}

def evaluate(args):
    kind, name, d, cond = args
    pp = CONDS[cond]
    rc = perturb(rc0, {name: d} if kind == "species" else {}, {name: d} if kind == "ts" else {}, pp.T)
    try:
        o = solve(rc, pp); return [o[k] for k in OBJ]
    except Exception:
        return [np.nan] * len(OBJ)

if __name__ == "__main__":
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    kT = 8.617e-5 * 973.15
    d = 0.01
    if os.path.exists(f"{RES}/s8_energy_local.csv") and "--local" not in sys.argv:
        loc = pd.read_csv(f"{RES}/s8_energy_local.csv")
    else:
        jobs = [(kind, nm, s * d, cond) for cond in CONDS for kind, names in [("species", SP), ("ts", TS)] for nm in names for s in (+1, -1)]
        with Pool(procs) as pool:
            out = pool.map(evaluate, jobs, chunksize=4)
        rows = []
        for j in range(0, len(jobs), 2):
            kind, nm, _, cond = jobs[j]; yp, ym = np.array(out[j]), np.array(out[j + 1])
            for o, vp, vm in zip(OBJ, yp, ym):
                dln = (np.log(vp) - np.log(vm)) / (2 * d)            # d ln(obj) / dG  [1/eV]
                rows.append(dict(condition=cond, kind=kind, name=nm, objective=o,
                                 X=-dln * kT,                           # Campbell: d ln(obj)/d(-G/kT)
                                 dlog10_per_0p1eV=dln * 0.1 / np.log(10)))
        loc = pd.DataFrame(rows); loc.to_csv(f"{RES}/s8_energy_local.csv", index=False)
    pd.set_option("display.width", 250)
    for cond in CONDS:
        g = loc[(loc.condition == cond) & (loc.objective == "S_coke")].copy()
        g["absX"] = g.X.abs(); g = g.sort_values("absX", ascending=False)
        print(f"\n== {cond}: sensitivity of coke selectivity (X = dln S/d(-G/kT); last column = change in log10 S per +0.1 eV)")
        print(g.head(16)[["kind", "name", "X", "dlog10_per_0p1eV"]].to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
    # tornado plots
    for cond in CONDS:
        fig, ax = plt.subplots(1, 3, figsize=(15, 6))
        for a, o in zip(ax, OBJ):
            g = loc[(loc.condition == cond) & (loc.objective == o)].copy()
            g["absX"] = g.X.abs(); g = g.sort_values("absX", ascending=False).head(14).iloc[::-1]
            colors = ["#c44e52" if k == "ts" else "#4c72b0" for k in g.kind]
            a.barh(range(len(g)), g.dlog10_per_0p1eV, color=colors)
            a.set_yticks(range(len(g))); a.set_yticklabels([("TS " if k == "ts" else "E(") + n + ("" if k == "ts" else ")") for k, n in zip(g.kind, g.name)], fontsize=8)
            a.axvline(0, color="k", lw=0.6); a.set_xlabel(f"Δlog10({o}) per +0.1 eV"); a.set_title(o)
        fig.suptitle(f"Local energy sensitivity, {cond} (red: transition-state energies, blue: species energies)", fontsize=10)
        fig.tight_layout(); fig.savefig(f"{FIG}/s8_tornado_{cond.replace('=', '')}.png"); plt.close(fig)

    # ---- global: Sobol over all energies, +-0.15 eV
    from SALib.sample import sobol as sobol_sample
    from SALib.analyze import sobol as sobol_analyze
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 256
    names = [f"E({s})" for s in SP] + [f"TS {r}" for r in TS]
    problem = dict(num_vars=len(names), names=names, bounds=[[-0.15, 0.15]] * len(names))
    X = sobol_sample.sample(problem, N, calc_second_order=False)
    import signal
    def _timeout(signum, frame): raise TimeoutError
    def f(x):
        ds = dict(zip(SP, x[:len(SP)])); dt = dict(zip(TS, x[len(SP):]))
        rc = perturb(rc0, ds, dt, 973.15)
        signal.signal(signal.SIGALRM, _timeout); signal.alarm(20)       # per-sample limit: stiff outliers are skipped
        try:
            o = solve(rc, CONDS["base"]); return [o["S_coke"], o["r_DRM"]]
        except Exception:
            return [np.nan, np.nan]
        finally:
            signal.alarm(0)
    globals()["f"] = f
    print(f"\nSobol: {len(X)} evaluations")
    with Pool(procs) as pool:
        Y = np.array(pool.map(f, X, chunksize=16))
    ok = np.isfinite(Y).all(axis=1); print("failed:", (~ok).sum()); Y[~ok] = np.nanmedian(Y[ok], axis=0)
    rows = []
    for j, metric in enumerate(["log10 S_coke", "log10 r_DRM"]):
        Si = sobol_analyze.analyze(problem, np.log10(np.clip(Y[:, j], 1e-12, None)), calc_second_order=False, print_to_console=False)
        for i, nm in enumerate(names):
            rows.append(dict(metric=metric, name=nm, S1=Si["S1"][i], S1_conf=Si["S1_conf"][i], ST=Si["ST"][i], ST_conf=Si["ST_conf"][i]))
    sob = pd.DataFrame(rows); sob.to_csv(f"{RES}/s8_energy_sobol_N{N}.csv", index=False)
    pd.DataFrame(X, columns=names).assign(S_coke=Y[:, 0], r_DRM=Y[:, 1]).to_csv(f"{RES}/s8_energy_sobol_raw_N{N}.csv", index=False)
    print("\nSobol total indices, log10 coke selectivity, base condition, all energies +-0.15 eV:")
    print(sob[sob.metric == "log10 S_coke"].sort_values("ST", ascending=False).head(15)[["name", "S1", "ST", "ST_conf"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    fig, ax = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    for a, metric in zip(ax, ["log10 S_coke", "log10 r_DRM"]):
        g = sob[sob.metric == metric]; x = np.arange(len(g)); w = 0.38
        a.bar(x - w / 2, g.S1, w, yerr=g.S1_conf, label="S1", color="#4c72b0", capsize=1.5)
        a.bar(x + w / 2, g.ST, w, yerr=g.ST_conf, label="ST", color="#c44e52", capsize=1.5)
        a.set_ylabel(metric); a.axhline(0, color="k", lw=0.6)
    ax[0].legend(); ax[-1].set_xticks(np.arange(len(names))); ax[-1].set_xticklabels(names, rotation=90, fontsize=8)
    fig.suptitle(f"Sobol indices in energy space (N={N}, {len(X)} runs, every species / TS energy ±0.15 eV, base condition)", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/s8_energy_sobol_N{N}.png")
