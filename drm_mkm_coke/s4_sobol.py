"""Sobol global sensitivity of coke rate / coke selectivity / DRM rate (flux model, Ea_cc = 2.2 eV)."""
from coke_tools import *
from multiprocessing import Pool
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 256
PROCS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
EA_CC = 2.2
rc0 = load_model(CokeParams(Ea_cc=EA_CC))
K_VARS = {"k2": 2, "k6": 6, "k9": 9, "k11": 11, "k13": 13, "k15": 15, "k16": 16, "k17": 17, "k18": 18, "k20": 20}
names = ["rm", "rint", "Zsm", "Zms", "pCH4", "pH2", "pCO", "pH2O"] + list(K_VARS)
lo, hi = 0.75, 1.25
bounds = [[2, 10], [0.25, 1.0], [1, 4], [1, 4], [lo * 0.5, hi * 0.5], [lo * 0.04, hi * 0.04], [lo * 0.04, hi * 0.04], [lo * 0.02, hi * 0.02]]
for n in K_VARS.values():
    bounds.append([lo * rc0.k[2 * (n - 1)], hi * rc0.k[2 * (n - 1)]])
problem = dict(num_vars=len(names), names=names, bounds=bounds)

def f(x):
    rm, rint, Zsm, Zms, pCH4, pH2, pCO, pH2O, *ks = x
    k = rc0.k.copy()
    for (nm, n), kv in zip(K_VARS.items(), ks):
        fw = 2 * (n - 1); k[fw] = kv
        if np.isfinite(rc0.Keq[n - 1]): k[fw + 1] = kv / rc0.Keq[n - 1]
    rc = RateConstants(k=k, Keq=rc0.Keq, Ea=rc0.Ea, k0=rc0.k0)
    pp = pp_with(pCH4=pCH4, pCO2=pCH4, pH2=pH2, pCO=pCO, pH2O=pH2O, rm=rm, rsm=rint, Zsm=Zsm, Zms=Zms)
    try:
        o = solve(rc, pp); return (o["r_DRM"], o["r_coke"], o["S_coke"])
    except Exception:
        return (np.nan, np.nan, np.nan)

if __name__ == "__main__":
    X = sobol_sample.sample(problem, N, calc_second_order=False)
    print("evaluations:", len(X))
    with Pool(PROCS) as pool:
        Y = np.array(pool.map(f, X, chunksize=16))
    ok = np.isfinite(Y).all(axis=1); print("failed:", (~ok).sum())
    Y[~ok] = np.nanmedian(Y[ok], axis=0)   # keep Saltelli structure intact
    rows = []
    for j, metric in enumerate(["r_DRM", "r_coke", "S_coke"]):
        y = Y[:, j]
        if metric != "r_DRM": y = np.log10(np.clip(y, 1e-12, None))   # log for the heavy-tailed coke outputs
        Si = sobol_analyze.analyze(problem, y, calc_second_order=False, print_to_console=False)
        for i, nm in enumerate(names):
            rows.append(dict(metric=metric, param=nm, S1=Si["S1"][i], S1_conf=Si["S1_conf"][i], ST=Si["ST"][i], ST_conf=Si["ST_conf"][i]))
    df = pd.DataFrame(rows); df.to_csv(f"{RES}/s4_sobol_N{N}.csv", index=False)
    pd.DataFrame(X, columns=names).assign(r_DRM=Y[:, 0], r_coke=Y[:, 1], S_coke=Y[:, 2]).to_csv(f"{RES}/s4_sobol_raw_N{N}.csv", index=False)
    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for a, metric, ttl in zip(ax, ["r_DRM", "r_coke", "S_coke"], ["DRM rate", "log10 coke rate", "log10 coke selectivity"]):
        g = df[df.metric == metric]; x = np.arange(len(g)); w = 0.38
        a.bar(x - w / 2, g.S1, w, yerr=g.S1_conf, label="S1", color="#4c72b0", capsize=2)
        a.bar(x + w / 2, g.ST, w, yerr=g.ST_conf, label="ST", color="#c44e52", capsize=2)
        a.set_ylabel(ttl); a.axhline(0, color="k", lw=0.6)
    ax[0].legend(); ax[-1].set_xticks(np.arange(len(names))); ax[-1].set_xticklabels(names, rotation=0)
    fig.suptitle(f"Sobol indices (N={N}, {len(X)} runs, Ea_cc={EA_CC} eV, ±25% on k and pressures, geometry over full range)", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/s4_sobol_N{N}.png")
    pd.set_option("display.width", 200)
    print(df.pivot(index="param", columns="metric", values=["S1", "ST"]).round(3).to_string())
