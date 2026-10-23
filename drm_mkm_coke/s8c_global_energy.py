"""Global energy sensitivity done properly for a 34-dimensional, heavy-tailed problem:
(1) Morris elementary-effects screening over all species / TS energies (+-0.10 eV),
(2) Sobol over the energies that survive screening."""
from coke_tools import *
from energy_perturb import perturb, species_to_vary, ts_to_vary
from multiprocessing import Pool
from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
import signal, sys, time
rc0 = load_model(CokeParams(Ea_cc=2.2)); SP = species_to_vary(); TS = ts_to_vary()
NAMES = [f"E({s})" for s in SP] + [f"TS {r}" for r in TS]
PP = ProcessParams(); BASE_SS = solve(rc0, PP)["_ss"]
def _to(s, f): raise TimeoutError
def f(x):
    ds = {s: v for s, v in zip(SP, x[:len(SP)]) if v != 0}; dt = {r: v for r, v in zip(TS, x[len(SP):]) if v != 0}
    rc = perturb(rc0, ds, dt, PP.T)
    signal.signal(signal.SIGALRM, _to); signal.alarm(40)
    try:
        o = solve_robust(rc, PP, theta_guess=BASE_SS); return [max(np.log10(max(o["S_coke"], 1e-12)), -6.0), max(np.log10(max(o["r_DRM"], 1e-12)), -6.0)]
    except Exception:
        return [np.nan, np.nan]
    finally:
        signal.alarm(0)
def f_cont(x):
    """Continuation: ramp the perturbation 0 -> 1 in steps, each solved from the previous steady state."""
    guess = BASE_SS
    signal.signal(signal.SIGALRM, _to); signal.alarm(120)
    try:
        for lam in (0.25, 0.5, 0.75, 1.0):
            xl = lam * np.asarray(x)
            ds = {s: v for s, v in zip(SP, xl[:len(SP)]) if v != 0}; dt = {r: v for r, v in zip(TS, xl[len(SP):]) if v != 0}
            o = solve_robust(perturb(rc0, ds, dt, PP.T), PP, theta_guess=guess); guess = o["_ss"]
        return [max(np.log10(max(o["S_coke"], 1e-12)), -6.0), max(np.log10(max(o["r_DRM"], 1e-12)), -6.0)]
    except Exception:
        return [np.nan, np.nan]
    finally:
        signal.alarm(0)

def run(X, procs):
    Y = []; t0 = time.time()
    with Pool(procs) as pool:
        for n, y in enumerate(pool.imap(f, X, chunksize=4)):
            Y.append(y)
            if n % 500 == 0: print(f"  {n}/{len(X)} {time.time() - t0:.0f}s", flush=True)
    return np.array(Y)
def f_sub(x, idx):
    full = np.zeros(len(NAMES)); full[idx] = x; return f(full)

if __name__ == "__main__":
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    # ---- (1) screening: use the clean Morris table from s8d if present, else run Morris here
    problem = dict(num_vars=len(NAMES), names=NAMES, bounds=[[-0.10, 0.10]] * len(NAMES))
    if os.path.exists(f"{RES}/s8d_morris_clean.csv"):
        mo = pd.read_csv(f"{RES}/s8d_morris_clean.csv")
    else:
        X = morris_sample.sample(problem, N=80, num_levels=6, seed=1)
        print("Morris evaluations:", len(X)); Y = run(X, procs)
        ok = np.isfinite(Y).all(axis=1); print("failed:", (~ok).sum())
        rows = []
        for j, metric in enumerate(["log10 S_coke", "log10 r_DRM"]):
            y = Y[:, j].copy(); y[~ok] = np.nanmedian(y[ok])
            Si = morris_analyze.analyze(problem, X, y, conf_level=0.95, num_levels=6, print_to_console=False)
            for i, nm in enumerate(NAMES):
                rows.append(dict(metric=metric, name=nm, mu=Si["mu"][i], mu_star=Si["mu_star"][i], sigma=Si["sigma"][i], mu_star_conf=Si["mu_star_conf"][i]))
        mo = pd.DataFrame(rows); mo.to_csv(f"{RES}/s8c_morris.csv", index=False)
    pd.set_option("display.width", 200)
    top = mo[mo.metric == "log10 S_coke"].sort_values("mu_star", ascending=False)
    print(top.head(12)[["name", "mu", "mu_star", "sigma"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # ---- (2) Sobol over the screened energies
    keep = top.head(10)["name"].tolist(); idx = [NAMES.index(n) for n in keep]
    sub = dict(num_vars=len(keep), names=keep, bounds=[[-0.10, 0.10]] * len(keep))
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 512
    Xs = sobol_sample.sample(sub, N, calc_second_order=True, seed=2)
    print(f"\nSobol on {len(keep)} screened energies: {len(Xs)} evaluations")
    from functools import partial
    globals()["g_sub"] = partial(f_sub, idx=idx)
    Y = []; t0 = time.time()
    with Pool(procs) as pool:
        for n, y in enumerate(pool.imap(g_sub, Xs, chunksize=4)):
            Y.append(y)
            if n % 1000 == 0: print(f"  {n}/{len(Xs)} {time.time() - t0:.0f}s", flush=True)
    Y = np.array(Y); ok = np.isfinite(Y).all(axis=1); print("failed in first pass:", (~ok).sum())
    bad = np.where(~ok)[0]
    if len(bad):
        Xfull = np.zeros((len(bad), len(NAMES))); Xfull[:, idx] = Xs[bad]
        with Pool(procs) as pool:
            Yrep = np.array(pool.map(f_cont, Xfull, chunksize=1))
        Y[bad] = Yrep
        ok = np.isfinite(Y).all(axis=1); print("failed after continuation repair:", (~ok).sum())
    np.savez(f"{RES}/s8c_sobol_raw_N{N}.npz", X=Xs, Y=Y, names=keep)
    rows = []
    for j, metric in enumerate(["log10 S_coke", "log10 r_DRM"]):
        y = Y[:, j].copy(); y[~ok] = np.nanmedian(y[ok])
        Si = sobol_analyze.analyze(sub, y, calc_second_order=True, print_to_console=False)
        for i, nm in enumerate(keep):
            rows.append(dict(metric=metric, name=nm, S1=Si["S1"][i], S1_conf=Si["S1_conf"][i], ST=Si["ST"][i], ST_conf=Si["ST_conf"][i]))
        if metric == "log10 S_coke":
            S2 = pd.DataFrame(Si["S2"], index=keep, columns=keep); S2.to_csv(f"{RES}/s8c_sobol_S2_Scoke.csv")
            pairs = [(keep[a], keep[b], Si["S2"][a, b]) for a in range(len(keep)) for b in range(a + 1, len(keep))]
            pairs.sort(key=lambda t: -abs(t[2])); print("\nstrongest pairwise interactions (S2) for log10 coke selectivity:")
            for a, b, v in pairs[:6]: print(f"   {a:12s} x {b:12s}  S2 = {v:+.3f}")
    so = pd.DataFrame(rows); so.to_csv(f"{RES}/s8c_sobol_screened_N{N}.csv", index=False)
    print("\nSobol indices, log10 coke selectivity:")
    print(so[so.metric == "log10 S_coke"].sort_values("ST", ascending=False)[["name", "S1", "S1_conf", "ST", "ST_conf"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("sum S1 =", so[so.metric == "log10 S_coke"].S1.sum().round(3))
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for a, metric in zip(ax, ["log10 S_coke", "log10 r_DRM"]):
        g = so[so.metric == metric]; x = np.arange(len(g)); w = 0.38
        a.bar(x - w / 2, g.S1, w, yerr=g.S1_conf, label="S1", color="#4c72b0", capsize=2)
        a.bar(x + w / 2, g.ST, w, yerr=g.ST_conf, label="ST", color="#c44e52", capsize=2)
        a.set_ylabel(metric); a.axhline(0, color="k", lw=0.6)
    ax[0].legend(); ax[-1].set_xticks(np.arange(len(keep))); ax[-1].set_xticklabels(keep, rotation=30, fontsize=8)
    fig.suptitle(f"Sobol indices over the {len(keep)} screened energies (N={N}, {len(Xs)} runs, ±0.10 eV, base condition)", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/s8c_sobol_screened.png")
