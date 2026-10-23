"""Morris screening recomputed from elementary effects of complete steps only
(failed evaluations are dropped instead of being replaced by a median), with
failure diagnostics."""
from coke_tools import *
from energy_perturb import perturb, species_to_vary, ts_to_vary
from multiprocessing import Pool
from SALib.sample import morris as morris_sample
import signal, sys, time
rc0 = load_model(CokeParams(Ea_cc=2.2)); SP = species_to_vary(); TS = ts_to_vary()
NAMES = [f"E({s})" for s in SP] + [f"TS {r}" for r in TS]
PP = ProcessParams(); BASE_SS = solve(rc0, PP)["_ss"]
def _to(s, f): raise TimeoutError
def f(x):
    ds = {s: v for s, v in zip(SP, x[:len(SP)]) if v != 0}; dt = {r: v for r, v in zip(TS, x[len(SP):]) if v != 0}
    rc = perturb(rc0, ds, dt, PP.T); t0 = time.time()
    signal.signal(signal.SIGALRM, _to); signal.alarm(60)
    try:
        o = solve_robust(rc, PP, theta_guess=BASE_SS); return [max(np.log10(max(o["S_coke"], 1e-12)), -6.0), max(np.log10(max(o["r_DRM"], 1e-12)), -6.0), 0, time.time() - t0]
    except TimeoutError:
        return [np.nan, np.nan, 1, time.time() - t0]
    except Exception:
        return [np.nan, np.nan, 2, time.time() - t0]
    finally:
        signal.alarm(0)
if __name__ == "__main__":
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    R = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    D = len(NAMES); levels = 6
    problem = dict(num_vars=D, names=NAMES, bounds=[[-0.10, 0.10]] * D)
    X = morris_sample.sample(problem, N=R, num_levels=levels, seed=1)
    print("evaluations:", len(X)); t0 = time.time(); Y = []
    with Pool(procs) as pool:
        for n, y in enumerate(pool.imap(f, X, chunksize=2)):
            Y.append(y)
            if n % 500 == 0: print(f"  {n}/{len(X)} {time.time() - t0:.0f}s", flush=True)
    Y = np.array(Y)
    np.savez(f"{RES}/s8d_morris_raw.npz", X=X, Y=Y, names=NAMES)
    code = Y[:, 2]; print(f"failed: {int((code > 0).sum())}  (timeouts {int((code == 1).sum())}, solver errors {int((code == 2).sum())}); median solve time {np.median(Y[:, 3]):.2f}s, 95% {np.percentile(Y[:, 3], 95):.1f}s")
    # elementary effects per complete step
    rows = []
    for j, metric in enumerate(["log10 S_coke", "log10 r_DRM"]):
        EE = {nm: [] for nm in NAMES}
        for tr in range(R):
            blk = slice(tr * (D + 1), (tr + 1) * (D + 1))
            Xt, Yt = X[blk], Y[blk, j]
            for s in range(D):
                dx = Xt[s + 1] - Xt[s]; k = int(np.argmax(np.abs(dx)))
                if np.isfinite(Yt[s]) and np.isfinite(Yt[s + 1]):
                    EE[NAMES[k]].append((Yt[s + 1] - Yt[s]) / dx[k] * 0.1)     # effect per +0.1 eV
        for nm in NAMES:
            e = np.array(EE[nm])
            rows.append(dict(metric=metric, name=nm, n_steps=len(e), mu=e.mean() if len(e) else np.nan,
                             mu_star=np.abs(e).mean() if len(e) else np.nan, sigma=e.std() if len(e) else np.nan,
                             mu_star_conf=1.96 * np.abs(e).std() / np.sqrt(len(e)) if len(e) > 1 else np.nan))
    mo = pd.DataFrame(rows); mo.to_csv(f"{RES}/s8d_morris_clean.csv", index=False)
    pd.set_option("display.width", 200)
    for metric in ["log10 S_coke", "log10 r_DRM"]:
        print(f"\nMorris (clean), {metric}: mean |effect| per +0.1 eV, averaged over the ±0.10 eV cube")
        print(mo[mo.metric == metric].sort_values("mu_star", ascending=False).head(14)[["name", "n_steps", "mu", "mu_star", "sigma", "mu_star_conf"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for a, metric in zip(ax, ["log10 S_coke", "log10 r_DRM"]):
        g = mo[mo.metric == metric]
        a.scatter(g.mu_star, g.sigma, s=18, color=["#c44e52" if n.startswith("TS") else "#4c72b0" for n in g.name])
        for _, r in g.sort_values("mu_star", ascending=False).head(10).iterrows():
            a.annotate(r["name"], (r.mu_star, r.sigma), fontsize=7, xytext=(3, 3), textcoords="offset points")
        a.set(xlabel=f"μ*: mean |Δ{metric}| per +0.1 eV", ylabel="σ of the elementary effects (interaction / nonlinearity)")
    fig.suptitle("Morris screening of all species (blue) and transition-state (red) energies, ±0.10 eV, base condition", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/s8c_morris.png")
