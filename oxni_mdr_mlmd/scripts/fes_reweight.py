"""Free-energy surface from an OPES COLVAR file by reweighting, plus a TST rate estimate.

  python fes_reweight.py COLVAR --cv cn --T 873 --skip 0.2 --is 3.9:4.1 --ts-range 3.0:3.9
  python fes_reweight.py COLVAR --cv nco ncc --T 873 --bins 60 60 --png fes2d.png

w_t = exp(V_bias(t)/kT) (all bias columns ending in '.bias' are summed); F = -kT ln P_w.
1D: dF_act = max F in ts-range - min F in is-range; k_TST = kT/h exp(-dF_act/kT).
This TST estimate ignores recrossing; check with committor shooting or OPES-flooding.
"""
import argparse

import numpy as np

KB = 8.617333262e-5  # eV/K
H = 4.135667696e-15  # eV s


def read_colvar(path):
    with open(path) as fh:
        names = fh.readline().split()[2:]
    data = np.loadtxt(path, comments="#")
    return names, data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("colvar")
    ap.add_argument("--cv", nargs="+", required=True)
    ap.add_argument("--T", type=float, default=873.0)
    ap.add_argument("--skip", type=float, default=0.2, help="fraction discarded as transient")
    ap.add_argument("--bins", nargs="+", type=int, default=[80])
    ap.add_argument("--is", dest="is_", default=None, help="lo:hi of the initial-state basin (1D)")
    ap.add_argument("--ts-range", default=None, help="lo:hi searched for the barrier top (1D)")
    ap.add_argument("--png", default=None)
    a = ap.parse_args()
    names, d = read_colvar(a.colvar)
    d = d[int(a.skip * len(d)):]
    kT = KB * a.T
    bias = sum(d[:, i] for i, n in enumerate(names) if n.endswith(".bias") and n.startswith("opes"))
    w = np.exp((bias - bias.max()) / kT)
    X = [d[:, names.index(c)] for c in a.cv]
    if len(X) == 1:
        h, e = np.histogram(X[0], bins=a.bins[0], weights=w, density=True)
        x = 0.5 * (e[1:] + e[:-1])
        F = -kT * np.log(np.where(h > 0, h, np.nan))
        F -= np.nanmin(F)
        for xi, fi in zip(x, F):
            print("%8.3f %8.3f" % (xi, fi))
        if a.is_ and a.ts_range:
            lo, hi = map(float, a.is_.split(":"))
            tlo, thi = map(float, a.ts_range.split(":"))
            f_is = np.nanmin(F[(x >= lo) & (x <= hi)])
            f_ts = np.nanmax(F[(x >= tlo) & (x <= thi)])
            dF = f_ts - f_is
            print("# dF_act = %.3f eV ; k_TST = %.3e s^-1 at %g K" % (dF, kT / H * np.exp(-dF / kT), a.T))
        if a.png:
            import matplotlib.pyplot as plt
            plt.plot(x, F)
            plt.xlabel(a.cv[0])
            plt.ylabel("F (eV)")
            plt.savefig(a.png, dpi=200, bbox_inches="tight")
    else:
        bins = a.bins if len(a.bins) == 2 else [a.bins[0]] * 2
        h, ex, ey = np.histogram2d(X[0], X[1], bins=bins, weights=w, density=True)
        F = -kT * np.log(np.where(h > 0, h, np.nan))
        F -= np.nanmin(F)
        np.savetxt((a.png or "fes2d") + ".dat", np.nan_to_num(F, nan=99.0), fmt="%.4f")
        if a.png:
            import matplotlib.pyplot as plt
            plt.contourf(0.5 * (ex[1:] + ex[:-1]), 0.5 * (ey[1:] + ey[:-1]), F.T, levels=20)
            plt.colorbar(label="F (eV)")
            plt.xlabel(a.cv[0])
            plt.ylabel(a.cv[1])
            plt.savefig(a.png, dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    main()
