"""Compare runs (e.g. MEO-confined vs bare Ni) from run_md.py logs or trajectories.

  python plot_runs.py ../runs/S2_meo_log.csv ../runs/S2_bare_log.csv --labels OxNi bare \
      --cols "CO(g)" "CO*" C_in_Cn largest_C_cluster C_dissolved --out ../runs/S2_compare.png
  python plot_runs.py ../runs/S1_873K.extxyz --cols "CH4(g)" "CH3*" "H2(g)" "CO(g)" --out s1.png
A .extxyz input is re-analysed frame by frame with species.py (use after changing cutoffs).
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def load(path):
    if path.endswith(".csv"):
        rows = list(csv.DictReader(open(path)))
        return {k: np.array([float(r[k]) if r[k] not in ("", None) else np.nan for r in rows])
                for k in rows[0] if k not in ("stage", "other_species")}
    from ase.io import iread
    import species as sp
    rows = []
    for f in iread(path):
        spc, met = sp.analyse(f)
        rows.append(dict(time_ps=f.info.get("time_ps", np.nan), **met, **spc))
    keys = set().union(*rows)
    return {k: np.array([r.get(k, 0) for r in rows], float) for k in keys}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--labels", nargs="+", default=None)
    ap.add_argument("--cols", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data = [load(p) for p in a.inputs]
    labels = a.labels or [os.path.basename(p) for p in a.inputs]
    n = len(a.cols)
    fig, axs = plt.subplots(n, 1, figsize=(6, 1.9 * n), sharex=True, squeeze=False)
    for ax, col in zip(axs[:, 0], a.cols):
        for d, lab in zip(data, labels):
            if col in d:
                ax.plot(d["time_ps"], d[col], lw=1.2, label=lab)
        ax.set_ylabel(col)
    axs[0, 0].legend(frameon=False, fontsize=8)
    axs[-1, 0].set_xlabel("time (ps)")
    fig.tight_layout()
    fig.savefig(a.out, dpi=200)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
