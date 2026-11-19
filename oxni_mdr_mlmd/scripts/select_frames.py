"""Active learning: pick MD frames for DFT labelling by committee force disagreement.

  python select_frames.py ../runs/S1_873K.extxyz --models mace:ft_s1.model mace:ft_s2.model \
      mace:ft_s3.model --lo 0.10 --hi 0.60 --n 40 --every 5 --out ../al/iter1_candidates.extxyz

sigma_i = sqrt( mean_m |F_i^m - <F_i>|^2 ) per atom; a frame is a candidate when
lo < max_i sigma_i < hi (below lo: already learnt; above hi: probably unphysical).
Candidates are spread over the trajectory (greedy, at least `--gap` frames apart, highest
sigma first). Each written frame carries sigma_max, the index of the most uncertain atom
(atom_sigma_max) and per-atom sigma, so you can also cut a smaller DFT cell around it.
Use a committee of fine-tuned models with different seeds; two different foundation
models (e.g. mace:small + mace:medium) only work as a smoke test.
"""
import argparse
import json
import os
import sys

import numpy as np
from ase.io import read, write

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calculators import get_calculator  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("traj")
    ap.add_argument("--models", nargs="+", required=True)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--dtype", default="float64")
    ap.add_argument("--every", type=int, default=1, help="use every k-th frame")
    ap.add_argument("--lo", type=float, default=0.10)
    ap.add_argument("--hi", type=float, default=0.60)
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--gap", type=int, default=2)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    frames = read(a.traj, index="::%d" % a.every)
    calcs = [get_calculator(m, a.device, a.dtype) for m in a.models]
    sig_max, sig_atoms, arg = [], [], []
    for f in frames:
        F = []
        for c in calcs:
            x = f.copy()
            x.calc = c
            F.append(x.get_forces())
        F = np.array(F)
        s = np.sqrt(((F - F.mean(axis=0)) ** 2).sum(axis=2).mean(axis=0))
        sig_atoms.append(s)
        sig_max.append(float(s.max()))
        arg.append(int(s.argmax()))
    sig_max = np.array(sig_max)

    cand = np.where((sig_max > a.lo) & (sig_max < a.hi))[0]
    chosen = []
    for k in cand[np.argsort(-sig_max[cand])]:
        if all(abs(k - c) >= a.gap for c in chosen):
            chosen.append(int(k))
        if len(chosen) == a.n:
            break
    out = []
    for k in sorted(chosen):
        f = frames[k].copy()
        f.info.update(sigma_max=sig_max[k], atom_sigma_max=arg[k], src_frame=int(k * a.every))
        f.arrays["sigma_F"] = sig_atoms[k]
        out.append(f)
    if out:
        write(a.out, out)
    stats = dict(n_frames=len(frames), n_candidates=int(len(cand)), n_written=len(out),
                 sigma_max_min=float(sig_max.min()), sigma_max_median=float(np.median(sig_max)),
                 sigma_max_max=float(sig_max.max()), frac_below_lo=float((sig_max <= a.lo).mean()),
                 frac_above_hi=float((sig_max >= a.hi).mean()), models=a.models)
    json.dump(stats, open(os.path.splitext(a.out)[0] + "_stats.json", "w"), indent=1)
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
