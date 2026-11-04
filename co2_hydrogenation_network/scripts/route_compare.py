"""Which species / elementary reactions do two jet-fuel routes share?

    python route_compare.py ../data/net_C16_full  ../data/route_compare_full.json
    python route_compare.py ../data/net_C16_shape ../data/route_compare_shape.json --codes ../data/route_codes_shape.npy

Route A  CH4 dry reforming -> syngas -> Fischer-Tropsch -> hydrocracking / hydroisomerisation / oligomerisation
         enters the network at CO* (syngas CO adsorbed) with H2 available; uses every family except the
         methanol-mediated ones (formate/methanol C1 steps, RWGS from CO2, SMS methylation, hydrocarbon pool).
Route B  direct CO2 hydrogenation (RWGS + FT on Fe-K, or oxide/zeolite via methanol): enters at CO2 + H2 and
         may use every family.
Reachability: level-synchronous hypergraph BFS (both directions of every reversible step; arene -> HCP one-way),
H2, H*, H+, H2O, [H], H. always available.
"""
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

src, out = Path(sys.argv[1]), Path(sys.argv[2])
codes_path = sys.argv[sys.argv.index("--codes") + 1] if "--codes" in sys.argv else None

cls, cn, lab = [], [], []
with gzip.open(src / "species.tsv.gz", "rt") as fh:
    next(fh)
    for line in fh:
        i, c, n, s, _ = line.rstrip("\n").split("\t")
        cls.append(c); cn.append(int(n)); lab.append(s)
NS = len(cls)
cls = np.array(cls); cn = np.array(cn); lab = np.array(lab, dtype=object)
pool = cls == "pool"

fam, A = [], []
with gzip.open(src / "reactions.tsv.gz", "rt") as fh:
    next(fh)
    for line in fh:
        _, f, r, p, _ = line.split("\t", 4)
        rs, ps = r.split(","), p.split(",")
        fam.append(f)
        A.append((int(rs[0]), int(rs[1]) if len(rs) > 1 else -1, int(ps[0]), int(ps[1]) if len(ps) > 1 else -1))
A = np.array(A, dtype=np.int64)
fam = np.array(fam)
NR = len(fam)
print("species", NS, "reactions", NR, flush=True)

# C1 steps that belong to the CO2 / methanol entry only
B_ONLY_C1 = {"CO₂", "CO₂^δ−*", "HCOO*", "H₂COO*", "H₂CO*", "CH₃O*", "CH₃OH/SMS", "DME", "烃池 HCP"}
b_only_ids = np.array([i for i in range(NS) if cls[i] == "core" and lab[i] in B_ONLY_C1])
inv = np.isin(A, b_only_ids).any(axis=1)
b_only_rx = np.isin(fam, ["Z:meth", "A:meth", "A:hcp"]) | ((fam == "C1") & inv)
allowed = {"A": ~b_only_rx, "B": np.ones(NR, dtype=bool)}
start = {"A": [i for i in range(NS) if cls[i] == "core" and lab[i] == "CO*"],
         "B": [i for i in range(NS) if cls[i] == "core" and lab[i] == "CO₂"]}
one_way = fam == "A:hcp"


def reach(route):
    ok = allowed[route]
    fw = ok
    bw = ok & ~one_way
    din = np.concatenate([A[fw][:, 0:2], A[bw][:, 2:4]])
    dout = np.concatenate([A[fw][:, 2:4], A[bw][:, 0:2]])
    lev = np.full(NS, -1, dtype=np.int32)
    lev[pool] = 0
    lev[start[route]] = 0
    fired = np.zeros(len(din), dtype=bool)
    L = 0
    while True:
        l0 = np.where(din[:, 0] >= 0, lev[np.maximum(din[:, 0], 0)], 0)
        l1 = np.where(din[:, 1] >= 0, lev[np.maximum(din[:, 1], 0)], 0)
        ready = ~fired & (l0 >= 0) & (l1 >= 0)
        if not ready.any():
            break
        fired |= ready
        o = dout[ready].ravel()
        o = o[o >= 0]
        o = o[lev[o] < 0]
        if not len(o):
            continue
        L += 1
        lev[o] = L
    return lev >= 0, lev


res = {}
R = {}
for route in ("A", "B"):
    reachable, lev = reach(route)
    # a reaction is usable in a route if its family is allowed and one side is fully reachable
    def side_ok(cols):
        m = np.ones(NR, dtype=bool)
        for c in cols:
            v = A[:, c]
            m &= (v < 0) | reachable[np.maximum(v, 0)]
        return m
    usable = allowed[route] & (side_ok([0, 1]) | (side_ok([2, 3]) & ~one_way))
    R[route] = (reachable, usable, lev)
    print(route, "reachable", int(reachable.sum()), "usable reactions", int(usable.sum()), flush=True)

ra, ua, la = R["A"]
rb, ub, lb = R["B"]
real = ~pool
jet = real & (cn >= 8) & (cn <= 16)
by_c = {}
for c in range(0, int(cn.max()) + 1):
    m = real & (cn == c)
    by_c[c] = {"both": int((m & ra & rb).sum()), "A_only": int((m & ra & ~rb).sum()), "B_only": int((m & rb & ~ra).sum()),
               "neither": int((m & ~ra & ~rb).sum())}
by_cls = {}
for c in sorted(set(cls.tolist())):
    m = real & (cls == c)
    if m.any():
        by_cls[c] = {"both": int((m & ra & rb).sum()), "A_only": int((m & ra & ~rb).sum()), "B_only": int((m & rb & ~ra).sum())}
by_fam = {}
for f in sorted(set(fam.tolist())):
    m = fam == f
    by_fam[f] = {"both": int((m & ua & ub).sum()), "A_only": int((m & ua & ~ub).sum()), "B_only": int((m & ub & ~ua).sum()),
                 "total": int(m.sum())}
res = {
    "network": src.name, "n_species": int(real.sum()), "n_reactions": NR,
    "species": {"A": int((real & ra).sum()), "B": int((real & rb).sum()), "both": int((real & ra & rb).sum()),
                "A_only": int((real & ra & ~rb).sum()), "B_only": int((real & rb & ~ra).sum())},
    "jet_C8_C16": {"total": int(jet.sum()), "A": int((jet & ra).sum()), "B": int((jet & rb).sum()), "both": int((jet & ra & rb).sum())},
    "reactions": {"A": int(ua.sum()), "B": int(ub.sum()), "both": int((ua & ub).sum()),
                  "A_only": int((ua & ~ub).sum()), "B_only": int((ub & ~ua).sum())},
    "B_only_species": [str(lab[i]) for i in np.where(real & rb & ~ra)[0][:40]],
    "A_only_species": [str(lab[i]) for i in np.where(real & ra & ~rb)[0][:40]],
    "by_carbon": by_c, "by_class": by_cls, "by_family": by_fam,
    "jet_steps": {"A_median": float(np.median(la[jet & ra])), "B_median": float(np.median(lb[jet & rb]))},
}
out.write_text(json.dumps(res, indent=1, ensure_ascii=False))
print(json.dumps({k: res[k] for k in ("species", "jet_C8_C16", "reactions", "B_only_species", "A_only_species", "jet_steps")}, ensure_ascii=False, indent=1))
if codes_path:
    code = np.full(NR, 3, dtype=np.uint8)           # 0 shared, 1 B only, 2 A only, 3 neither
    code[ua & ub] = 0
    code[ub & ~ua] = 1
    code[ua & ~ub] = 2
    np.save(codes_path, code)
    np.save(codes_path.replace(".npy", "_species.npy"), (ra.astype(np.uint8) | (rb.astype(np.uint8) << 1)))
