"""按重原子上限 N 截取子网络（反应中所有物种都 <= N），统计节点数与反应数。

    python subnet_stats.py output_n12/ -o results/subnet_scale_n12.csv
"""

import argparse
import collections
import csv
import gzip
import os


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="run_enumeration.py 的输出目录")
    ap.add_argument("-o", "--out", help="写出 CSV 的路径（默认只打印）")
    args = ap.parse_args()

    heavy = [int(r["heavy_atoms"]) for r in csv.DictReader(open(os.path.join(args.run_dir, "species.csv")))]
    gz = os.path.join(args.run_dir, "reactions.csv.gz")
    fh = gzip.open(gz, "rt", newline="") if os.path.exists(gz) else \
        open(os.path.join(args.run_dir, "reactions.csv"), newline="")
    rx = collections.Counter()
    for r in csv.DictReader(fh):
        ids = (r["reactant_ids"] + " " + r["product_ids"]).split()
        rx[(max(heavy[int(i)] for i in ids), r["type"])] += 1

    cnt = collections.Counter(heavy)
    rows, nodes, tot = [], 0, collections.Counter()
    for n in range(max(heavy) + 1):
        nodes += cnt[n]
        for t in ("b1f0", "b0f1", "b1f1"):
            tot[t] += rx[(n, t)]
        rows.append({"max_heavy": n, "nodes": nodes, "reactions": sum(tot.values()),
                     "b1f0": tot["b1f0"], "b0f1": tot["b0f1"], "b1f1": tot["b1f1"]})
    for r in rows:
        print("{max_heavy:3d} {nodes:10d} {reactions:12d} {b1f0:10d} {b0f1:10d} {b1f1:12d}".format(**r))
    if args.out:
        with open(args.out, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
