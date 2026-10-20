"""命令行入口：枚举费托合成反应网络。

    python run_enumeration.py -c config.yaml -o output/ [--max-heavy 8] [--nproc 32] [--compact] [--yarp]

输出：
    output/species.csv     物种表
    output/reactions.csv   反应表（有向；b0f1 与其逆反应 b1f0 分别列出）
    output/reactions.csv.gz  （--compact）gzip 反应表，只有 type、reactant_ids、product_ids，
                           行号（从 0 开始）即 rxn_id；大规模运行时用
    output/summary.json    统计
    output/species_yarp.csv  （--yarp）YARP Lewis 结构 SMILES 与 yarpecule hash
"""

import argparse
import copy
import json
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ftc_enum import Network  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", default=os.path.join(os.path.dirname(__file__), "config.yaml"))
    ap.add_argument("-o", "--out", default="output")
    ap.add_argument("--max-heavy", type=int, help="覆盖 config 中的 max_heavy")
    ap.add_argument("--max-branches", type=int, help="覆盖 config 中的 max_branches")
    ap.add_argument("--max-species", type=int, help="覆盖 config 中的 max_species")
    ap.add_argument("--yarp", action="store_true", help="枚举完成后用 YARP 生成 Lewis 结构并校验去重")
    ap.add_argument("--nproc", type=int, default=os.cpu_count(), help="并行进程数（大分子段枚举与 YARP 注释）")
    ap.add_argument("--compact", action="store_true", help="反应表写成 gzip、不带反应 SMILES（大规模运行）")
    ap.add_argument("--log-every", type=int, default=1000)
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    cfg = copy.deepcopy(cfg)
    if args.max_heavy is not None:
        cfg["max_heavy"] = args.max_heavy
        cfg["limits"][-1]["max_heavy"] = max(cfg["limits"][-1]["max_heavy"], args.max_heavy)
    if args.max_branches is not None:
        cfg["max_branches"] = args.max_branches
    if args.max_species is not None:
        cfg["max_species"] = args.max_species

    net = Network(cfg, args.out, nproc=args.nproc, compact=args.compact)
    for smi in cfg["seeds"]:
        net.add_seed(smi)
    summary = net.run(log_every=args.log_every)
    print(json.dumps({k: v for k, v in summary.items() if k != "config"}, indent=2, ensure_ascii=False))

    if args.yarp:
        from ftc_enum.yarp_annotate import annotate
        res = annotate(net.by_id, os.path.join(args.out, "species_yarp.csv"), nproc=args.nproc,
                       smiles=net.smiles)
        print(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
