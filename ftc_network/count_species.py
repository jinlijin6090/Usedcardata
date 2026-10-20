"""在给定 config 的过滤条件下，直接数出每个重原子数 N 的合法物种数（不枚举反应）。

用于在跑完整网络之前估计规模：网络闭包得到的物种集合就是满足过滤条件、
且能从种子到达的物种，通常等于这里数出来的集合。

    python count_species.py -c config.yaml --max-heavy 12 [--max-branches 0]
"""

import argparse
import copy
import os
import sys
import time

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ftc_enum.network import _tier, make_filter  # noqa: E402
from ftc_enum.species import VALENCE, _connectivity_key, make_species  # noqa: E402


def grow_trees(nmax, max_O, max_branches, forbid_oo):
    """按原子逐个生长重原子树（C/O），规范键去重。返回 {N: [(elems, edges)]}。

    max_O(n)、max_branches(n) 给出 n 个重原子时的上限。O 数和支链数随加原子单调不减，
    所以只要上限随 n 不增（config 的默认设置就是这样），逐层剪枝就是安全的。
    """
    def key(elems, edges):
        return _connectivity_key(elems, [0] * len(elems), edges)

    def ok(elems, edges):
        deg = [0] * len(elems)
        cdeg = [0] * len(elems)
        for (i, j) in edges:
            deg[i] += 1
            deg[j] += 1
            if elems[j] == "C":
                cdeg[i] += 1
            if elems[i] == "C":
                cdeg[j] += 1
            if forbid_oo and elems[i] == "O" and elems[j] == "O":
                return False
        if any(deg[i] > VALENCE[e] for i, e in enumerate(elems)):
            return False
        mb = max_branches(len(elems))
        if mb is not None:
            br = sum(max(0, cdeg[i] - 2) for i, e in enumerate(elems) if e == "C")
            if br > mb:
                return False
        return True

    layers = {1: [(["C"], [])] + ([(["O"], [])] if max_O(1) >= 1 else [])}
    for n in range(2, nmax + 1):
        seen = {}
        for elems, edges in layers[n - 1]:
            nO = elems.count("O")
            for i in range(len(elems)):
                for el in ("C", "O"):
                    if el == "O" and nO + 1 > max_O(n):
                        continue
                    ne, nb = elems + [el], edges + [(i, len(elems))]
                    if not ok(ne, nb):
                        continue
                    k = key(ne, nb)
                    if k not in seen:
                        seen[k] = (ne, nb)
        layers[n] = list(seen.values())
    return layers


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--config", default=os.path.join(os.path.dirname(__file__), "config.yaml"))
    ap.add_argument("--max-heavy", type=int, default=10)
    ap.add_argument("--max-branches", type=int)
    args = ap.parse_args()
    with open(args.config) as fh:
        cfg = copy.deepcopy(yaml.safe_load(fh))
    cfg["max_heavy"] = args.max_heavy
    if args.max_branches is not None:
        cfg["max_branches"] = args.max_branches
    cfg["limits"][-1]["max_heavy"] = max(cfg["limits"][-1]["max_heavy"], args.max_heavy)
    flt = make_filter(cfg)

    def max_O(n):
        t = _tier(cfg, n)
        return t["max_O"] if t and t.get("max_O") is not None else n

    t0 = time.time()
    def max_br(n):
        t = _tier(cfg, n)
        return t.get("max_branches", cfg.get("max_branches")) if t else cfg.get("max_branches")

    trees = grow_trees(args.max_heavy, max_O, max_br, cfg.get("forbid_O_O", True))
    total = 2  # H, H2
    print(f"N= 0 species=       2 cumulative=         2")
    for n in range(1, args.max_heavy + 1):
        t = _tier(cfg, n)
        fv_max = t.get("max_free_valence")
        pi_max = t.get("max_pi")
        r_max = None if fv_max is None or pi_max is None else fv_max + 2 * pi_max
        keys = set()
        for elems, edges in trees[n]:
            deg = [0] * n
            for (i, j) in edges:
                deg[i] += 1
                deg[j] += 1
            caps = [VALENCE[e] - deg[i] for i, e in enumerate(elems)]
            budget = r_max if r_max is not None else sum(caps)
            # 枚举残余价分布 r（sum r <= budget），H 数 = 剩余价
            for r in _residuals(caps, budget):
                hs = [caps[i] - r[i] for i in range(n)]
                s = make_species(elems, hs, edges, flt=flt)
                if s is not None:
                    keys.add(s.key)
        total += len(keys)
        print(f"N={n:2d} species={len(keys):8d} cumulative={total:10d}  ({time.time() - t0:.0f} s)", flush=True)


def _residuals(caps, budget):
    """所有满足 0<=r_i<=caps_i 且 sum(r)<=budget 的向量（稀疏枚举）。"""
    n = len(caps)
    yield [0] * n
    for total in range(1, budget + 1):
        # 把 total 个单位分到原子上：先选出现的位置组合，再分配
        for parts in _compositions(total, n, caps):
            yield parts


def _compositions(total, n, caps):
    def rec(i, left, cur):
        if left == 0:
            yield cur + [0] * (n - i)
            return
        if i == n:
            return
        for v in range(min(left, caps[i]), -1, -1):
            yield from rec(i + 1, left - v, cur + [v])
    yield from rec(0, total, [])



if __name__ == "__main__":
    main()
