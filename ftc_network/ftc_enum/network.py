"""反应网络闭包枚举。

从种子物种出发，反复对"新物种"做单分子反应、以及和已处理物种的双分子反应，
直到不再产生新物种（或达到 max_species 安全上限）。

双分子配对规则（config.partners）：
  * A + A（同种两分子）总是计算；
  * 若较小的一方重原子数 <= pool_max_heavy（小分子池：H、H2、CO、HCO、OH、CHx ...），计算；
  * 若两者重原子数都 <= full_max_heavy（小分子段），全部两两计算；
  * 其余（大 + 大）不计算——否则配对数 ~ 物种数^2，无法做到 C25。

每个反应物组合（单个物种，或一对物种）只会被处理一次，所以反应去重只需要在
一次处理内部进行，不需要全局的反应集合（大规模运行时省内存）。

执行分两段（nproc > 1 时）：
  1. 小分子段（重原子 <= full_max_heavy）按"轮"并行：同一轮中第 k 个物种与
     已处理物种以及本轮前 k 个物种配对，每对仍只计算一次；
  2. 大分子只和小分子池 + 自身反应，按"轮"批量交给多进程并行处理。
     若大分子反应又产生了新的小分子，串行处理这些小分子（和所有已处理物种配对），
     再继续下一轮。
nproc = 1 时全部按重原子数顺序串行。最终网络与处理顺序无关。
"""

import csv
import heapq
import json
import os
import time
from multiprocessing import Pool

from .rules import bimolecular, unimolecular
from .species import make_species, species_from_smiles


def _tier(cfg, heavy):
    for t in cfg["limits"]:
        if heavy <= t["max_heavy"]:
            return t
    return None


def make_filter(cfg):
    max_heavy = cfg["max_heavy"]
    forbid_oo = cfg.get("forbid_O_O", True)
    max_branches = cfg.get("max_branches")

    def flt(s):
        if s.heavy > max_heavy:
            return False
        t = _tier(cfg, s.heavy)
        if t is None:
            return False
        if t.get("max_O") is not None and s.nO > t["max_O"]:
            return False
        if t.get("max_free_valence") is not None and s.free_valence > t["max_free_valence"]:
            return False
        if t.get("max_pi") is not None and s.n_pi > t["max_pi"]:
            return False
        mb = t.get("max_branches", max_branches)
        if mb is not None and s.branches > mb:
            return False
        if forbid_oo:
            for (i, j) in s.edges:
                if s.elems[i] == "O" and s.elems[j] == "O":
                    return False
        return True

    return flt


class Reactor:
    """对一个物种 s 生成它的全部反应（单分子 + 与给定 partners 的双分子）。可在子进程中运行。"""

    def __init__(self, cfg):
        self.cfg = cfg
        self.flt = make_filter(cfg)
        r = cfg["reactions"]
        self.do = {k: bool(r.get(k, True)) for k in ("b1f0", "b0f1", "b1f1")}
        self.group_transfer = r.get("b1f1_bimolecular_group_transfer", False)
        self.open_shell_acceptor = r.get("b1f1_bimolecular_open_shell_acceptor", True)
        self.self_rxn = cfg["partners"].get("self_reaction", True)
        self.max_heavy = cfg["max_heavy"]

    def _collect(self, gen, reactants, out):
        r_keys = tuple(sorted(x.key for x in reactants))
        seen = set()
        for rtype, raw in gen:
            prods = []
            for g in raw:
                p = make_species(*g, flt=self.flt)
                if p is None:
                    break
                prods.append(p)
            else:
                p_keys = tuple(sorted(p.key for p in prods))
                if p_keys == r_keys or p_keys in seen:
                    continue  # 简并交换（CH4 + CH3 -> CH3 + CH4）或重复
                seen.add(p_keys)
                out.append((rtype, r_keys, prods))

    def react(self, s, partners):
        out = []
        self._collect(unimolecular(s, b1f0=self.do["b1f0"], b1f1=self.do["b1f1"]), [s], out)
        partners = list(partners)
        if self.self_rxn:
            partners.append(s)
        for b in partners:
            # 结合产物超过重原子上限时，b0f1 不可能有产物，直接跳过
            b0f1 = self.do["b0f1"] and s.heavy + b.heavy <= self.max_heavy
            self._collect(bimolecular(s, b, b0f1=b0f1, b1f1=self.do["b1f1"],
                                      group_transfer=self.group_transfer,
                                      open_shell_acceptor=self.open_shell_acceptor),
                          [s, b], out)
        return out


def eligible(h1, h2, pool_max, full_max):
    """两个物种（重原子数 h1、h2）之间是否计算双分子反应。"""
    return min(h1, h2) <= pool_max or max(h1, h2) <= full_max


# 子进程全局状态
_W = {}


def _worker_init(cfg, pool_species):
    _W["reactor"] = Reactor(cfg)
    _W["pool"] = pool_species


def _worker_react(s):
    return _W["reactor"].react(s, _W["pool"])


def _worker_init_round(cfg, processed, batch):
    _W["reactor"] = Reactor(cfg)
    _W["processed"] = processed
    _W["batch"] = batch
    p = cfg["partners"]
    _W["lim"] = (p["pool_max_heavy"], p["full_max_heavy"])


def _worker_react_round(k):
    """处理本轮第 k 个物种：与已处理物种 + 本轮前 k 个物种配对（每对只算一次）。"""
    s = _W["batch"][k]
    cand = _W["processed"] + _W["batch"][:k]
    partners = [b for b in cand if eligible(s.heavy, b.heavy, *_W["lim"])]
    return _W["reactor"].react(s, partners)


class Network:
    def __init__(self, cfg, out_dir, nproc=1):
        self.cfg = cfg
        self.reactor = Reactor(cfg)
        self.flt = self.reactor.flt
        self.out_dir = out_dir
        self.nproc = max(1, nproc or 1)
        self.species = {}              # key -> Species
        self.by_id = []
        self.smiles = []               # id -> 展示用 SMILES
        self.processed_by_heavy = {}   # heavy -> [Species]
        self.n_rxn = {"b1f0": 0, "b0f1": 0, "b1f1": 0}
        self.queue = []
        p = cfg["partners"]
        self.pool_max = p["pool_max_heavy"]
        self.full_max = p["full_max_heavy"]
        if self.pool_max > self.full_max:
            raise ValueError("partners.pool_max_heavy 不能大于 partners.full_max_heavy")
        self.max_species = cfg.get("max_species")
        os.makedirs(out_dir, exist_ok=True)
        self._rxn_fh = open(os.path.join(out_dir, "reactions.csv"), "w", newline="")
        self._rxn_w = csv.writer(self._rxn_fh)
        self._rxn_w.writerow(["rxn_id", "type", "molecularity", "reactant_ids", "product_ids",
                              "reaction_smiles"])
        self.truncated = False

    # ------------------------------------------------------------ species
    def _register(self, s):
        old = self.species.get(s.key)
        if old is not None:
            return old
        if self.max_species is not None and len(self.by_id) >= self.max_species:
            self.truncated = True
            return None
        s.id = len(self.by_id)
        self.species[s.key] = s
        self.by_id.append(s)
        self.smiles.append(s.smiles())
        heapq.heappush(self.queue, (s.heavy, s.id))
        return s

    def add_seed(self, smi):
        s = species_from_smiles(smi)
        if not self.flt(s):
            raise ValueError(f"种子 {smi} 被过滤条件排除，请检查 limits 设置")
        return self._register(s)

    # ------------------------------------------------------------ reactions
    def _write(self, results):
        for rtype, r_keys, prods in results:
            reg = [self._register(p) for p in prods]
            if any(x is None for x in reg):
                continue
            r_ids = sorted(self.species[k].id for k in r_keys)
            p_ids = sorted(x.id for x in reg)
            rid = sum(self.n_rxn.values())
            self.n_rxn[rtype] += 1
            smi = ".".join(self.smiles[i] for i in r_ids) + ">>" + \
                  ".".join(self.smiles[i] for i in p_ids)
            self._rxn_w.writerow([rid, rtype, len(r_keys), " ".join(map(str, r_ids)),
                                  " ".join(map(str, p_ids)), smi])

    def _partners(self, s):
        out = []
        for h, lst in self.processed_by_heavy.items():
            if eligible(h, s.heavy, self.pool_max, self.full_max):
                out.extend(lst)
        return out

    def _mark_processed(self, s):
        self.processed_by_heavy.setdefault(s.heavy, []).append(s)

    def _pool(self):
        return [s for h, lst in sorted(self.processed_by_heavy.items()) if h <= self.pool_max for s in lst]

    def _log(self, t0, force=False):
        done = sum(len(v) for v in self.processed_by_heavy.values())
        if force or done - self._last_log >= self._log_every:
            self._last_log = done
            print(f"  processed {done:9d} / found {len(self.by_id):9d} species, "
                  f"{sum(self.n_rxn.values()):11d} rxns, {time.time() - t0:9.1f} s", flush=True)

    def run(self, log_every=1000):
        t0 = time.time()
        self._log_every = log_every or float("inf")
        self._last_log = 0
        workers, workers_pool_size = None, -1
        try:
            while self.queue:
                # 0) 还没有处理过大分子时，小分子按轮并行（每轮内部按顺序配对，结果与串行相同）
                if self.nproc > 1 and self.queue[0][0] <= self.full_max and \
                        not any(h > self.full_max for h in self.processed_by_heavy):
                    batch = []
                    while self.queue and self.queue[0][0] <= self.full_max:
                        _, sid = heapq.heappop(self.queue)
                        batch.append(self.by_id[sid])
                    processed = [x for lst in self.processed_by_heavy.values() for x in lst]
                    with Pool(self.nproc, initializer=_worker_init_round,
                              initargs=(self.cfg, processed, batch)) as rp:
                        for s, res in zip(batch, rp.imap(_worker_react_round, range(len(batch)))):
                            self._write(res)
                            self._mark_processed(s)
                            self._log(t0)
                    continue
                # 1) 小分子（以及只用单进程时的所有物种）：串行
                while self.queue and (self.queue[0][0] <= self.full_max or self.nproc == 1):
                    _, sid = heapq.heappop(self.queue)
                    s = self.by_id[sid]
                    self._write(self.reactor.react(s, self._partners(s)))
                    self._mark_processed(s)
                    self._log(t0)
                if not self.queue:
                    break
                # 2) 大分子：本轮所有待处理的大分子一起并行，只和小分子池 + 自身反应
                batch = []
                while self.queue:
                    _, sid = heapq.heappop(self.queue)
                    batch.append(self.by_id[sid])
                pool = self._pool()
                if workers is None or len(pool) != workers_pool_size:
                    if workers is not None:
                        workers.close()
                        workers.join()
                    workers = Pool(self.nproc, initializer=_worker_init, initargs=(self.cfg, pool))
                    workers_pool_size = len(pool)
                for s, res in zip(batch, workers.imap(_worker_react, batch, chunksize=16)):
                    self._write(res)
                    self._mark_processed(s)
                    self._log(t0)
        finally:
            if workers is not None:
                workers.close()
                workers.join()
        if log_every:
            self._log(t0, force=True)
        self._rxn_fh.close()
        self.write_species()
        summary = self.summary(time.time() - t0)
        with open(os.path.join(self.out_dir, "summary.json"), "w") as fh:
            json.dump(summary, fh, indent=2, ensure_ascii=False)
        return summary

    # ------------------------------------------------------------ output
    def write_species(self):
        with open(os.path.join(self.out_dir, "species.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["species_id", "smiles", "formula", "heavy_atoms", "nC", "nH", "nO",
                        "free_valence", "n_pi", "branches", "connectivity_key"])
            for s in self.by_id:
                w.writerow([s.id, self.smiles[s.id], s.formula(), s.heavy, s.nC, s.nH, s.nO,
                            s.free_valence, s.n_pi, s.branches, s.key])

    def summary(self, elapsed):
        by_heavy = {}
        for s in self.by_id:
            by_heavy[s.heavy] = by_heavy.get(s.heavy, 0) + 1
        return {
            "n_species": len(self.by_id),
            "n_reactions": sum(self.n_rxn.values()),
            "reactions_by_type": self.n_rxn,
            "species_by_heavy_atoms": dict(sorted(by_heavy.items())),
            "truncated_by_max_species": self.truncated,
            "elapsed_s": round(elapsed, 1),
            "config": self.cfg,
        }
