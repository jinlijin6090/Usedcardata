"""用 YARP 为每个物种生成 Lewis 结构（键级 / 形式电荷 / 自由基），并用 yarpecule hash 校验去重。

物种只由连接关系决定，这里把连接关系转成 YARP 的邻接矩阵输入
(adj_mat, geo, elements, q, atom_info)，由 YARP 的 lewis_struct 给出最优 Lewis 结构。
"""

import csv
from multiprocessing import Pool

import numpy as np


def _adjacency(elems, hs, edges):
    el = list(elems)
    bonds = list(edges)
    for i, h in enumerate(hs):
        for _ in range(h):
            el.append("H")
            bonds.append((i, len(el) - 1))
    n = len(el)
    adj = np.zeros((n, n), dtype=int)
    for i, j in bonds:
        adj[i, j] = adj[j, i] = 1
    return adj, el


def _one(args):
    sid, elems, hs, edges = args
    from yarp.yarpecule.yarpecule import yarpecule
    try:
        adj, el = _adjacency(elems, hs, edges)
        y = yarpecule((adj, np.zeros((len(el), 3)), el, 0, {}))
        y.get_smiles()
        bmat = y.bond_mats[0]
        n_rad = int(sum(int(bmat[i][i]) % 2 for i in range(len(bmat))))
        fc = [int(x) for x in y.fc]
        return sid, y.canon_smi, str(y.hash), n_rad, int(sum(abs(x) for x in fc)), ""
    except Exception as e:  # 个别结构 YARP 无法处理时记录下来，不中断
        return sid, "", "", -1, -1, repr(e)[:200]


def annotate(species, path, nproc=1, smiles=None):
    from .species import unpack
    jobs = []
    for s in species:
        if s.elems is None:   # 已处理的大分子瘦身过，从紧凑表示还原连接关系
            elems, hs, edges = unpack(s.packed)
            jobs.append((s.id, elems, hs, edges))
        else:
            jobs.append((s.id, s.elems, s.hs, s.edges))
    if nproc > 1:
        with Pool(nproc) as pool:
            rows = pool.map(_one, jobs, chunksize=64)
    else:
        rows = [_one(j) for j in jobs]
    rows.sort()
    hash_owner, dup, failed = {}, [], 0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species_id", "smiles", "yarp_smiles", "yarp_hash", "yarp_n_radical_e",
                    "yarp_sum_abs_formal_charge", "error"])
        for (sid, ysmi, yh, nrad, fc, err) in rows:
            if err:
                failed += 1
            elif yh in hash_owner:
                dup.append((hash_owner[yh], sid))
            else:
                hash_owner[yh] = sid
            w.writerow([sid, smiles[sid] if smiles else "", ysmi, yh, nrad, fc, err])
    return {"yarp_annotated": len(rows) - failed, "yarp_failed": failed,
            "yarp_hash_duplicates": dup[:20], "n_yarp_hash_duplicates": len(dup)}
