"""基元反应规则（按 σ 连接关系计 break / form，和 YARP 的邻接矩阵语义一致）。

b1f0  断 1 个 σ 键：均裂（C–H、C–C、C–O、O–H、H–H）与 β-断裂
      （断键后键级由 Lewis 结构重新分配，例如 丙基 -> 乙烯 + CH3）。
b0f1  成 1 个 σ 键：自由基结合，以及自由基/卡宾加成到 π 键上
      （H + CO -> HCO，CH3 + C2H4 -> 丙基，CH2 + CH2 -> C2H4 ...）。
b1f1  断 1 个、成 1 个 σ 键，且两个键共用一个原子（基团转移）：
      分子内：H 迁移（含酮-烯醇互变）、重原子基团迁移（迁移距离不限于 1,2）；
      分子间：H 抽提 / 歧化 / H2 + R -> RH + H，可选 SH2 取代与重原子片段转移。

每个生成器返回 (rtype, products) 的迭代器，products 是 (elems, hs, edges) 的列表，
由调用方用 make_species 规范化和过滤。
"""


def _component(sp, start, cut):
    """去掉边 cut 后，包含 start 的连通分量（原子下标列表）。"""
    a, b = cut
    comp, stack = {start}, [start]
    while stack:
        v = stack.pop()
        for w in sp.nbrs[v]:
            if (v == a and w == b) or (v == b and w == a) or w in comp:
                continue
            comp.add(w)
            stack.append(w)
    return sorted(comp)


def _sub(sp, atoms, hs_delta=None):
    """取子图，返回 (elems, hs, edges) 以及 原下标 -> 新下标 的映射。"""
    remap = {a: k for k, a in enumerate(atoms)}
    elems = [sp.elems[a] for a in atoms]
    hs = [sp.hs[a] + (hs_delta.get(a, 0) if hs_delta else 0) for a in atoms]
    edges = [(remap[i], remap[j]) for (i, j) in sp.edges if i in remap and j in remap]
    return elems, hs, edges, remap


def _whole(sp, hs_delta=None):
    elems, hs, edges, _ = _sub(sp, range(len(sp.elems)), hs_delta)
    return elems, hs, edges


def _join(g1, g2, i, k):
    """把两个图拼起来，并在 g1 的 i 与 g2 的 k 之间成键。"""
    e1, h1, b1 = g1
    e2, h2, b2 = g2
    off = len(e1)
    return (list(e1) + list(e2), list(h1) + list(h2),
            list(b1) + [(x + off, y + off) for (x, y) in b2] + [(i, k + off)])


def _split(sp, cut):
    """断边 cut，返回 (frag_i, frag_j)，每个是 (elems, hs, edges, remap)。"""
    i, j = cut
    ci = _component(sp, i, cut)
    cj = _component(sp, j, cut)
    return _sub(sp, ci), _sub(sp, cj)


def _acceptors(sp):
    return [k for k in range(len(sp.elems)) if sp.residual(k) > 0]


def _h_donors(sp):
    return [i for i in range(len(sp.elems)) if sp.hs[i] > 0]


# ---------------------------------------------------------------- 单分子

def unimolecular(sp, b1f0=True, b1f1=True):
    if b1f0:
        # X–H 均裂（含 H2 -> H + H）
        for i in _h_donors(sp):
            yield "b1f0", [_whole(sp, {i: -1}), (["H"], [0], [])]
        # 重原子–重原子 断键（均裂 / β-断裂）
        for cut in sp.edges:
            (ei, hi, bi, _), (ej, hj, bj, _) = _split(sp, cut)
            yield "b1f0", [(ei, hi, bi), (ej, hj, bj)]
    if b1f1:
        acc = _acceptors(sp)
        # 分子内 H 迁移：断 i–H，成 k–H
        for i in _h_donors(sp):
            if sp.elems[i] == "H":
                continue
            for k in acc:
                if k != i:
                    yield "b1f1", [_whole(sp, {i: -1, k: +1})]
        # 分子内重原子基团迁移：断 i–j，j 侧接到 i 侧的 k 上（k != i）
        for cut in sp.edges:
            for (a, b) in (cut, cut[::-1]):
                stay = _component(sp, a, cut)   # a 所在的一侧
                move = _component(sp, b, cut)   # b 所在一侧整体迁移
                es, hs_, bs, rs = _sub(sp, stay)
                em, hm, bm, rm = _sub(sp, move)
                for k in stay:
                    if k == a:
                        continue
                    # k 在原分子中的残余价 >0 才能多成一个键
                    if sp.residual(k) <= 0:
                        continue
                    yield "b1f1", [_join((es, hs_, bs), (em, hm, bm), rs[k], rm[b])]


# ---------------------------------------------------------------- 双分子

def _transfer(donor, acceptor, group_transfer=False, open_shell_acceptor=True):
    """donor 失去一个基团，acceptor 接受它。只产生 b1f1。

    group_transfer=False 时只做 H 转移（H 抽提、歧化、H2 + R -> RH + H）；
    True 时还包括 SH2 取代与重原子片段（CH3、OH ...）的分子间转移。
    open_shell_acceptor=True 时接受方必须是开壳层（自由价 > 0）。
    """
    if open_shell_acceptor and acceptor.free_valence == 0:
        return
    acc = _acceptors(acceptor)
    if not acc:
        return
    gacc = _whole(acceptor)
    # H 转移：断 i–H（donor），成 H–k（acceptor）
    for i in _h_donors(donor):
        d_after = _whole(donor, {i: -1})
        for k in acc:
            ea, ha, ba = gacc
            ha = list(ha)
            ha[k] += 1
            yield "b1f1", [d_after, (ea, ha, ba)]
        # SH2 取代：断 i–H，donor 剩余部分在 i 处接到 acceptor 的 k 上，H· 离去
        if group_transfer and donor.elems[i] != "H" and acceptor.elems != ("H",):
            for k in acc:
                yield "b1f1", [_join(d_after, gacc, i, k), (["H"], [0], [])]
    if not group_transfer:
        return
    # 重原子片段转移：断 i–j，j 侧片段接到 acceptor 的 k 上
    for cut in donor.edges:
        for (a, b) in (cut, cut[::-1]):
            stay = _component(donor, a, cut)
            move = _component(donor, b, cut)
            es, hs_, bs, _ = _sub(donor, stay)
            em, hm, bm, rm = _sub(donor, move)
            for k in acc:
                yield "b1f1", [(es, hs_, bs), _join((em, hm, bm), gacc, rm[b], k)]


def bimolecular(a, b, b0f1=True, b1f1=True, group_transfer=False, open_shell_acceptor=True):
    """a + b（a 可以就是 b，表示同种物种两分子反应）。"""
    if b0f1:
        ga, gb = _whole(a), _whole(b)
        acc_b = _acceptors(b)
        for i in _acceptors(a):
            for k in acc_b:
                yield "b0f1", [_join(ga, gb, i, k)]
    if b1f1:
        yield from _transfer(a, b, group_transfer, open_shell_acceptor)
        if a is not b:
            yield from _transfer(b, a, group_transfer, open_shell_acceptor)
