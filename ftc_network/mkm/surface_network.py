"""表面费托反应网络（用于随机参数微观动力学）。

物种沿用 ftc_enum 的连接关系表示；"自由价"解释为与表面形成的键数：
CH3* (1)、CH2* (2)、CH* (3)、C* (4)、CO* 写成 [C]=O (2)、HCO* (1) ...
每个吸附物种占 1 个位点；自由价为 0 的含 C 物种（烷烃、1-烯烃、醛、醇）与 H2O
视为生成即离开表面的气相产物；H2 只通过解离吸附进入。

物种空间（surface_ok）：
  * H*、O*、OH*、H2O；
  * C1：任意 H 数、至多 1 个 O，自由价 <= 4；
  * Cn (2 <= n <= n_max)：直链；只有一端（活化端 E）可以有剩余价：
      - 除 E、E 上的 O、以及与 E 相邻的碳（Cβ）外，所有原子饱和；
      - Cβ 至多 1 个剩余价，且必须与 E 形成 π 键（烯基、亚乙烯基一类）；
      - O 至多 1 个，只能是连在 E 上的端基 O（酰基、羟烷基、烷氧基、醛、醇）；
      - 至多 1 个 π 键，自由价 <= 3。

反应（全部可逆，按"结合方向"记录一次）：
  b0f1 / b1f0：A + B <=> C（成/断一个 σ 键），链物种只与 C1、H、O、OH 结合；
  b1f1：分子内 H 迁移 A <=> A'；
  另加 CO(g) + * <=> CO*、H2(g) + 2* <=> 2H*（在 kinetics 中处理）。
"""

from ftc_enum.rules import _join, _split, _whole
from ftc_enum.species import make_species, species_from_smiles


def _active_end_ok(s, E, r):
    el = s.elems
    allowed = {E} | {w for w in s.nbrs[E] if el[w] == "O"}
    beta = [w for w in s.nbrs[E] if el[w] == "C"]
    for i in range(len(el)):
        if i in allowed:
            continue
        if beta and i == beta[0]:
            if r[i] == 0 or (r[i] == 1 and r[E] >= 1):
                continue
            return False
        if r[i] != 0:
            return False
    return True


def surface_ok(s, n_max, fv_max_c1=4, fv_max_chain=3):
    el = s.elems
    if el == ("H",):
        return True
    if s.nC == 0:
        return s.nO == 1
    if s.nC > n_max or s.nO > 1:
        return False
    cs = [i for i, e in enumerate(el) if e == "C"]
    cdeg = {i: sum(1 for w in s.nbrs[i] if el[w] == "C") for i in cs}
    if any(d > 2 for d in cdeg.values()):
        return False
    os_ = [i for i, e in enumerate(el) if e == "O"]
    if os_ and len(s.nbrs[os_[0]]) != 1:
        return False   # 醚
    if s.nC == 1:
        return s.free_valence <= fv_max_c1
    if s.n_pi > 1 or s.free_valence > fv_max_chain:
        return False
    r = [s.residual(i) for i in range(len(el))]
    ends = [i for i in cs if cdeg[i] <= 1]
    if os_:
        E = s.nbrs[os_[0]][0]
        if E not in ends:
            return False
        cands = [E]
    else:
        cands = ends
    return any(_active_end_ok(s, E, r) for E in cands)


def is_gas(s):
    """自由价为 0 的物种（除 H2 外都是离开表面的产物）。"""
    return s.free_valence == 0


class Reaction:
    __slots__ = ("left", "right", "kind", "bond", "cls", "fine")

    def __init__(self, left, right, kind, bond):
        self.left = tuple(sorted(left))     # 结合方向的左边（物种键）
        self.right = tuple(sorted(right))
        self.kind = kind                    # "b0f1" 或 "b1f1"
        self.bond = bond                    # 形成的键，如 "C-H"、"C-C"
        self.cls = None                     # 反应类别（决定本征能垒的分布）
        self.fine = None                    # 细分类别（用于分析）

    def key(self):
        return frozenset([self.left, self.right])


class SurfaceNetwork:
    def __init__(self, n_max=10, h_shift=True):
        self.n_max = n_max
        self.h_shift = h_shift
        self.species = {}
        self.order = []
        self.rxns = {}
        self.co_key = species_from_smiles("[C-]#[O+]").key
        self.h_key = species_from_smiles("[H]").key
        self.h2_key = species_from_smiles("[H][H]").key
        self.h2o_key = species_from_smiles("O").key
        self.c_key = species_from_smiles("[C]").key
        self.o_key = species_from_smiles("[O]").key

    # -------------------------------------------------------------- helpers
    def _mk(self, g):
        s = make_species(*g)
        if s is None or not surface_ok(s, self.n_max):
            return None
        return self.species.get(s.key, s)

    def _add_species(self, s):
        if s.key not in self.species:
            self.species[s.key] = s
            self.order.append(s)
            return True
        return False

    def _add_rxn(self, left, right, kind, bond):
        if kind == "b0f1" and bond == "H-H":
            return   # H2 只通过解离吸附处理
        if sorted(left) == sorted(right):
            return
        r = Reaction(left, right, kind, bond)
        self.rxns.setdefault(r.key(), r)

    @staticmethod
    def _bond(s, i, t, k):
        return "-".join(sorted((s.elems[i], t.elems[k])))

    def is_pool(self, s):
        return s.nC <= 1 and not (is_gas(s))

    # -------------------------------------------------------------- generation
    def build(self, seeds=("[C-]#[O+]", "[H]")):
        queue = []
        for smi in seeds:
            s = species_from_smiles(smi)
            if self._add_species(s):
                queue.append(s)
        processed_pool, processed = [], []
        while queue:
            s = queue.pop(0)
            new = []
            # 气相产物（闭壳层）不再继续反应；它们参与的可逆反应已从另一侧记录
            if not is_gas(s):
                new += self._unimolecular(s)
                partners = processed if self.is_pool(s) else processed_pool
                for b in partners + ([s] if self.is_pool(s) else []):
                    new += self._associate(s, b)
                processed.append(s)
                if self.is_pool(s):
                    processed_pool.append(s)
            for p in new:
                if self._add_species(p):
                    queue.append(p)
        self._classify()
        return self

    def _unimolecular(self, s):
        out = []
        # b1f0：X–H 断裂（逆反应为 S−H + H* -> S）
        for i in range(len(s.elems)):
            if s.hs[i] > 0 and s.elems[i] != "H":
                a = self._mk(_whole(s, {i: -1}))
                h = self.species.get(self.h_key) or species_from_smiles("[H]")
                if a is not None:
                    self._add_rxn([a.key, h.key], [s.key], "b0f1", s.elems[i] + "-H" if s.elems[i] < "H" else "H-" + s.elems[i])
                    out.append(a)
        # b1f0：重原子之间断键
        for cut in s.edges:
            (ei, hi, bi, _), (ej, hj, bj, _) = _split(s, cut)
            a, b = self._mk((ei, hi, bi)), self._mk((ej, hj, bj))
            if a is not None and b is not None:
                self._add_rxn([a.key, b.key], [s.key], "b0f1", self._bond(s, cut[0], s, cut[1]))
                out += [a, b]
        # b1f1：分子内 H 迁移
        if self.h_shift:
            acc = [k for k in range(len(s.elems)) if s.residual(k) > 0]
            for i in range(len(s.elems)):
                if s.hs[i] == 0 or s.elems[i] == "H":
                    continue
                for k in acc:
                    if k != i:
                        a = self._mk(_whole(s, {i: -1, k: +1}))
                        if a is not None and a.key != s.key:
                            self._add_rxn([s.key], [a.key], "b1f1", "H-shift")
                            out.append(a)
        return out

    def _associate(self, s, t):
        out = []
        acc_s = [i for i in range(len(s.elems)) if s.residual(i) > 0]
        acc_t = [k for k in range(len(t.elems)) if t.residual(k) > 0]
        gs, gt = _whole(s), _whole(t)
        for i in acc_s:
            for k in acc_t:
                c = self._mk(_join(gs, gt, i, k))
                if c is not None:
                    self._add_rxn([s.key, t.key], [c.key], "b0f1", self._bond(s, i, t, k))
                    out.append(c)
        return out

    # -------------------------------------------------------------- classification
    def _classify(self):
        S = self.species
        gasC = lambda k: is_gas(S[k]) and S[k].nC >= 1
        for r in self.rxns.values():
            L, R = [S[k] for k in r.left], [S[k] for k in r.right]
            gas_left = any(gasC(k) for k in r.left)
            gas_right = any(gasC(k) for k in r.right)
            fine = None
            if r.kind == "b1f1":
                cls = "H迁移"
            elif r.bond in ("C-H",):
                if gas_right:
                    cls = "含氧产物生成" if R[0].nO else "终止:加氢成烷烃"
                elif gas_left:
                    g = [x for x in L if is_gas(x) and x.nC >= 1][0]
                    cls = "含氧产物生成" if g.nO else "终止:β-H消除成烯烃"
                else:
                    cls = "C1加氢(C–H)" if R[0].nC == 1 else "链端加氢/脱氢(C–H)"
            elif r.bond == "H-O":
                if r.right == (self.h2o_key,):
                    cls = "水生成"
                elif gas_right or gas_left:
                    cls = "含氧产物生成"
                elif R[0].nC == 0:
                    cls = "O加氢(O*→OH*)"
                else:
                    cls = "C–O–H 加氢"
            elif r.bond == "C-O":
                if r.right == (self.co_key,) and set(r.left) == {self.c_key, self.o_key}:
                    cls = "CO直接解离"
                elif gas_right or gas_left:
                    cls = "含氧产物生成"
                elif R[0].nC == 1:
                    cls = "H辅助C–O断裂(C1)"
                else:
                    cls = "链上C–O断裂"
            elif r.bond == "C-C":
                if gas_right:
                    cls = "终止:偶联成产物"
                elif gas_left:
                    cls = "β-C断裂(链裂解)"
                elif self.co_key in r.left:
                    cls = "CO插入"
                    other = [x for x in L if x.key != self.co_key]
                    fine = "CO插入:" + (("C1" if other[0].nC == 1 else "链")) if other else "CO插入:CO+CO"
                elif all(x.nC == 1 for x in L):
                    cls = "C1–C1偶联"
                    fine = "C1–C1偶联:" + "+".join(sorted(x.formula() for x in L))
                else:
                    cls = "链增长"
                    mono = [x for x in L if x.nC == 1]
                    fine = "链增长:+" + (mono[0].formula() if mono else "?")
            else:
                cls = "其他"
            r.cls = cls
            r.fine = fine or cls

    # -------------------------------------------------------------- summary
    def summary(self):
        import collections
        c = collections.Counter(r.cls for r in self.rxns.values())
        n_ads = sum(1 for s in self.order if not is_gas(s))
        return {"n_species": len(self.order), "n_adsorbates": n_ads,
                "n_gas_products": len(self.order) - n_ads, "n_reactions": len(self.rxns),
                "reactions_by_class": dict(c.most_common())}
