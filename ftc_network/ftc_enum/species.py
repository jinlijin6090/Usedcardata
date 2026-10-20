"""物种表示：只记录重原子之间的 σ 连接 + 每个重原子上的 H 数。

约定
----
* 物种是一棵树（无环），节点为重原子 C / O，H 以计数形式挂在重原子上。
* 氢原子 H· 和 H2 用一个 'H' 节点表示：H· = ('H', h=0)，H2 = ('H', h=1)。
* 残余价 r_i = 标准价 - σ 邻居数 - H 数。在树上用"叶子优先贪心"把残余价
  两两配成 π 键（对树上的 b-matching 是最优的），剩下的就是自由价
  （free valence，自由基电子 + 卡宾孤对，例如 CH3· = 1, CH2 = 2, CO 写成 [C]=O 时 = 2）。
* 键级 / Lewis 结构只用于展示与过滤；物种的身份完全由连接关系决定，
  与 YARP 一致（YARP 的 yarpecule 也是由邻接矩阵唯一确定）。
"""

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

VALENCE = {"C": 4, "O": 2, "H": 1}
_BT = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}


class Species:
    """一个规范化之后的物种。只读，构造请用 make_species()。"""

    __slots__ = ("elems", "hs", "edges", "nbrs", "key", "heavy", "nC", "nO", "nH",
                 "free_valence", "n_pi", "branches", "bond_orders", "radicals", "id")

    def __repr__(self):
        return f"Species({self.smiles()})"

    def residual(self, i):
        return VALENCE[self.elems[i]] - len(self.nbrs[i]) - self.hs[i]

    def formula(self):
        parts = []
        for el, n in (("C", self.nC), ("H", self.nH), ("O", self.nO)):
            if n:
                parts.append(el + (str(n) if n > 1 else ""))
        return "".join(parts)

    def smiles(self):
        """按贪心配出的键级生成的中性 SMILES（自由价写成自由基电子）。"""
        if self.elems == ("H",):
            return "[H]" if self.hs[0] == 0 else "[H][H]"
        m = Chem.RWMol()
        for i, el in enumerate(self.elems):
            a = Chem.Atom(el)
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(self.hs[i])
            a.SetNumRadicalElectrons(self.radicals[i])
            m.AddAtom(a)
        for (i, j) in self.edges:
            m.AddBond(i, j, _BT[self.bond_orders[(i, j)]])
        mol = m.GetMol()
        mol.UpdatePropertyCache(strict=False)
        return Chem.MolToSmiles(mol, canonical=True)


def _connectivity_key(elems, hs, edges):
    """连接关系的规范 SMILES（全部单键），作为物种的唯一键。

    H 数用同位素标签编码（isotope = H 数 + 1），保证未 sanitize 时 SMILES 也区分 H 数。
    """
    m = Chem.RWMol()
    for i, el in enumerate(elems):
        a = Chem.Atom(el)
        a.SetNoImplicit(True)
        a.SetIsotope(hs[i] + 1)
        m.AddAtom(a)
    for (i, j) in edges:
        m.AddBond(i, j, Chem.BondType.SINGLE)
    return Chem.MolToSmiles(m.GetMol(), canonical=True)


def _fold_hydrogens(elems, hs, edges):
    """把与重原子相连的 'H' 节点折叠成 H 计数；H–H 折叠成 H2。"""
    if "H" not in elems:
        return list(elems), list(hs), list(edges)
    if all(e == "H" for e in elems):
        # 只可能是 H· + H· -> H2（或单个 H 节点）
        n_h = len(elems) + sum(hs)
        if n_h == 1:
            return ["H"], [0], []
        if n_h == 2:
            return ["H"], [1], []
        return None
    hs = list(hs)
    drop = set()
    for (i, j) in edges:
        if elems[i] == "H" and elems[j] != "H":
            hs[j] += 1 + hs[i]
            drop.add(i)
        elif elems[j] == "H" and elems[i] != "H":
            hs[i] += 1 + hs[j]
            drop.add(j)
        elif elems[i] == "H" and elems[j] == "H":
            return None
    remap, new_elems, new_hs = {}, [], []
    for i, el in enumerate(elems):
        if i in drop:
            continue
        if el == "H":
            return None  # 孤立的 H 节点混在重原子分子里，不合法
        remap[i] = len(new_elems)
        new_elems.append(el)
        new_hs.append(hs[i])
    new_edges = [(remap[i], remap[j]) for (i, j) in edges if i not in drop and j not in drop]
    return new_elems, new_hs, new_edges


def make_species(elems, hs, edges, flt=None):
    """规范化 + 计算性质 + 过滤。返回 Species，或在不合法/被过滤时返回 None。

    flt: callable(Species) -> bool，在计算规范键之前调用（规范键最贵）。
    """
    folded = _fold_hydrogens(elems, hs, edges)
    if folded is None:
        return None
    elems, hs, edges = folded
    n = len(elems)
    nbrs = [[] for _ in range(n)]
    cdeg = [0] * n
    for (i, j) in edges:
        nbrs[i].append(j)
        nbrs[j].append(i)
        if elems[j] == "C":
            cdeg[i] += 1
        if elems[i] == "C":
            cdeg[j] += 1
    r = [VALENCE[elems[i]] - len(nbrs[i]) - hs[i] for i in range(n)]
    if min(r) < 0:
        return None  # 超价
    if len(edges) != n - 1:
        return None  # 不连通或有环

    # 叶子优先贪心：残余价配成 π 键
    order, parent, seen = [0], {0: None}, {0}
    for v in order:
        for w in nbrs[v]:
            if w not in seen:
                seen.add(w)
                parent[w] = v
                order.append(w)
    if len(order) != n:
        return None
    left = list(r)
    bo = {}
    for (i, j) in edges:
        bo[(min(i, j), max(i, j))] = 1
    n_pi = 0
    for v in reversed(order):
        p = parent[v]
        if p is None:
            continue
        x = min(left[v], left[p], 2)
        if x:
            bo[(min(v, p), max(v, p))] += x
            left[v] -= x
            left[p] -= x
            n_pi += x

    s = Species()
    s.elems = tuple(elems)
    s.hs = tuple(hs)
    s.edges = tuple((min(i, j), max(i, j)) for (i, j) in edges)
    s.nbrs = tuple(tuple(x) for x in nbrs)
    s.bond_orders = bo
    s.radicals = tuple(left)
    s.free_valence = sum(left)
    s.n_pi = n_pi
    s.nC = elems.count("C")
    s.nO = elems.count("O")
    s.heavy = s.nC + s.nO
    s.nH = sum(hs) + elems.count("H")
    s.branches = sum(c - 2 for i, c in enumerate(cdeg) if c > 2 and elems[i] == "C")
    s.id = None
    if flt is not None and not flt(s):
        return None
    s.key = _connectivity_key(s.elems, s.hs, s.edges)
    return s


def species_from_smiles(smi):
    """从 SMILES 读入物种（只取连接关系，键级会被重新分配）。"""
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        raise ValueError(f"无法解析 SMILES: {smi}")
    mol.UpdatePropertyCache(strict=False)
    mol = Chem.AddHs(mol)
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() != "H"]
    if not heavy:
        n_h = mol.GetNumAtoms()
        sp = make_species(["H"], [n_h - 1], [])
        if sp is None:
            raise ValueError(smi)
        return sp
    idx = {a: k for k, a in enumerate(heavy)}
    elems = [mol.GetAtomWithIdx(a).GetSymbol() for a in heavy]
    for el in elems:
        if el not in VALENCE:
            raise ValueError(f"不支持的元素 {el} in {smi}")
    hs = [sum(1 for nb in mol.GetAtomWithIdx(a).GetNeighbors() if nb.GetSymbol() == "H") for a in heavy]
    edges = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in idx and j in idx:
            edges.append((idx[i], idx[j]))
    sp = make_species(elems, hs, edges)
    if sp is None:
        raise ValueError(f"不合法的物种（超价/成环/不连通）: {smi}")
    return sp
