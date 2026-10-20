"""随机但热力学一致的能量模型。

物种能量（eV，越低越稳定；以气态原子为零点，空位点为 0）：
    E = -( n_CH·D_CH + n_OH·D_OH + n_CC·D_CC + n_CO·D_CO
           + n_CCπ·D_CCπ + n_COπ·D_COπ
           + fv_C·b_C + fv_O·b_O ) + ε_species
  fv_C、fv_O 是 C、O 上的表面键数（自由价）；b_C = E_C/4、b_O = E_O/2 即标度关系。
  锚点：CO* = -(D_CO,gas + E_ads,CO)，H* = -b_H，H2 = -D_HH，H2O = -2·D_OH。
反应能 = 产物能量 - 反应物能量，因此任何闭合路径的能量和都为 0。

过渡态：E_TS = max(E_左, E_右) + E_int，E_int 来自所属类别：
    每个样本先抽一个类均值 μ_c ~ U(lo, hi)，每个反应再加 N(0, σ_rxn)。
"""

import numpy as np

# 键能与结合能参数的取值范围（每个样本均匀抽样，eV）
THERMO_RANGES = {
    "D_CH": (4.10, 4.40),
    "D_CC": (3.50, 3.80),
    "D_CCpi": (2.50, 2.90),
    "D_CO": (3.60, 3.90),
    "D_COpi": (4.20, 4.60),
    "D_OH": (4.60, 4.90),
    "b_C": (1.60, 1.90),      # 每个 C 表面键；C* 结合能 = 4 b_C
    "b_O": (2.60, 3.00),      # 每个 O 表面键；O* 结合能 = 2 b_O
    "b_H": (2.45, 2.70),      # H*；H2 解离吸附能 = D_HH - 2 b_H
    "E_ads_CO": (1.30, 1.80),
}
FIXED = {"D_HH": 4.52, "D_CO_gas": 11.11}

# 各类反应本征能垒（E_TS 高出 max(初态, 终态) 的部分）的类均值范围
BARRIER_RANGES = {
    "CO直接解离": (0.8, 1.6),
    "H辅助C–O断裂(C1)": (0.5, 1.3),
    "链上C–O断裂": (0.5, 1.3),
    "C1加氢(C–H)": (0.3, 0.8),
    "链端加氢/脱氢(C–H)": (0.3, 0.8),
    "O加氢(O*→OH*)": (0.3, 0.9),
    "C–O–H 加氢": (0.3, 0.9),
    "水生成": (0.5, 1.1),
    "C1–C1偶联": (0.4, 1.1),
    "链增长": (0.4, 1.1),
    "CO插入": (0.6, 1.3),
    "终止:加氢成烷烃": (0.4, 0.9),
    "终止:β-H消除成烯烃": (0.4, 1.0),
    "终止:偶联成产物": (0.4, 1.1),
    "β-C断裂(链裂解)": (0.8, 1.6),
    "含氧产物生成": (0.4, 1.0),
    "H迁移": (1.2, 2.0),
    "其他": (0.5, 1.5),
}


def sample_params(rng, sigma_species=0.10, sigma_rxn=0.10):
    p = {k: rng.uniform(*v) for k, v in THERMO_RANGES.items()}
    p.update(FIXED)
    p["mu"] = {c: rng.uniform(*v) for c, v in BARRIER_RANGES.items()}
    p["sigma_species"] = sigma_species
    p["sigma_rxn"] = sigma_rxn
    return p


def bond_counts(s):
    """(n_CH, n_OH, n_CC, n_CO, n_CCpi, n_COpi, fv_C, fv_O)"""
    el = s.elems
    n_ch = sum(h for e, h in zip(el, s.hs) if e == "C")
    n_oh = sum(h for e, h in zip(el, s.hs) if e == "O")
    n_cc = n_co = n_ccpi = n_copi = 0
    for (i, j) in s.edges:
        pair = "".join(sorted((el[i], el[j])))
        extra = s.bond_orders[(i, j)] - 1
        if pair == "CC":
            n_cc += 1
            n_ccpi += extra
        elif pair == "CO":
            n_co += 1
            n_copi += extra
    fv_c = sum(r for e, r in zip(el, s.radicals) if e == "C")
    fv_o = sum(r for e, r in zip(el, s.radicals) if e == "O")
    return n_ch, n_oh, n_cc, n_co, n_ccpi, n_copi, fv_c, fv_o


def species_energies(net, p, rng):
    """返回 {物种键: 能量}。"""
    E = {}
    for s in net.order:
        k = s.key
        if k == net.h_key:
            E[k] = -p["b_H"]
        elif k == net.h2_key:
            E[k] = -p["D_HH"]
        elif k == net.co_key:
            E[k] = -(p["D_CO_gas"] + p["E_ads_CO"])
        else:
            n_ch, n_oh, n_cc, n_co, n_ccpi, n_copi, fv_c, fv_o = bond_counts(s)
            e = -(n_ch * p["D_CH"] + n_oh * p["D_OH"] + n_cc * p["D_CC"] + n_co * p["D_CO"]
                  + n_ccpi * p["D_CCpi"] + n_copi * p["D_COpi"] + fv_c * p["b_C"] + fv_o * p["b_O"])
            if s.free_valence > 0:              # 吸附物种加随机噪声；气相分子不加
                e += rng.normal(0.0, p["sigma_species"])
            E[k] = e
    return E


def ts_energies(net, rxns, E, p, rng):
    """按 rxns 的顺序返回过渡态能量数组与本征能垒数组。"""
    ts, eint = [], []
    for r in rxns:
        el = sum(E[k] for k in r.left)
        er = sum(E[k] for k in r.right)
        e_int = max(0.05, p["mu"][r.cls] + rng.normal(0.0, p["sigma_rxn"]))
        ts.append(max(el, er) + e_int)
        eint.append(e_int)
    return np.array(ts), np.array(eint)
