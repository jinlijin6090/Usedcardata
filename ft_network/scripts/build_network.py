#!/usr/bin/env python3
"""费托合成（FTS）生成直链烃的基元反应网络生成器。

按规则枚举 C1 … C_nmax 的全部表面中间体和基元步骤，并对每一步做
C/H/O 原子守恒与表面位点守恒校验。

用法:
    python build_network.py --nmax 20 --out ../data

表面物种记号（n = 碳数）:
    烃类中间体 HC(n, a, b):  R'-CH_a-CH_b*   (R' = 正构烷基, C_{n-2})
        b = 端位 α-碳上的 H 数 (0..3), a = β-碳上的 H 数 (0..2; n=2 时 0..3)
    含氧中间体 OXY(frag, p, q):  frag-CH_p-OH_q*  (CO 插入产物及其加氢物)
    烯丙基 ALLYL(n):  R'-CH-CH-CH2*  (η3)
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

G = ["C", "CH", "CH2", "CH3"]  # CH_x 片段


# --------------------------------------------------------------------------
# 物种
# --------------------------------------------------------------------------
class Species:
    __slots__ = ("name", "C", "H", "O", "sites", "n", "kind")

    def __init__(self, name, C, H, O, sites, n, kind):
        self.name, self.C, self.H, self.O = name, C, H, O
        self.sites, self.n, self.kind = sites, n, kind

    def __repr__(self):
        return self.name


SPECIES: dict[str, Species] = {}


def sp(name, C, H, O, sites=1, n=None, kind=""):
    if name not in SPECIES:
        SPECIES[name] = Species(name, C, H, O, sites, C if n is None else n, kind)
    return SPECIES[name]


def gas(name, C, H, O, kind="气相"):
    return sp(f"{name}(g)", C, H, O, sites=0, kind=kind)


def tail(k):
    """CH3(CH2)_k 链尾（k >= 0）。"""
    if k == 0:
        return "CH3"
    if k == 1:
        return "CH3CH2"
    return f"CH3(CH2){k}"


HC_KIND = {
    (2, 3): "正构烷烃(物理吸附)",
    (2, 2): "1-烷基 alkyl",
    (2, 1): "亚烷基 alkylidene",
    (2, 0): "次烷基 alkylidyne",
    (1, 3): "2-烷基(仲烷基)",
    (1, 2): "α-烯烃(π/di-σ 吸附)",
    (1, 1): "1-烯基 alkenyl",
    (1, 0): "亚乙烯基型 vinylidene",
    (0, 3): "2-亚烷基",
    (0, 2): "1-烯-2-基",
    (0, 1): "1-炔(吸附)",
    (0, 0): "炔基 alkynyl",
}
C2_KIND = {
    (3, 3): "乙烷(物理吸附)", (3, 2): "乙基", (3, 1): "亚乙基", (3, 0): "次乙基",
    (2, 2): "乙烯(π/di-σ 吸附)", (2, 1): "乙烯基", (2, 0): "亚乙烯基",
    (1, 1): "乙炔(吸附)", (1, 0): "乙炔基", (0, 0): "C2 二聚体",
}
C1_KIND = ["表面碳(碳化物)", "次甲基 CH", "亚甲基 CH2", "甲基 CH3"]


def c1(z):
    return sp(f"{G[z]}*", 1, z, 0, kind=C1_KIND[z])


def hc(n, a, b):
    """烃类中间体 R'-CH_a-CH_b*。"""
    if n == 1:
        return c1(b)
    if n == 2:
        a, b = max(a, b), min(a, b)
        return sp(f"{G[a]}{G[b]}*", 2, a + b, 0, kind=C2_KIND[(a, b)])
    assert 0 <= a <= 2 and 0 <= b <= 3
    head = tail(n - 2) if a == 2 else tail(n - 3) + G[a]   # β-CH2 并入链尾书写
    return sp(f"{head}{G[b]}*", n, 3 + 2 * (n - 3) + a + b, 0,
              kind=HC_KIND[(a, b)])


def frags(m):
    """碳数 m 的可继续链增长的片段（端位碳仍有空价键）。返回 [(species, b)]。"""
    if m == 1:
        return [(c1(z), z) for z in range(4)]
    sat = 3 if m == 2 else 2
    return [(hc(m, sat, b), b) for b in range(3)]


OXY_G = {(0, 0): "CO", (0, 1): "COH", (1, 0): "CHO", (1, 1): "CHOH",
         (2, 0): "CH2O", (2, 1): "CH2OH"}
OXY_KIND = {(0, 0): "酰基/烯酮型(CO 插入产物)", (0, 1): "羟基卡宾/烯醇型",
            (1, 0): "醛型", (1, 1): "羟烷基型", (2, 0): "烷氧基型", (2, 1): "醇型"}


def oxy(frag, p, q):
    """frag-CH_p-OH_q*：链端带 C–O 的含氧中间体。"""
    return sp(f"{frag.name[:-1]}{OXY_G[(p, q)]}*", frag.C + 1, frag.H + p + q, 1,
              kind=OXY_KIND[(p, q)])


def allyl(n):
    name = "CH2CHCH2*" if n == 3 else f"{tail(n - 4)}CHCHCH2*"
    return sp(name, n, 2 * n - 1, 0, kind="η3-烯丙基")


def alkane_g(n):
    return gas("CH4" if n == 1 else f"C{n}H{2 * n + 2}", n, 2 * n + 2, 0, "产物:正构烷烃")


def alkene_g(n):
    return gas(f"C{n}H{2 * n}", n, 2 * n, 0, "产物:直链 α-烯烃")


# --------------------------------------------------------------------------
# 反应
# --------------------------------------------------------------------------
REACTIONS: list[dict] = []
_seen = set()


def rxn(n, cat, lhs, rhs, note=""):
    """登记一个可逆基元步骤；自动补空位 * 并校验守恒。lhs/rhs 为物种列表。"""
    lhs, rhs = list(lhs), list(rhs)
    for el in "CHO":
        if sum(getattr(s, el) for s in lhs) != sum(getattr(s, el) for s in rhs):
            raise ValueError(f"{el} 不守恒: {lhs} -> {rhs}")
    d = sum(s.sites for s in lhs) - sum(s.sites for s in rhs)
    stars_l, stars_r = max(-d, 0), max(d, 0)
    key = (frozenset(Counter(s.name for s in lhs).items()),
           frozenset(Counter(s.name for s in rhs).items()))
    if key in _seen or (key[1], key[0]) in _seen:
        return
    _seen.add(key)

    def side(specs, stars):
        terms = []
        for name, k in Counter(s.name for s in specs).items():
            terms.append(name if k == 1 else f"{k} {name}")
        if stars:
            terms.append("*" if stars == 1 else f"{stars} *")
        return " + ".join(terms)

    REACTIONS.append({"n": n, "category": cat,
                      "equation": f"{side(lhs, stars_l)} ⇌ {side(rhs, stars_r)}",
                      "note": note})


def build(nmax):
    # ---- 通用小分子 ----
    CO, H, O, OH = sp("CO*", 1, 0, 1, n=1), sp("H*", 0, 1, 0), sp("O*", 0, 0, 1), sp("OH*", 0, 1, 1)
    H2O, CO2 = sp("H2O*", 0, 2, 1), sp("CO2*", 1, 0, 2, n=1)
    COOH, HCOO = sp("COOH*", 1, 1, 2, n=1, kind="羧基"), sp("HCOO*", 1, 1, 2, n=1, kind="甲酸根")
    COg, H2g, H2Og, CO2g = gas("CO", 1, 0, 1), gas("H2", 0, 2, 0), gas("H2O", 0, 2, 1), gas("CO2", 1, 0, 2)
    C, CH, CH2, CH3 = (c1(z) for z in range(4))
    # C1 含氧中间体 CH_pO / CH_pOH
    HCO, COH = sp("HCO*", 1, 1, 1, kind="甲酰基"), sp("COH*", 1, 1, 1, kind="羟基次甲基")
    HCOH, CH2O = sp("HCOH*", 1, 2, 1, kind="羟基卡宾"), sp("CH2O*", 1, 2, 1, kind="甲醛型")
    CH2OH, CH3O = sp("CH2OH*", 1, 3, 1, kind="羟甲基"), sp("CH3O*", 1, 3, 1, kind="甲氧基")

    # 1. 反应物吸附
    cat = "1 反应物吸附"
    rxn(0, cat, [COg], [CO], "CO 分子吸附")
    rxn(0, cat, [H2g], [H, H], "H2 解离吸附")

    # 2. CO 活化
    cat = "2a CO 直接解离"
    rxn(1, cat, [CO], [C, O], "碳化物机理的起点")
    rxn(1, cat, [CO, CO], [C, CO2], "Boudouard 歧化")
    cat = "2b CO 氢助活化(加氢)"
    rxn(1, cat, [CO, H], [HCO])
    rxn(1, cat, [CO, H], [COH])
    rxn(1, cat, [HCO, H], [CH2O])
    rxn(1, cat, [HCO, H], [HCOH])
    rxn(1, cat, [COH, H], [HCOH])
    rxn(1, cat, [CH2O, H], [CH3O])
    rxn(1, cat, [CH2O, H], [CH2OH])
    rxn(1, cat, [HCOH, H], [CH2OH])
    cat = "2c 氢助 C–O 断裂"
    rxn(1, cat, [HCO], [CH, O])
    rxn(1, cat, [COH], [C, OH])
    rxn(1, cat, [HCOH], [CH, OH])
    rxn(1, cat, [CH2O], [CH2, O])
    rxn(1, cat, [CH2OH], [CH2, OH])
    rxn(1, cat, [CH3O], [CH3, O])

    # 3. 除氧 / 水煤气变换
    cat = "3a 除氧生成 H2O"
    rxn(0, cat, [O, H], [OH])
    rxn(0, cat, [OH, H], [H2O])
    rxn(0, cat, [OH, OH], [H2O, O], "羟基歧化")
    rxn(0, cat, [H2O], [H2Og], "H2O 脱附")
    cat = "3b 除氧生成 CO2 / 水煤气变换"
    rxn(1, cat, [CO, O], [CO2])
    rxn(1, cat, [CO, OH], [COOH], "羧基路径")
    rxn(1, cat, [COOH], [CO2, H])
    rxn(1, cat, [COOH, O], [CO2, OH])
    rxn(1, cat, [COOH, OH], [CO2, H2O])
    rxn(1, cat, [HCO, O], [HCOO], "甲酸根路径")
    rxn(1, cat, [HCOO], [CO2, H])
    rxn(1, cat, [CO2], [CO2g], "CO2 脱附")

    # 4. C1 逐步加氢与甲烷
    cat = "4 C1 加氢(链引发单体 CHx 生成)"
    rxn(1, cat, [C, H], [CH])
    rxn(1, cat, [CH, H], [CH2])
    rxn(1, cat, [CH2, H], [CH3], "链引发物种 CH3*")
    CH4 = sp("CH4*", 1, 4, 0, kind="甲烷(物理吸附)")
    rxn(1, cat, [CH3, H], [CH4], "甲烷化(链终止 n=1)")
    rxn(1, cat, [CH4], [alkane_g(1)], "CH4 脱附")

    # 5. 链增长 / 终止, 逐碳数
    for n in range(2, nmax + 1):
        m = n - 1
        # 5a CHx–CHy 偶联 (碳化物机理)
        cat = "5a 链增长: C–C 偶联(碳化物机理)"
        for f, b in frags(m):
            for z in range(4):
                if m == 1:
                    prod = hc(2, b, z)
                else:
                    prod = hc(n, b, z)
                note = ""
                if z == 2 and ((m == 1 and b == 3) or (m > 1 and b == 2)):
                    note = "烷基机理(Brady–Pettit): CH2 插入金属–烷基键"
                elif z == 2 and b == 1 and m > 1:
                    note = "亚烷基 + CH2 → 吸附 α-烯烃"
                elif z == 1 and b == 1 and m > 1:
                    note = "亚烷基 + CH → 烯基"
                elif z == 1 and b == 0 and m > 1:
                    note = "次烷基 + CH"
                elif z == 0:
                    note = "与表面 C 偶联"
                rxn(n, cat, [f, c1(z)], [prod], note)

        # 5b CO 插入机理 (Pichler–Schulz): 片段(m) + CO → 含氧中间体 → C–O 断裂得 C_n
        cat = "5b 链增长: CO 插入"
        for f, b in frags(m):
            rxn(n, cat, [f, CO], [oxy(f, 0, 0)],
                "CO 插入金属–烷基键(Pichler–Schulz)" if f.kind in ("1-烷基 alkyl", "乙基", "甲基 CH3") else "")
        cat = "5c CO 插入中间体加氢"
        for f, b in frags(m):
            for q in (0, 1):
                for p in (0, 1):
                    rxn(n, cat, [oxy(f, p, q), H], [oxy(f, p + 1, q)], "C 端加氢")
            for p in (0, 1, 2):
                rxn(n, cat, [oxy(f, p, 0), H], [oxy(f, p, 1)], "O 端加氢")
        cat = "5d CO 插入中间体 C–O 断裂"
        for f, b in frags(m):
            for p in (0, 1, 2):
                prod = hc(2, b, p) if m == 1 else hc(n, b, p)
                rxn(n, cat, [oxy(f, p, 0)], [prod, O])
                rxn(n, cat, [oxy(f, p, 1)], [prod, OH])

        # 5e 烯醇(羟基卡宾)缩合机理 (Storch–Anderson)
        cat = "5e 链增长: 烯醇缩合(Storch–Anderson)"
        if m == 1:
            rxn(n, cat, [HCOH, HCOH], [oxy(CH, 0, 1), H2O], "两个羟基卡宾脱水缩合")
        else:
            sat_prev = frags(m - 1)[-1][0]          # C_{m-1} 的饱和端片段(甲基/烷基)
            enol = oxy(sat_prev, 0, 1)              # R-C(OH)*, 共 m 个碳
            alkylidyne = frags(m)[0][0]             # R-C*, m 个碳
            rxn(n, cat, [enol, HCOH], [oxy(alkylidyne, 0, 1), H2O], "RC(OH)* + HCOH* 脱水缩合")

        # 5f 烯基机理 (Maitlis): 烯基 + CH2 → 烯丙基 → 异构为高一级烯基
        if n >= 3:
            cat = "5f 链增长: 烯基机理(Maitlis)"
            alkenyl_prev = hc(2, 2, 1) if m == 2 else hc(m, 1, 1)
            rxn(n, cat, [alkenyl_prev, CH2], [allyl(n)], "烯基 + CH2 → 烯丙基")
            rxn(n, cat, [allyl(n)], [hc(n, 1, 1)], "烯丙基 ⇌ 烯基 (1,3-H 迁移)")
            rxn(n, cat, [allyl(n), H], [hc(n, 1, 2)], "烯丙基加氢 → 吸附 α-烯烃")

        # 5g C_n 中间体的加氢/脱氢
        cat = "5g C_n 中间体加氢/脱氢"
        if n == 2:
            for a in range(4):
                for b in range(a + 1):
                    for (a2, b2) in ((a + 1, b), (a, b + 1)):
                        if a2 <= 3 and b2 <= 3:
                            note = ""
                            if (a, b) == (3, 2):
                                note = "链终止: 乙基加氢 → 乙烷"
                            elif (a, b) == (2, 2):
                                note = "逆反应 = 乙基 β-H 消除 → 乙烯"
                            rxn(n, cat, [hc(2, a, b), H], [hc(2, a2, b2)], note)
        else:
            for a in range(3):
                for b in range(4):
                    if b < 3:
                        note = "链终止: 烷基加氢(还原消除) → 正构烷烃" if (a, b) == (2, 2) else "α-C 加氢"
                        rxn(n, cat, [hc(n, a, b), H], [hc(n, a, b + 1)], note)
                    if a < 2:
                        note = "逆反应 = 烷基 β-H 消除 → α-烯烃(链终止)" if (a, b) == (1, 2) else "β-C 加氢"
                        rxn(n, cat, [hc(n, a, b), H], [hc(n, a + 1, b)], note)

        # 5h 产物脱附 (逆反应 = 烯烃再吸附 → 二次反应)
        cat = "5h 产物脱附/再吸附"
        if n == 2:
            rxn(n, cat, [hc(2, 3, 3)], [alkane_g(2)], "乙烷脱附")
            rxn(n, cat, [hc(2, 2, 2)], [alkene_g(2)], "乙烯脱附; 逆向 = 再吸附")
        else:
            rxn(n, cat, [hc(n, 2, 3)], [alkane_g(n)], "正构烷烃脱附")
            rxn(n, cat, [hc(n, 1, 2)], [alkene_g(n)], "α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发)")


SMALL_FV = {"H*": 1, "O*": 2, "OH*": 1, "H2O*": 0, "CO*": 0, "CO2*": 0, "COOH*": 1, "HCOO*": 1}
SMALL_RAD = {"H*": "氢原子 H·", "O*": "氧原子 ·O·", "OH*": "羟基自由基 ·OH", "COOH*": "羧基自由基 ·COOH",
             "HCOO*": "甲酰氧基自由基 HCOO·"}
RAD = {
    "表面碳(碳化物)": "碳原子", "次甲基 CH": "次甲基自由基 ·CH (卡拜)", "亚甲基 CH2": "亚甲基 :CH2 (卡宾)",
    "甲基 CH3": "甲基自由基 ·CH3", "甲酰基": "甲酰基自由基 HĊ=O", "羟基次甲基": "羟基卡拜 ·C–OH",
    "羟基卡宾": "羟基卡宾 H–C̈–OH", "甲醛型": "甲醛 CH2=O (闭壳层)", "羟甲基": "羟甲基自由基 ·CH2OH",
    "甲氧基": "甲氧基自由基 CH3O·",
    "1-烷基 alkyl": "伯烷基自由基 RCH2ĊH2", "乙基": "乙基自由基 CH3ĊH2",
    "亚烷基 alkylidene": "烷基卡宾 RCH2C̈H", "亚乙基": "甲基卡宾 CH3C̈H",
    "次烷基 alkylidyne": "烷基卡拜 RCH2C⋮", "次乙基": "甲基卡拜 CH3C⋮",
    "2-烷基(仲烷基)": "仲烷基自由基 RĊHCH3", "α-烯烃(π/di-σ 吸附)": "α-烯烃 (闭壳层) / 1,2-双自由基",
    "乙烯(π/di-σ 吸附)": "乙烯 (闭壳层) / 1,2-双自由基",
    "1-烯基 alkenyl": "烯基自由基 RCH=ĊH", "乙烯基": "乙烯基自由基 CH2=ĊH",
    "亚乙烯基型 vinylidene": "亚乙烯基卡宾 RCH=C:", "亚乙烯基": "亚乙烯基卡宾 CH2=C:",
    "2-亚烷基": "二烷基卡宾 R–C̈–CH3", "1-烯-2-基": "2-烯基自由基 RĊ=CH2",
    "1-炔(吸附)": "1-炔 (闭壳层)", "乙炔(吸附)": "乙炔 (闭壳层)",
    "炔基 alkynyl": "炔基自由基 RC≡C·", "乙炔基": "乙炔基自由基 HC≡C·", "C2 二聚体": "C2 双原子碳",
    "η3-烯丙基": "烯丙基自由基 (共振稳定)",
    "酰基/烯酮型(CO 插入产物)": "含氧端基: 酰基 –Ċ=O", "羟基卡宾/烯醇型": "含氧端基: 羟基卡宾 –C̈–OH",
    "醛型": "含氧端基: 醛基 –CH=O", "羟烷基型": "含氧端基: α-羟烷基 –ĊHOH",
    "烷氧基型": "含氧端基: 烷氧基 –CH2O·", "醇型": "含氧端基: 醇 –CH2OH",
}


def free_valence(s):
    if s.sites == 0:
        return 0
    if s.name in SMALL_FV:
        return SMALL_FV[s.name]
    return 2 * s.C + 2 - s.H


def radical_form(s):
    if s.sites == 0:
        return "闭壳层分子"
    if s.name in SMALL_RAD:
        return SMALL_RAD[s.name]
    if free_valence(s) == 0:
        return "闭壳层分子(物理吸附)"
    r = RAD.get(s.kind, "")
    if r.startswith("含氧端基") and s.H == 2 * s.C + 1:
        r += " (单自由基)"
    return r


# --------------------------------------------------------------------------
# 输出
# --------------------------------------------------------------------------
def write(out: Path, nmax: int):
    out.mkdir(parents=True, exist_ok=True)
    for i, r in enumerate(REACTIONS, 1):
        r["id"] = f"R{i:04d}"
    with open(out / "reactions.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["id", "n", "category", "equation", "note"])
        w.writeheader()
        w.writerows(REACTIONS)
    order = sorted(SPECIES.values(), key=lambda s: (s.n, s.sites == 0, s.O, s.H, s.name))
    with open(out / "species.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["species", "carbon_number", "C", "H", "O", "phase", "type",
                    "free_valence", "radical_form"])
        for s in order:
            w.writerow([s.name, s.n, s.C, s.H, s.O, "gas" if s.sites == 0 else "surface", s.kind,
                        free_valence(s), radical_form(s)])
    sl = [f"# 物种全表 (C1–C{nmax})：{len(SPECIES)} 个", "",
          "自由价 = 该片段脱离表面后未饱和的价键总数（= 2C+2−H，含 C–C 多重键折算；"
          "小分子按常规价态）。自由价 1 = 单自由基，2 = 双自由基/卡宾，3 = 卡拜型，依此类推；"
          "0 = 闭壳层分子。", ""]
    cur = None
    for s in order:
        if s.n != cur:
            cur = s.n
            sl += ["", f"## {'通用小分子' if cur == 0 else 'C' + str(cur)}", "",
                   "| 物种 | 分子式 | 相 | 类型 | 自由价 | 游离态对应物(自由基) |", "|---|---|---|---|---|---|"]
        fm = (f"C{s.C}" if s.C > 1 else "C" * s.C) + (f"H{s.H}" if s.H > 1 else "H" * s.H) \
            + (f"O{s.O}" if s.O > 1 else "O" * s.O)
        sl.append(f"| `{s.name}` | {fm} | {'气相' if s.sites == 0 else '表面'} | {s.kind} | "
                  f"{free_valence(s)} | {radical_form(s)} |")
    (out / "species_full.md").write_text("\n".join(sl) + "\n", encoding="utf-8")

    lines = [f"# 费托合成直链烃基元反应网络全表 (C1–C{nmax})", "",
             f"由 `scripts/build_network.py --nmax {nmax}` 自动生成；"
             f"共 {len(SPECIES)} 个物种、{len(REACTIONS)} 个可逆基元步骤，"
             "每一步均通过 C/H/O 原子守恒与表面位点守恒校验。`*` 为表面空位。", ""]
    cur_n, cur_cat = None, None
    for r in REACTIONS:
        if r["n"] != cur_n:
            cur_n, cur_cat = r["n"], None
            title = {0: "通用步骤（与碳数无关）", 1: "C1 化学"}.get(cur_n, f"C{cur_n - 1} → C{cur_n}")
            lines += ["", f"## {title}"]
        if r["category"] != cur_cat:
            cur_cat = r["category"]
            lines += ["", f"### {cur_cat}", "", "| 编号 | 基元步骤 | 说明 |", "|---|---|---|"]
        lines.append(f"| {r['id']} | `{r['equation']}` | {r['note']} |")
    (out / "reactions_full.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nmax", type=int, default=20, help="最大碳数 (默认 20)")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "data")
    args = ap.parse_args()
    build(args.nmax)
    write(args.out, args.nmax)
    cats = Counter(r["category"] for r in REACTIONS)
    print(f"nmax={args.nmax}: {len(SPECIES)} 物种, {len(REACTIONS)} 基元步骤")
    for c in sorted(cats):
        print(f"  {c}: {cats[c]}")


if __name__ == "__main__":
    main()
