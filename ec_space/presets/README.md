# 不同溶剂体系怎么跑：EC、DMC、EA、DME，以及 EC/DMC 混合

## 1. 先弄清楚：换溶剂，哪些变、哪些不变

| 溶剂 | SMILES | 分子式 | 重原子 | 两个分子的重原子 | 两个分子的原子包络 |
|---|---|---|---|---|---|
| EC 碳酸乙烯酯 | `O=C1OCCO1` | C₃H₄O₃ | 6 | 12 | C₆H₈O₆ |
| DMC 碳酸二甲酯 | `COC(=O)OC` | C₃H₆O₃ | 6 | 12 | C₆H₁₂O₆ |
| EA 乙酸乙酯 | `CCOC(C)=O` | C₄H₈O₂ | 6 | 12 | C₈H₁₆O₄ |
| DME 1,2-二甲氧基乙烷 | `COCCOC` | C₄H₁₀O₂ | 6 | 12 | C₈H₂₀O₄ |
| EC + DMC 混合 | 两者都有 | — | 6 + 6 | 12 | EC+EC、EC+DMC（C₆H₁₀O₆）、DMC+DMC 的并集：C ≤ 6，O ≤ 6，H ≤ 12 |

DME 在电池里指 1,2-二甲氧基乙烷（glyme）。如果你要的是 CO₂ 化学里的二甲醚 `COC`，在预设里写 `solvents = DIMETHYL_ETHER`。

这四种溶剂恰好都是 **6 个重原子**，所以"最多两个溶剂分子那么大"都是 **≤ 12 个重原子**。

**不变的部分**：硬枚举的定义是"C/H/O 能写出的全部结构"，所以下面这些对所有溶剂**完全相同**：
- 闭壳层分子；
- 全部自由基、离子、Li⁺ 配位物；
- 所有单分子反应和与 e⁻ / Li⁺ / H 的反应。

EC 的空间里本来就包含 DMC、EA、DME，以及它们的全部碎片。

**随溶剂变化的有三点**：

1. **与溶剂分子的双分子反应**（新增的 `solvent` 阶段，只有这一步需要按溶剂重算）。对空间里每一个定域、不含 Li 的自由基 / 负离子 X，和每一种溶剂分子 S：

   | 反应族 | 写法 | 例子 |
   |---|---|---|
   | 从溶剂抽 H | X• + S ⇌ XH + S• | •CH₃ + DME → CH₄ + CH₃OĊHCH₂OCH₃ |
   | 夺溶剂质子 | X⁻ + S ⇌ XH + S⁻ | CH₃O⁻ + EA → CH₃OH + ⁻CH₂COOEt（烯醇负离子） |
   | 亲核酰基取代 | X⁻ + S(C(=O)–OR) ⇌ X–C(=O)… + RO⁻ | EtO⁻ + DMC → EMC + CH₃O⁻（酯交换）；CH₃O⁻ + EC → CH₃OCO₂CH₂CH₂O⁻（开环）；烯醇负离子 + EA → 乙酰乙酸乙酯（Claisen） |
   | SN2 烷基转移 | X⁻ + S(R–OR′) ⇌ X–R + R′O⁻ | RO⁻ + DMC → ROCH₃ + CH₃OCO₂⁻（DMC 甲基化）；RO⁻ + DME → ROCH₃ + ⁻OCH₂CH₂OCH₃（醚断裂）；RO⁻ + EC → ROCH₂CH₂OCO₂⁻ |

   - EC/DMC 混合时，两种溶剂都作为伙伴。所以 EC 与 DMC 之间互相抽 H、互相酯交换（例如 CH₃O⁻ 进攻 EC 开环，得到 DMDOHC 前体），都自动包括在内。
   - 产物超过 `max_heavy`、或不属于配置的物种类别的反应，不写出。
   - 两边都能生成的互逆反应只保留一条；退化交换（例如 CH₃O⁻ + DMC → DMC + CH₃O⁻）去掉。

2. **起点**：网页和路径搜索从哪个（哪几个）溶剂分子出发。EC/DMC 预设里两个都是起点，路径会从更近的那个算起。

3. **（可选）组成包络**：见下面的方案 B。

## 2. 方案 A（推荐）：共用一个基础空间，每种溶剂只补一个阶段

所有预设都把 `run.out` 指向同一个目录 `ec_space_out/`。

第一次：完整运行，约 1 400 核时，和原来一样：

```bash
bash run.sh presets/EC.ini
```

之后每加一种溶剂，只跑 solvent、stats、viewer 三个阶段。基础数据直接复用，不再重算：

```bash
bash run.sh presets/DMC.ini    --stages solvent,stats,viewer
bash run.sh presets/EA.ini     --stages solvent,stats,viewer
bash run.sh presets/DME.ini    --stages solvent,stats,viewer
bash run.sh presets/EC_DMC.ini --stages solvent,stats,viewer
```

在提交文件里也只改这一行。Windows / macOS 用 `python run.py presets/xxx.ini ...`，参数相同。

每种溶剂单独存放、互不覆盖：

| 内容 | 路径 |
|---|---|
| 溶剂反应 | `ec_space_out/solvent/<TAG>/reactions/L{n}/part-*.tsv.gz` |
| 新物种（产物大于 `derive_max_heavy` 时） | `ec_space_out/solvent/<TAG>/species/L{n}/` |
| 溶剂反应位点、溶剂说明 | `SOLVENTS.json`（在 `solvent/<TAG>/` 下） |
| 汇总 | `ec_space_out/SUMMARY_<TAG>.md`、`stats_<TAG>.json` |
| 网页 | `ec_space_out/viewer_<TAG>/index.html` |

`<TAG>` 为 `EC`、`DMC`、`EA`、`DME`、`EC+DMC`。

查询时可以只看某个溶剂体系：

```bash
python -m ecspace.query ec_space_out --solvent EC+DMC --reactions-of "COC(=O)OC"
```

## 3. 方案 B：组成包络，只要"原子最多来自两个溶剂分子"的物种

如果只关心某个溶剂体系、不需要开放的 C/H/O 空间：
1. 在对应预设里取消 `max_C / max_O / max_H` 三行的注释（数值已按上表填好）；
2. 把 `run.out` 改成该溶剂自己的目录；
3. 完整运行一次。

需要知道的取舍：
- 空间小很多，所以更快；但不能和别的溶剂共用基础数据。
- `max_H` 只用于派生阶段的筛选（过滤 H 数超标的母体），枚举阶段仍会生成全部 H 数。
- EC/DMC 混合用并集包络（C ≤ 6，O ≤ 6，H ≤ 12）。它同时覆盖 EC+EC、EC+DMC、DMC+DMC。

## 4. 没有包括的（如需要可以再加）

- 自由基加成到溶剂的 C=O / C=C、碳正离子与溶剂的反应（溶剂作为亲核体进攻 R⁺）。
- 以 Li⁺ 配位物、离域的 M•⁻ / M•⁺ 为反应物的溶剂反应。它们与 Li⁺、e⁻ 之间的转化已经在基础网络里，所以拓扑上仍然连通。
- 任意两个非溶剂大物种之间的双分子反应（规模是平方级）。
- F、P（LiPF₆ / FEC）。不在 C/H/O 空间内。
