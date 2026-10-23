# 文献 DFT 数值提取说明（供本地 Claude 按此执行）

目的：为 Ni/CeO₂ 干重整积炭微动力学模型的 20 个基元步和 14 个表面物种，
从文献中系统提取 DFT 能垒、反应能和吸附能，带完整出处，用于
(1) 生成每个参数的文献区间（替换 `../energy_ranges_literature.json` 中标为 LIT 的估计值），
(2) 生成每个催化剂体系的过渡态能量表（用于控制相图）。

**原则**
- 过渡态能量与物种能量一律独立记录，不用 BEP 推算，没有算过的就留空。
- 每条记录只对应一篇论文里的一个数字，不做任何换算以外的加工；单位统一为 eV（kJ/mol ÷ 96.485，kcal/mol ÷ 23.06）。
- 必须记录计算设置：泛函、U 值、色散校正、覆盖度、slab 模型、k 点可省。
- 能垒注明定义：相对直接反应物吸附态（IS）还是相对气相参考；反应能注明终态定义。
- 不确定的字段写 `?`，不要猜。

## 1. 提取模板

文件 `extraction_template.csv`，一行一个数字，列定义：

| 列 | 含义 | 示例 |
|---|---|---|
| doi | 论文 DOI | 10.1021/acscatal.5c07743 |
| year | 年份 | 2026 |
| first_author | 第一作者 | Naveen |
| location | 数字所在位置 | Table S2 / Fig. 3b / text p.5 |
| system | 催化剂体系 | Ni(111) / Ni(211) / Ni4-CeO2(111) / CeO2(111) / Ni13 cluster |
| site_type | 位点类型 | terrace / step / interface / support / subsurface |
| functional | 泛函 | PBE / RPBE / BEEF-vdW / PBE+U |
| U_eV | U 值（无则空） | 5.0 (Ce 4f) |
| dispersion | 色散校正 | none / D3 / vdW-DF |
| coverage_ML | 覆盖度 | 1/9 / 1/4 |
| quantity_type | 数字类型 | TS_barrier / reaction_energy / adsorption_energy / TS_absolute |
| reaction_id | 对应反应编号 | R2 / R9 / R18 / R22 |
| species | 对应物种（吸附能时填） | C* / OHv / Ov |
| reaction_text | 论文原文写法 | CH4* → CH3* + H* |
| value_eV | 数值，eV | 0.98 |
| reference_state | 参考态 | IS adsorbed / gas CH4 / clean slab + 1/2 H2 |
| direction | 正向或逆向 | forward / reverse |
| temperature_K | 若为自由能，温度 | 973 |
| is_free_energy | 是否自由能 | 0 / 1 |
| notes | 备注 | ZPE included; 4-layer slab |

反应编号对照（与 `../mkmc/network.py` 一致）：
R1 CH4 吸附；R2 CH4*→CH3*+H*；R3 CH3*→CH2*+H*；R4 CH2*→CH*+H*；R5 CH*→C*+H*；
R6 2H*→H2；R7 CH*+O*→CHO*；R8 CHO*→CO*+H*；R9 C*+O*→CO*；R10 CO 脱附；
R11 Ov+*→O*+v（界面氧传输）；R12 CO2+v→CO2v；R13 CO2v→CO+Ov；R14 CO2+Ov→CO3v；
R15 CO3v+*→CO2v+O*；R16 H*+Ov→OHv（氢溢流）；R17 OHv+H*→H2O；
R18 C*+C*→C2*（C–C 偶联或石墨成核）；R19 C2*+O*→CO*+C*；R20 C2*→沉积碳；
R21 Cg*+O*→CO*（石墨碳气化）；R22 C*→C(次表面/体相)；R23 体相碳析出；R24 碳岛边缘生长。

## 2. 检索关键词（按步骤分组）

- R1–R5：`CH4 dissociation Ni(111) DFT barrier`、`methane activation Ni(211) step DFT`、`CH4 dehydrogenation Ni/CeO2 interface DFT`、`CHx dehydrogenation nickel first principles`
- R6：`hydrogen desorption Ni(111) DFT barrier`、`H2 dissociative adsorption nickel`
- R7–R9：`carbon oxidation C+O Ni(111) DFT`、`CHO formation Ni DFT`、`CO formation coke removal Ni O* OH* DFT`
- R10：`CO adsorption energy Ni(111) DFT overbinding`
- R11、R15：`oxygen transport Ni/CeO2 interface DFT`、`lattice oxygen spillover ceria nickel DFT`、`oxygen vacancy formation Ni/CeO2(111) interface DFT+U`
- R12–R14：`CO2 dissociation CeO2(111) oxygen vacancy DFT+U`、`carbonate formation ceria DFT`、`CO2 activation reduced ceria`
- R16、R17：`hydrogen spillover Ni/CeO2 DFT`、`hydroxyl formation CeO2(111) DFT+U`、`water formation ceria hydroxyl DFT`
- R18–R21：`carbon dimer formation Ni(111) DFT`、`C-C coupling nickel coke DFT`、`graphene nucleation Ni step edge DFT`、`carbon island nucleation Ni(111) Gao Ding`、`coke gasification graphene edge oxygen Ni DFT`
- R22–R24：`subsurface carbon Ni(111) DFT`、`carbon diffusion nickel bulk first principles`、`carbon solubility nickel DFT`、`carbon nanotube growth Ni DFT mechanism`
- 物种能量：`adsorption energy C H O CO Ni(111) DFT`、`OH adsorption CeO2(111) DFT+U`、`oxygen vacancy formation energy CeO2(111) U dependence`

种子文献：Naveen 等 ACS Catal. 2026, 16, 2535 的 SI 表 S1–S3（已汇总 R2、R6、R13 多组值）及其参考文献 10–17、26、37、38、45–51。

## 3. 优先级

按敏感性结果，先做：R2、R13、R9、R6、R18、R11 的过渡态；OHv、C*、Ov、O* 的能量；R22 次表面碳能量与溶解度。其余步骤各找 3–5 个值即可。

## 4. 聚合

抽取完成后运行 `python aggregate_ranges.py`（待写），它会：
- 按 reaction_id / species 分组，给出 10%、50%、90% 分位数作为区间，并列出离群值；
- 按 system 分组生成每个催化剂体系的过渡态表 `ts_by_system.csv`；
- 写出新的 `../energy_ranges_literature.json`，每个参数附记录数和出处列表。
