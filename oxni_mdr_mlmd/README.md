# OxNi 甲烷干重整：通用 ML 势 MD 模型与工作流

对象：He et al., ACS Catal. 2021, 11, 12409。约 15 nm Ni 颗粒外覆 1–2 nm 多元氧化物（MEO：Ni²⁺/Mg/Al/Si/O），
用于 CH₄ + CO₂ 干重整（MDR），873 K 下几乎不积碳。
**机理研究和动态模拟的设计见 [`MECHANISM_PLAN.md`](MECHANISM_PLAN.md)；势函数和工作流见 [`PLAN.md`](PLAN.md)。** 本文件只说明代码和数据。

## 目录

| 路径 | 内容 |
|---|---|
| `scripts/meo_lib.py` | 建模库：Ni(211)/(111)/(100) slab、随机岩盐 MEO 膜（按电中性配置阳离子/O 空位）、开孔、CH₄/CO₂ 气体、C*/O*/H* 吸附原子、球状 Ni 颗粒 + 多孔 MEO 壳 + Mg(Al)O 载体 |
| `scripts/build_models.py` | 生成 S0–S3、P1–P2 全部模型 → `structures/*.extxyz`，渲染图 → `renders/` |
| `scripts/run_md.py` | 分阶段 MD（opt / Langevin / Bussi / NVE，可线性升降温），按 tag 冻结，气体上方反射墙，在线统计物种、积碳和 Ni 结构指标，可选 PLUMED |
| `scripts/calculators.py` | 统一接口：MACE（含本地微调模型）、UMA、SevenNet、MatterSim、ORB、DeePMD/DPA、EMT |
| `scripts/species.py` | 图论物种识别（CH₄(g)、CH₃*、CO*、CO₃*…）、C–C 团簇、溶解 C、晶格 O 进入含碳物种（MvK 标记）、Ni 剥离/氧化/析出 |
| `scripts/add_species.py` | 给退火后的结构加气体或吸附原子 |
| `scripts/select_frames.py` | 主动学习：委员会力偏差选帧，供 DFT 标注 |
| `scripts/benchmark_umlip.py` | 通用势横评：Ni(111) 吸附能、反应能，及与 DFT 帧的 parity |
| `scripts/make_plumed.py`, `fes_reweight.py` | OPES 输入模板（一个 CH₄ 的 C–H 断裂、一个 CO₂ 的解离、C+O vs C+C 二维）与 FES/TST 速率。CH₄/CO₂ 模式要用不带气体的结构：脚本把一个分子放进开口轴线上的参考圆柱，并写出 `*_start.extxyz` 作为 MD 的起始结构 |
| `scripts/plot_runs.py` | 对比多条轨迹（MEO vs 裸 Ni） |
| `protocols/*.json` | 冒烟、MEO 退火、873/1073/1400 K 反应、积碳、烧结、OPES |
| `dft_templates/` | VASP 标注/AIMD INCAR、KPOINTS、vasprun→extxyz、MACE 微调脚本（默认基础模型是不加 U 的 `mace-matpes-pbe-0`） |
| `plumed/` | 由 `make_plumed.py` 生成的 S1/S2 输入示例，以及 CH₄/CO₂ 输入对应的起始结构 |
| `data/benchmark_umlip.json` | MACE 基础模型基准结果（无 DFT 对照，见下） |

Tag 约定：0 Ni，1 Ni（冻结），2 MEO，3 载体，4 载体（冻结），5 气体分子，6 预置吸附原子。

## 模型（`structures/SUMMARY.txt` 有完整组成）

| 文件 | 原子数 | 说明 |
|---|---|---|
| `S0_smoke_*` | 627 / 336 | 小模型，CPU 测试 |
| `S1_ni211_meo_pore_{nogas,gas}` / `S1_ref_ni211_bare_*` | 1271, 1367 / 660, 756 | Ni(211) 24.4×27.4 Å + 1.05 nm MEO + 半径 6 Å 的孔；12 CH₄ + 12 CO₂ |
| `S2_ni211_meo_pore_CO` / `S2_ref_ni211_bare_CO` | 1279 / 668 | 孔内（或同面积圆盘内）预置 4 C* + 4 O*，用于积碳竞争 |
| `S3_ni111_meo_pore_gas` / `S3_ref_*` | 1471 / 756 | Ni(111) 平台版本 |
| `P1_np3nm_meo_on_MgAlO` / `P1_ref_*` | 9258 / 5496 | 3 nm Ni + 1 nm 多孔 MEO 壳，置于 Mg(Al)O(100) 上 |
| `P2_pair_np2.5nm_meo` / `P2_ref_*` | 12952 / 7386 | 两个 2.5 nm 颗粒，用于研究合并 |

MEO 默认组成为等摩尔 Mg:Ni:Al:Si（论文没有给出比例，这只是工作假设），岩盐晶格 a = 4.20 Å，形式电荷为零。
初始结构是理想的外延膜，**生产前先用 `protocols/meo_anneal.json` 退火**。

## 已验证的内容（本容器：4 核 CPU，无 GPU，无 PLUMED）

- `build_models.py`：全部模型生成并渲染，MEO/载体形式电荷为 0。
- `run_md.py` + MACE-MP-0 small：S0 模型 FIRE 预优化 + 873 K Langevin 100 步，运行正常，1.8 s/步（0.024 ns/天，CPU）。
- `species.py`、`plot_runs.py`、`add_species.py`、`select_frames.py`（small+medium 两模型委员会）均已在冒烟输出上跑通（输出未保留，按用户要求目前只交付建模部分）。
- `benchmark_umlip.py`：已对 MACE-MP-0 small/medium 运行，结果见下（MACE-OMAT-0 等其余模型留到 GPU 环境再跑）。
- `make_plumed.py`：四个 CH₄/CO₂ 示例的起始结构已用 numpy 复算，分子在参考圆柱之内（墙的自变量在范围内，与其余原子的最小距离大于 4.7 Å，集合变量取反应物的值）。PLUMED 本身没有运行。
- **未在本容器测试**：UMA/SevenNet/MatterSim/ORB/DPA 的调用分支、PLUMED 输入的实际运行（没有 PLUMED kernel）、`fes_reweight.py`、VASP 模板、`finetune_mace.sh`（只核对了所用参数在 mace-torch 0.3.17 里存在）。首次使用时先用 `plumed driver` 对一段轨迹检查 CV。

## 基础模型基准（Ni(111) p(3×3)，单位 eV）

吸附能的参考是 CH₄、H₂、H₂O、CO、CO₂ 气相分子（定义见 `benchmark_umlip.py` 文件头）。这里**没有 DFT 对照**，只用来说明通用势在这个化学体系上的可靠程度。

| 量 | MACE small | MACE medium |
|---|---:|---:|
| H\* | -0.73 | -0.62 |
| C\* | 2.81 | 3.15 |
| O\* | 6.01 | 5.06 |
| OH\* | 3.84 | 3.37 |
| CH\* | 1.28 | 1.50 |
| CH2\* | 0.65 | 0.71 |
| CH3\* | -0.42 | -0.40 |
| CO\* | -0.06 | -0.51 |
| CH4(g)->CH3\*+H\* | -1.15 | -1.02 |
| CH3\*->CH2\*+H\* | 0.33 | 0.50 |
| CH2\*->CH\*+H\* | -0.09 | 0.18 |
| CH\*->C\*+H\* | 0.79 | 1.03 |
| CO2(g)->CO\*+O\* | 6.90 | 5.36 |
| CO\*->C\*+O\* | 5.12 | 4.55 |
| C\*+O\*->CO(g) | -5.05 | -4.04 |
| CO2(g)+H2(g)->CO(g)+H2O(g) | 0.95 | 0.81 |
| CH4(g)+CO2(g)->2CO(g)+2H2(g) | 4.72 | 4.98 |

如何解读：以 H₂O−H₂ 为参考的 O* 吸附能物理上应在 0 附近略负（O/Ni(111) 相对 ½O₂ 约 −2.5 到 −3 eV，H₂O 生成能约 −2.5 eV）。
MP 系模型给出 +5 到 +6 eV，即认为 O 几乎不在 Ni 上吸附。这与"含 O 的 Ni 体系在训练集中按 GGA+U 计算"的能量标度问题一致。
干重整的电子能反应能约 +3.2 eV（实验 ΔH = 2.56 eV 加零点能修正）。
所以：**先用统一的、不加 U 的 DFT 数据微调，再做任何反应能量学**（PLAN.md §2）。

## 运行

依赖：`ase numpy scipy matplotlib torch mace-torch`；生产还需 GPU、PLUMED（含 OPES）、DFT 软件。完整命令见 PLAN.md §7。
