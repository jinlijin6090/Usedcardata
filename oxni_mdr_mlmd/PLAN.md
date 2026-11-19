# OxNi 甲烷干重整：用通用机器学习势做"类 AIMD"反应动力学研究计划

参考文献：L. He, …, W. Shen, A.-H. Lu, *Robust and Coke-free Ni Catalyst Stabilized by 1–2 nm-Thick
Multielement Oxide for Methane Dry Reforming*, ACS Catal. 2021, 11, 12409–12416.

---

## 0. 从论文里提炼出的、值得用模拟回答的问题

| 论文中的观察 | 关键数据 | 对应的模拟问题 |
|---|---|---|
| ~15 nm 球状 Ni，表面以高指数面为主，外覆 1–2 nm MEO（Ni²⁺、Mg、Al、Si、O；类岩盐，条纹间距 2.2 Å） | XPS：MEO 中 Ni²⁺ 854.4 eV，界面 Ni²⁺ 856.0 eV | **Q5** MEO 在 873–1173 K 下是否保持混合、多孔？Ni²⁺ 会不会被还原析出（熵稳定是否成立）？ |
| 873 K 是最严重的积碳温度，但 OxNi 300 h 几乎不积碳 | IGA：OxNi 消碳 0.51 > 积碳 0.24 g g<sub>Ni</sub>⁻¹ h⁻¹；参比 Ni 积碳 2.39 > 消碳 1.79 | **Q1** 为什么 MEO 限域的 Ni 上 C–C 偶联（成焦）比 C(或 CHₓ)+O 气化慢？ |
| CH₄ 活化是速控步，台阶比平台快 | 15 nm 兼顾活性与抗氧化 | **Q2** CH₄ 首个 C–H 断裂在哪里最快：暴露的台阶/高指数 Ni、Ni/MEO 周界、还是 MEO 中的 Ni²⁺？ |
| 外层碱性位促进 CO₂ 吸附与消碳 | CO₂-TPD；IGA 交替通 CH₄/CO₂ | **Q3** CO₂ 在 Ni、周界、Mg–O 碱性位上的活化路径；MEO 晶格 O 是否参与（Mars–van Krevelen）？ |
| 1073–1173 K 原位 TEM 中粒径/形貌不变；参比 Ni 在 CH₄ 下迅速长成大的有棱角颗粒 | XRD 晶粒 15 nm 不变 | **Q4** MEO 壳如何抑制 Ni 原子剥离、表面扩散和颗粒合并；表面 C/H 覆盖如何驱动重构？ |

主线：**Q1（积碳–消碳平衡）** 最能体现"动态过程"，也最直接对应论文的核心结论，建议作为主要成果；Q2–Q4 作为机理支撑，Q5 作为结构基础。

---

## 1. 先说清楚：MD 能直接看到什么，看不到什么

"用通用 ML 势做 AIMD" 的实际含义是：**ML 势 MD（MLMD）做生产，真正的 DFT-AIMD 只做验证和产生训练数据**。DFT-AIMD 只能做 ~200–300 原子、~10 ps；MLMD 可以做 10³–10⁴ 原子、ns 量级。但 ns 依然很短：

单个位点上一个元反应的平均等待时间（前因子取 10¹³ s⁻¹，仅用于估量级）：

| 能垒 Eₐ | 873 K | 1073 K | 直接 MD 能否看到 |
|---|---|---|---|
| 0.3 eV | 5 ps | 3 ps | 能（H、CO 扩散） |
| 0.5 eV | 0.08 ns | 0.02 ns | 能 |
| 0.7 eV | 1.1 ns | 0.2 ns | 勉强（多位点、多条轨迹） |
| 0.9 eV | 16 ns | 1.7 ns | 873 K 基本看不到 |
| 1.1 eV | 0.2 µs | 15 ns | 看不到 |
| 1.5 eV | 46 µs | 1 µs | 看不到 |

DFT 文献里 Ni 上 CH₄ 首个 C–H 断裂、C+O→CO、C–C 偶联、C 进入次表层的能垒多数在 0.7–2 eV 之间，所以：

- **直接 MLMD 可以研究**：吸附/脱附、H/C/O 扩散、H₂ 生成、台阶上较快的 C–H 断裂（需要高压气体加速碰撞）、MEO 层与 Ni 颗粒的结构动力学、烧结前驱过程（原子剥离、形貌涨落、颗粒接触）。
- **必须用增强采样**：C+O/CH+O 气化、C–C 偶联（成焦）、CO₂ 在平台上的解离、C 溶入次表层。用 OPES/Metadynamics 在真实温度下算**自由能垒 ΔF‡(T)**。它自带非谐、构型熵和动态位点重构的效应，这正是"动态催化"相对静态 NEB 的价值。
- **升温加速（1400 K）** 只用来发现事件和采集训练数据，结论必须回到 873/1073 K 用 OPES 复核。
- 最后用 ΔF‡ 搭**微观动力学模型/kMC**，得到 TOF 和积碳、消碳速率，再和实验的 g g<sub>Ni</sub>⁻¹ h⁻¹ 数据比较。

---

## 2. 关键技术决策（先定下来，后面不要再改）

### 2.1 DFT 参考层级

- **推荐 PBE+D3(BJ)，自旋极化，全局不加 U**。体系里同时有金属 Ni⁰ 和氧化物中的 Ni²⁺，Ni 还会在两者之间转化。按元素加 U 会让"Ni 附近有没有 O"变成能量标度的突变，ML 势学不好，反应能也会失真。
- 只在少量构型上做敏感性检验：把 MEO 中的 Ni 定义成单独的 POTCAR 物种，只给它加 U（模板见 `dft_templates/INCAR_label`），看 Ni²⁺ 磁矩、电荷和关键反应能变化多少。
- 如果最终选择 UMA 的 `oc20` 头作为基础模型，就改用 RPBE(-D3)，和它的预训练数据保持一致。
- 训练用**力一致的自由能（TOTEN）**，ISMEAR=0、SIGMA=0.05；整个数据集的 ENCUT、smearing、D3、U 必须完全一致。

### 2.2 基础模型：先横评，再微调（微调不可省）

候选模型（都能通过 `scripts/calculators.py` 统一调用）：

| 模型 | 训练数据 | 和本体系的匹配度 |
|---|---|---|
| MACE-MP-0 small/medium | MPtrj（PBE/PBE+U，体相为主） | 氧化物可以，**Ni–O 吸附很差（见下）** |
| MACE-MPA-0 / MACE-OMAT-0 | sAlex+MPtrj / OMat24（含非平衡构型） | 高温结构更稳，但仍是 PBE+U 体系 |
| MACE-MATPES-PBE-0 | MatPES（PBE，**不加 U**） | 金属 Ni 与氧化物同处一个能量标度，是最一致的微调起点；缺吸附物数据 |
| MACE-MH-1（多头） | OMat/OMol/OC20/MatPES | 跨体相、分子、表面；需选头（如 `omat_pbe`） |
| UMA-S 1.2（fairchem，`uma-s-1p2p1`）`oc20` / `oc22` / `omat` 头 | OC20（金属表面+吸附质，RPBE，无自旋、无氧化物）/ OC22（氧化物，PBE+U）/ OMat24 | 表面反应最接近；每次只能微调一个头 |
| SevenNet-Omni、MatterSim、ORB-v3（只用能量守恒版本）、DPA 系列 | 各自多数据集 | 一起横评 |

上表后四行是 2026 年 11 月按各软件包的发布说明补充的，这些模型尚未在本仓库测试。没有任何一个现成的头在同一 DFT 水平上同时覆盖金属 Ni、吸附物和含 Ni²⁺ 的氧化物。

**本仓库已经实测到的问题**（`data/benchmark_umlip.json`，Ni(111) p(3×3)，见第 6 节）：

- MACE-MP-0 small 和 medium 给出的 O* 吸附能（以 H₂O−H₂ 为参考）是 **+5 到 +6 eV**。物理上这个值应在 0 附近略负（O/Ni(111) 相对 ½O₂ 约 −2.5 到 −3 eV，H₂O 生成能约 −2.5 eV）。也就是说模型认为 O 几乎不吸附在 Ni 上，CO₂→CO*+O* 被算成 +5 到 +7 eV。这和 MP 训练集"含 O 的 Ni 体系都用 GGA+U"的能量标度问题一致。
- 同一个量，small 和 medium 之间相差可达 1.5 eV。干重整反应能算成 +4.7 到 +5.0 eV，实验值约 +3.2 eV（ΔH = 2.56 eV 加零点能修正）。
- 委员会力偏差（small vs medium，在 S0 模型的 873 K 帧上）：MEO 区最大（O、Si 原子最高 0.8–0.9 eV/Å）；Ni/MEO 界面 2.6 Å 内的 Ni 平均 0.32 eV/Å，远离界面的 Ni 平均 0.22 eV/Å；气相分子不超过约 0.33 eV/Å。

结论：**任何通用势都不能直接用于本体系的反应能量学**。它们的价值在于：(1) 作为微调起点（数据效率高）；(2) 在微调之前做结构退火和构型探索。

### 2.3 软件栈

- 原型/小体系：ASE + 本仓库脚本（`run_md.py` 已带分阶段协议、冻结区、气体反射墙、在线物种统计、PLUMED 接口）。
- 生产：LAMMPS（MACE 通过 ML-IAP/KOKKOS 或 symmetrix 接口，`fix plumed`）。GPU 必需。
- 增强采样：PLUMED ≥ 2.9，编译时打开 OPES 模块（conda-forge 版本默认包含）。
- DFT：VASP（模板已给），或 CP2K（做更大的 AIMD 验证）。
- 微调：`mace_run_train` 多头微调（`dft_templates/finetune_mace.sh`）；若用 DPA-3 则用 DP-GEN2。
- 分析：本仓库 `species.py`（图论物种识别、积碳指标、Ni 结构指标）、OVITO（CNA、团簇）、ReacNetGenerator（反应网络）。

### 2.4 用更快的势函数：能换来什么，怎么用

**速度换来的是时间，而能直接看到的能垒只随时间对数增长。** 模拟时间每增加 10 倍，873 K 下能直接观察到的能垒只提高 kT·ln10 ≈ 0.17 eV（1073 K 下约 0.21 eV）。

单条轨迹里能看到约 3 次事件的最大能垒（10 个等价位点，前因子 10¹³ s⁻¹，只看量级）：

| 单条轨迹长度 | 873 K | 1073 K |
|---|---|---|
| 1 ns（MACE-medium 级别，约 1000 原子） | 0.78 eV | 0.96 eV |
| 10 ns | 0.96 eV | 1.18 eV |
| 100 ns（快 100 倍） | 1.13 eV | 1.39 eV |
| 1 µs | 1.30 eV | 1.60 eV |

所以快 10–100 倍时：
- **能多做的**：台阶上 CH₄ 活化、C 扩散、部分 C–C 偶联可以在 1073 K 直接观察并做统计；烧结和颗粒合并可以跑到 10–100 ns；颗粒可以放大到 5–8 nm，甚至直接模拟 15 nm 全尺寸颗粒（Ni 约 15 万原子 + 1.5 nm MEO 壳 + 载体，共约 30 万原子，可跑 ns 量级）；OPES 可以开更多副本，收敛更快。
- **做不到的**：873 K 下 C+O→CO（1.5 eV 以上）仍然看不到，增强采样省不掉。速度也不解决精度问题。§2.2 里 Ni–O 能量标度出错是数据的问题，而小模型通常比大模型更不准。

可选的快速模型：

| 类别 | 例子 | 相对 MACE-medium 的速度（粗估，以实测为准） | 是否适合本体系 |
|---|---|---|---|
| 同一模型的推理加速 | float32、cuEquivariance、LAMMPS-Kokkos、torch.compile | 2–5 倍 | 应该先做，不影响精度 |
| 更小的通用模型 | MACE-MP-0 small 等 | 约 2–3 倍 | 只适合预探索。small 已实测有同样的 Ni–O 问题 |
| **蒸馏/专用小模型** | Allegro 或 NequIP（小参数）、DeePMD se_e2_a + 模型压缩、ACE（pacemaker）/GRACE-FS、MTP | 约 10–100 倍 | **推荐**：用同一套 DFT 数据训练，再用精确势标注大量增广构型 |
| 非保守模型（直接预测力，没有能量梯度约束） | 部分 ORB、eqV2 的 direct 变体 | 快 | **不推荐**：能量不守恒，自由能和重加权都不成立 |
| 经典反应力场 | ReaxFF | 100 倍以上 | 不建议用于主结果：Ni/Mg/Al/Si/O/C/H 七元素没有一致的参数，需要重新拟合 |

推荐做法：**两级势（teacher–student）**
1. **精确势 A**：微调后的 MACE/UMA（§2.2）。所有要报告的能量数字都由它（或 DFT）给出。
2. **快速势 B**：Allegro / 压缩 DeePMD / ACE 三选一（用 DP-GEN 生态就选 DeePMD，需要 CPU 集群就选 ACE），在 A 的 DFT 训练集上训练。再用 A 给 B 跑出的大量 MD 构型打标签扩充数据。这一步只花 A 的推理时间，不增加 DFT 成本。B 本身也要走主动学习：B 的委员会偏差大时，先让 A 标注；A 的委员会偏差也大时，再交给 DFT。
3. **验证 B**：在反应构型留出集上，相对 A 的力 MAE ≲ 50 meV/Å，关键 NEB 能垒相差 < 0.1 eV；和 A 跑同样的短 MD，比较 RDF、扩散系数和事件统计。
4. **分工**：B 负责长时间无偏 MD、大颗粒、多副本和 OPES 采样；A 负责抽检和重加权。对 B 采到的构型计算 ΔE = E<sub>A</sub> − E<sub>B</sub>，用自由能微扰修正：ΔF<sub>A</sub> = ΔF<sub>B</sub> − kT ln⟨exp(−ΔE/kT)⟩<sub>B</sub>。ΔE 的涨落要在约 1–2 kT 以内才可靠。B 跑出的新事件都要用 A 复核。
5. 进阶：多时间步积分（内层每步用 B，外层每 n 步用 A 修正力）。

其他提速办法：
- 不含 H 的体系（烧结、MEO 退火）用 1 fs 步长。
- 含 H 的体系把 H 质量换成 D 后可以用 1 fs，但只用于结构和烧结研究，不用于速率，因为它会改变动力学同位素效应。
- 统计事件时，多条独立副本比一条长轨迹更划算，并行效率也更高。

---

## 3. 模型体系（已建好，`structures/`，渲染图在 `renders/`）

| 编号 | 体系 | 原子数 | 用途 |
|---|---|---|---|
| S0 | 小 Ni(211)+MEO（1 nm）+孔 + 6 CH₄/6 CO₂；裸 Ni(211) 对照 | 627 / 336 | CPU 冒烟测试 |
| **S1** | **Ni(211) 24×27 Å，5 个台阶单元 + 1.05 nm MEO 膜，中间开 1.2 nm 孔露出台阶 Ni + 12 CH₄/12 CO₂**；裸 Ni(211) 对照；以及不含气体的版本（先退火） | 1367 / 756 | Q2、Q3：反应 MD、CH₄/CO₂ 活化的 OPES |
| **S2** | 同 S1，不含气体，孔内预置 4 C* + 4 O*；裸 Ni 在同一圆盘内放同样的 C/O | 1279 / 668 | **Q1：气化 vs 成焦的二维自由能面** |
| S3 | Ni(111) 平台版本的 S1 | 1471 / 756 | 平台与台阶对比 |
| **P1** | **3 nm 球状 Ni（高指数表面，底部截出 (100) 接触面）+ 1 nm 多孔 MEO 壳（6 个孔）+ Mg(Al)O(100) 载体**；裸颗粒对照 | 9258 / 5496 | Q4、Q5：稳定性、Ni 迁移、MEO 完整性 |
| P2 | 两个 2.5 nm 颗粒并排（有壳 / 无壳） | 12952 / 7386 | Q4：颗粒合并（烧结） |

设计依据和局限：

- **为什么用 Ni(211)**：论文强调球状颗粒的高指数/台阶表面；211 同时含台阶和 (111) 平台。3 nm 球形颗粒（P1）本身就是天然的高指数表面。
- **为什么开孔**：论文的 N₂ 吸附说明 MEO 多孔，"表面金属 Ni 被 MEO 隔离、限域"。孔半径（默认 6 Å）是**关键可调参数**，建议扫描 4/6/9 Å，研究暴露 Ni 团簇（ensemble）尺寸对 C–C 偶联的影响。这一点直接检验论文"隔离 Ni 原子、避免开放晶面上快速形成 C–C 键"的说法。
- **MEO 组成**：论文只给了元素种类，没有给比例。默认按"高熵"假设取等摩尔 Mg:Ni:Al:Si，用阳离子空位保证形式电中性（约 27% 阳离子空位）。建议另做富 Mg（例如 0.45/0.2/0.2/0.15）和富 Ni 两组作敏感性分析，并用 EDS/XPS 比例校准。
- **初始 MEO 是理想的随机岩盐外延膜**（应变 ≤ 3%），**必须先退火**（`protocols/meo_anneal.json`：只放开 MEO 在 2000 K 熔化 → 淬到 1073 K → 整体在 1073 K 退火 100 ps → 降到 873 K），模拟"原位演化"的 MEO，再用 `add_species.py` 加气体或吸附质。退火后要检查：条纹间距是否约 2.2 Å、Ni²⁺ 是否留在层内、孔是否保持。
- **尺寸**：实验颗粒 15 nm（约 16 万个 Ni 原子 + 约 10 万个 MEO 原子），无法直接模拟。局部化学用 S 系列（局部环境等效）；颗粒层面的动力学用 2.5–3 nm 颗粒，**只比较有壳/无壳的趋势**，不要外推绝对烧结速率。
- **气体压力**：1 atm、873 K 下，一个 25×27×24 Å 的气体层里平均只有约 0.1 个分子。模型里放 24 个分子相当于 10²–10³ bar，所以反应 MD 的结果要按"每次碰撞的反应概率"来读，不能当作真实时间速率。

---

## 4. 分阶段工作流

### 阶段 0（1–2 周）：环境、收敛性、基线
1. 在 GPU 机器上装好 MACE/fairchem/LAMMPS/PLUMED/VASP；跑通 `build_models.py`、冒烟测试（第 7 节命令）。
2. DFT 收敛测试：ENCUT、k 点、slab 厚度、真空层；分子（CH₄、CO₂、CO、H₂、H₂O）、体相 Ni、MgO、NiO 的晶格常数与磁矩。
3. 定下 2.1 节的所有 DFT 参数，此后冻结。

### 阶段 1（2–3 周）：DFT 基准集 + 通用势横评
- **B1** 气相：干重整、RWGS、Boudouard 反应能。
- **B2** 吸附能：H、C、O、OH、CH₁₋₃、CO、CO₂ 在 Ni(111)、Ni(211)、Ni/MEO 孔边缘、MgO(100)/MEO 表面（含 CO₂→碳酸盐）。
- **B3** 能垒（CI-NEB）：CH₄→CH₃+H（平台/台阶/周界/Ni²⁺），CH→C+H，CO₂→CO+O，C+O→CO，CH+O→CHO，C+C→C₂，C 表面→次表层。
- **B4** 873 K 下 5 ps DFT-AIMD（约 200 原子的 S1 局部模型）的力对比。
- 运行 `benchmark_umlip.py --calcs ... --dft-ref dft_ref.json --parity b4.extxyz`。
- 合格线（建议）：吸附能 MAE < 0.1 eV，能垒误差 < 0.15 eV，反应构型力 MAE < 0.1 eV/Å。预期没有一个通用势能达标，据此选出误差最小、最稳定的 1–2 个作为微调起点。

### 阶段 2（4–6 周）：主动学习微调（核心工作量）
1. **初始数据**（约 500–800 帧）：B1–B4、NEB 路径上的所有图像、带随机扰动/应变的结构、MEO 熔淬快照、Ni/MEO 界面、孤立分子、单原子（用于 E0）。DFT 单点用 100–250 原子的小 cell（用同一套 builder 生成 12–15 Å 的版本），少量全尺寸 S 帧。
2. **循环**（每轮约 1 周）：用当前模型做探索 MD（873/1073/1400 K；带气体；不同 C/O 覆盖度；颗粒）→ `select_frames.py` 用 4 个不同种子的微调模型做委员会，选出 0.1 < σ<sub>F,max</sub> < 0.6 eV/Å 的帧 → DFT → 重新训练 → 下一轮。**OPES 探索到的过渡态附近构型一定要加入训练集**，反应区是通用势最弱的地方。
3. 收敛判据：新探索轨迹中超过 95% 的帧 σ<sub>F,max</sub> < 0.1 eV/Å；B3 能垒误差 < 0.1–0.15 eV；和 DFT-AIMD 的 RDF、振动态密度（VDOS）、H/C 均方位移一致。
4. 预计 5–8 轮，总计 2000–4000 帧 DFT。

### 阶段 3（4–8 周）：生产模拟

| 子任务 | 模型 / 协议 | 主要输出 | 对应问题 |
|---|---|---|---|
| 3a MEO 结构与熵稳定 | S1/P1 nogas + `meo_anneal.json`；另做 SQS 体相 MEO vs 二元氧化物 | 混合焓 ΔH<sub>mix</sub>、−TΔS<sub>conf</sub>（等摩尔 4 元 R ln4）、层内 Ni²⁺ 去留、孔径演化、2.2 Å 间距 | Q5 |
| 3b 烧结 / 重构 | P1、P2 + `sinter_1073K/1173K.json`，各 ≥ 3 条独立轨迹 | Ni 回转半径、平均配位数、剥离 Ni 原子数、颗粒合并时间、MEO 被破坏的程度；有/无表面 C、H 覆盖 | Q4 |
| 3c 反应 MD | S1/S3（退火后加气体）+ `mdr_873K/1073K.json`；加速版用 1400 K | 物种时间序列、首次活化位点统计（台阶/周界/Ni²⁺）、反应网络 | Q2、Q3 |
| 3d 自由能垒 | `make_plumed.py --mode ch4/co2` + `opes_873K.json`（同样做 1073 K） | ΔF‡(873/1073 K)：CH₄ 在台阶/平台/周界/Ni²⁺ 上；CO₂ 在 Ni/周界/MgO 上 | Q2、Q3 |
| **3e 积碳竞争** | **S2（MEO vs 裸 Ni）+ `--mode coke`（CV：N(C–O)、N(C–C)）+ `coke_873K.json` 无偏 MD** | **二维自由能面 F(N<sub>CO</sub>, N<sub>CC</sub>)；ΔF‡(C+O→CO) vs ΔF‡(C+C→C₂)；C 溶入次表层；晶格 O 参与（`ncolat`）** | **Q1** |

3e 的核心预期：如果论文的机理成立，MEO 限域的小 Ni 区域上，C 原子的迁移和聚集受到孔边缘阻挡，C–C 偶联自由能垒应该升高，或者 C₂ 在热力学上不稳定；同时孔边缘的晶格 O 或 O* 提供更近的氧，C+O 气化能垒应该降低。裸 Ni 上应该相反。要扫描孔径和 C 覆盖度，看这个"交叉"出现在什么条件下。

### 阶段 4（3–4 周）：尺度上推与实验对照
1. **微观动力学模型 / kMC**：以 ΔF‡(T) 和吸附自由能为输入，计算 TOF、CH₄/CO₂ 转化率随温度的变化（对照 Fig. 2a）、H₂/CO 比（RWGS 使其 < 1，对照 Fig. S7）、积碳与消碳速率；做速控度（DRC）分析。
2. **换算到实验单位**：r<sub>C</sub>[g<sub>C</sub> g<sub>Ni</sub>⁻¹ h⁻¹] = TOF<sub>C</sub>[s⁻¹] × 3600 × D × f<sub>exp</sub> × M<sub>C</sub>/M<sub>Ni</sub>。D 是 15 nm 颗粒的分散度（约 0.07），f<sub>exp</sub> 是 MEO 下仍暴露的表面 Ni 比例（从 P1 退火结构或 H₂ 化学吸附得到）。对照实验：OxNi 0.24（积碳）/0.51（消碳），参比 Ni 2.39/1.79。**定性复现"OxNi 消碳 > 积碳、参比 Ni 积碳 > 消碳"就是一个很强的结果。**
3. 其他可对照的实验量：原位 TEM/XRD 粒径不变；Raman 显示少量无定形碳；XPS Ni²⁺（用 Bader 电荷/磁矩作代理）；MEO 条纹 2.2 Å、Ni(111) 2.0 Å。

---

## 5. 风险与对策

| 风险 | 对策 |
|---|---|
| 通用势中 Ni–O 能量标度错误（已实测） | 统一无 U 的 DFT 数据微调；所有含 O 的吸附能都要做 DFT 抽检 |
| ML 势不显含磁性/长程静电（金属–氧化物界面有电荷转移） | DFT 标签统一磁矩初始化；界面误差偏大时尝试更大截断半径或带电荷/长程项的模型；用 Bader/磁矩做后验 |
| 高压气体、温度加速可能改变机理 | 速率只用 OPES 在真实 T 下得到的 ΔF‡；无偏 MD 只用于发现事件和统计位点 |
| 生产段恒温器扰动反应放热的耗散 | 平衡段用 Langevin，生产段用弱耦合的 Bussi（τ = 100 fs）；关键事件用 NVE 片段复核 |
| 经典核运动忽略零点能（C–H 断裂受影响约 0.1–0.2 eV） | 对关键能垒加谐振 ZPE 修正；有条件时用 i-PI 做 PIMD 检验 |
| MEO 组成、孔径未知 | 做 2–3 种组成 × 3 种孔径的敏感性扫描，结论只报告稳健的趋势 |
| 统计误差 | 每个条件至少 3–5 条独立轨迹；OPES 检查 ΔF 随时间的收敛并做分块误差分析 |

---

## 6. 计算量与时间（粗估，请用 `*_timing.json` 实测更新）

- 本容器实测（4 核 CPU，MACE-MP-0 small，float32）：627 原子 1.8 s/步，即 0.024 ns/天。**CPU 不能做生产**。
- GPU（A100/H100 量级）粗估：约 1000 原子时 1–10 ns/天，约 1 万原子时 0.1–1 ns/天，取决于模型大小和 cuEquivariance 等加速。
- 总量：探索与主动学习约 500 GPU·h；生产（3a–3e，含多条轨迹和约 20 个 OPES 自由能面）约 2000–5000 GPU·h。DFT 标签 2000–4000 帧 × 30–100 核时/帧（100–250 原子，自旋极化）≈ 10–30 万核时。
- 周期：约 4–6 个月（阶段 0–1：1 个月；阶段 2：1.5 个月；阶段 3：2 个月；阶段 4：1 个月）。

本仓库对 MACE-MP-0 small/medium 两个基础模型的实测结果（Ni(111) p(3×3)，单位 eV；吸附能参考 CH₄、H₂、H₂O、CO、CO₂ 气相）：见 `data/benchmark_umlip.json` 和 README 中的表。

---

## 7. 本仓库已经可以运行的内容

```bash
cd oxni_mdr_mlmd/scripts
python build_models.py                                   # 生成全部模型与渲染图
# 冒烟测试（CPU，约 5 分钟）
python run_md.py ../structures/S0_smoke_ni211_meo_pore_gas.extxyz --calc mace:small \
       --protocol ../protocols/smoke.json --out ../runs/smoke/S0_meo
# 1) MEO 退火（GPU）
python run_md.py ../structures/S1_ni211_meo_pore_nogas.extxyz --calc mace:medium-omat-0 \
       --device cuda --protocol ../protocols/meo_anneal.json --out ../runs/S1_anneal
# 2) 加气体 / 预置 C*、O*
python add_species.py ../runs/S1_anneal_final.extxyz --gas 12 12 --out ../structures/S1a_gas.extxyz
python add_species.py ../runs/S1_anneal_final.extxyz --adatoms C=4,O=4 --out ../structures/S2a_CO.extxyz
# 3) 反应 MD / 积碳 MD
python run_md.py ../structures/S1a_gas.extxyz --calc mace:<微调模型> --device cuda \
       --protocol ../protocols/mdr_873K.json --out ../runs/S1a_873K
# 4) 增强采样
python make_plumed.py ../structures/S2a_CO.extxyz --mode coke --out ../plumed/S2a_coke.dat
python run_md.py ../structures/S2a_CO.extxyz --calc mace:<微调模型> --device cuda \
       --protocol ../protocols/opes_873K.json --plumed ../plumed/S2a_coke.dat --out ../runs/S2a_opes
python fes_reweight.py ../runs/COLVAR --cv nco ncc --T 873 --png ../runs/S2a_fes.png
# 5) 主动学习选帧、基准、作图
python select_frames.py ../runs/S1a_873K.extxyz --models mace:m1.model mace:m2.model mace:m3.model \
       --out ../al/iter1.extxyz
python benchmark_umlip.py --calcs mace:medium-omat-0 uma:uma-s-1p1/oc20 --dft-ref ../data/dft_ref.json
python plot_runs.py ../runs/S2a_meo_log.csv ../runs/S2a_bare_log.csv --labels OxNi bare \
       --cols "CO(g)" C_in_Cn largest_C_cluster C_dissolved --out ../runs/S2_compare.png
```

---

## 8. 接下来两周的清单

1. [ ] 确定 DFT 软件和参数（2.1），完成收敛测试。
2. [ ] 准备 GPU 环境，用 `benchmark_umlip.py` 跑 MACE-OMAT-0、MACE-MPA-0、UMA-S（oc20/omat）、SevenNet、DPA-3。
3. [ ] 计算 B1–B3 基准集（约 40 个吸附构型 + 8–10 条 NEB），写入 `data/dft_ref.json`。
4. [ ] 用选出的基础模型跑 `meo_anneal.json`（S1、S3、P1 各一次），检查 MEO 结构是否合理；从退火轨迹中抽取第一批 DFT 帧。
5. [ ] 和实验同事确认：MEO 的 EDS 半定量比例、孔径/BET、H₂ 化学吸附得到的暴露 Ni 比例。这三个数直接决定模型参数。

## 参考
- He et al., ACS Catal. 2021, 11, 12409（本研究对象）
- Bengaard et al., J. Catal. 2002, 209, 365（Ni 上台阶位点、蒸汽重整与积碳）
- Helveg et al., Nature 2004, 427, 426（Ni 台阶上碳纳米纤维生长的原位观测）
- Batatia et al., MACE-MP-0 基础模型（arXiv:2401.00096）；Meta FAIR, UMA（2025）
- Invernizzi & Parrinello, J. Phys. Chem. Lett. 2020, 11, 2731（OPES）；PLUMED consortium, Nat. Methods 2019
- Bussi, Donadio & Parrinello, J. Chem. Phys. 2007, 126, 014101（CSVR 恒温器）
