# 交接文档：费托反应网络枚举与 α 控制步骤分析

本文汇总到目前为止的全部工作、数据、结论和下一步，方便在本地（有 GPU 的机器）继续。
代码都在 `ftc_network/` 下，分支 `claude/ftc-reaction-network-enum-t2ihw9`。

## 0. 最终要回答的问题

用一个比较完备的费托（FT）反应网络，结合微观动力学模型（MKM），找出哪一步（或哪一类反应）
控制链增长概率 α。当前阶段不做 DFT：先用随机但热力学一致的能量把整套流程跑通，
看结论是由网络结构决定（随机参数下也稳定），还是必须依赖准确能量（之后可用 UMA 等 MLIP）。

## 1. YARP 调研结论

* PyPI 上的 `yarp` 是同名的响应式编程库，不是反应网络工具。要用的是 Savoie 组的
  YARP 3.0.1：`git+https://github.com/Savoie-Research-Group/yarp.git`。
* 实测 `enumerate_products` 的问题：
  * concerted 模式忽略 `n_form`（强制 n_break = n_form），`b1f0` 和 `b1f1` 结果一样；
  * sequential 模式的 `b0f1` 返回空（连 CH3 + CH3 都没有）；
  * 自由基加成到 π 键（H + CO → HCO）不在 `b0f1` 里，只在 `b1f1` 里出现；
  * 会生成带电产物（如 `[C-]#[OH+2]`）。
* **CO 的 C–O 解理在 YARP 里的分类**：按 YARP 的定义（数 σ 连接）是 `b1f0`，逆反应 C + O → CO 是 `b0f1`。
  但 YARP 默认根本不生成它：sequential 的 `break_bonds` 默认 `break_higher_order=False`，只断能写成单键的键，
  CO 在 YARP 里是 `[C-]#[O+]` 三键；`enumerate_products` 和 `yarp-init` 的 YAML 都没有暴露这个参数；
  concerted 每个键最多降一级且强制成键数等于断键数，相当于 b3f0，做不到。
  同理 C=C、C=O 的完全断裂（C2H4 → 2CH2、HCO → CH + O、CH2O → CH2 + O）都不会生成；
  单键 C–O（CH3OH → CH3 + OH）会生成。
* 因此自己写了一层很薄的枚举（`ftc_enum/`），沿用 YARP 的约定（b/f 数 σ 连接，键级由 Lewis 结构重新确定），
  YARP 只用来算 Lewis 结构和校验去重（`--yarp`）。

## 2. 气相网络枚举工具（`ftc_enum/`、`run_enumeration.py`）

* 物种 = 重原子（C、O）组成的树 + 每个重原子上的 H 数。H· 和 H2 用一个 'H' 节点表示。
  不成环、不允许 O–O、不带电荷；不区分立体异构和单重态/三重态。
* 键级：把剩余价在树上贪心配成 π 键（对树是最优的），剩下的叫自由价（自由基电子 + 卡宾）。
  例如 CH3 = 1、CH2 = 2、CO 写成 `[C]=O` 为 2。
* 规范键：连接关系的规范 SMILES，H 数用同位素标签编码（早期版本没区分 H 数导致 CH4 和 C 混同，已修复并加测试）。
* 规则：
  * `b1f0`：断 1 个 σ 键（均裂、β-断裂）；
  * `b0f1`：成 1 个 σ 键（自由基结合、加成到 π 键）；
  * `b1f1`：断 1 成 1、两键共用一个原子（分子内 H 迁移、重原子基团迁移，迁移距离不限于 1,2；分子间 H 抽提/歧化）。
    分子间 `b1f1` 默认要求接受方是开壳层，分子间重原子基团转移默认关闭。
* 双分子配对：A + A 总是算；任意物种 + 小分子池（≤ 3 个重原子）；≤ 5 个重原子的物种之间全部两两配对；
  大 + 大不算（否则配对数 ~ 物种数²）。
* 分档限制（`config.yaml` 的 `limits`）：

| 重原子数 | 自由价 | π 键 | O 原子 |
|---|---|---|---|
| 1 | ≤ 4 | 不限 | ≤ 1 |
| 2–5 | ≤ 2 | 不限 | ≤ 2 |
| 6–25 | ≤ 1 | ≤ 1 | ≤ 1 |

  支链上限全局为 2（`max_branches`，可按档覆盖）。≤ 5 个重原子时 2 个支链从不起作用（最多就是新戊烷）。
* 输出：`species.csv`、`reactions.csv`（或 `--compact` 的 gzip 版）、`summary.json`；加 `--yarp` 时另有 `species_yarp.csv`。
* `--nproc 1`（全串行）与 `--nproc ≥ 2` 的处理顺序不同：网络按规范键完全相同，但 species_id 编号不同，
  且约 2.7% 的开壳层物种（N = 9：347 / 12,773）的 SMILES 写成另一种共振式（如 `C#C[C]=O` 与 `[CH]=C=C=O`）。
  自由价、π 键数不受影响；但按 SMILES 做子结构匹配的 `tools/classify_species.py` 计数会小幅变化
  （N = 9：炔 232 → 225、累积烯 63 → 71、乙烯酮 53 → 57）。本文的统计都来自 `--nproc 4`。
  若需要与处理顺序无关的 SMILES，可在注册物种时先用 `species_from_key` 按规范原子顺序重建再生成 SMILES（尚未实现）。
* 性能：同一次调用中产物规范化做缓存、子进程只回传一次物种对象、处理完的大分子瘦身（保留原子顺序的紧凑表示）、
  gzip 输出。N = 9 从 7.5 分钟降到 64.5 秒（4 核），结果与优化前逐条一致。
  代码审查发现"瘦身后从规范键重建会改变原子顺序、从而改变 id 编号"，已修复为保留原顺序，并与优化前代码逐字节比对一致。
* 测试：`python -m pytest -q tests`（15 个，全部通过）。

## 3. 规模（实测，默认配置，≤ 2 个支链，4 核）

| | N = 9 | N = 12 |
|---|---:|---:|
| 物种 | 12,773 | 180,623 |
| 反应 | 7,288,960 | 44,138,680 |
| b1f0 / b0f1 / b1f1 | 57,212 / 57,212 / 7,174,536 | 771,282 / 580,397 / 42,787,001 |
| 耗时 | 64.5 s | 15.6 min |
| 压缩反应表 | 33 MB | 212 MB |

* 物种数都与 `count_species.py` 直接枚举的结果完全一致：闭包得到的物种集合只由过滤条件决定，
  所以"先列出全部物种、再连线"与逐层迭代在物种层面等价。
* 按重原子上限截取的子网络见 `results/subnet_scale_n12.csv` 和 `results/scale_vs_max_heavy.png`：
  ≤ 5 时反应数随节点数增长略快于平方（约为节点数的 2.4 次方，分子变大后每对物种的反应数也增加；
  N = 5：1,135 个节点、479 万个反应，平均每个节点 4,218 个）；
  之后每个新节点只带来约 220 个反应，N = 12 时平均 244 个。
* 从 N = 10 起 b1f0 多于 b0f1：例如 C10 → C4 + C6 会枚举，但 C4 + C6 的结合按配对规则不算。
* 每层物种数（恰好 N 个重原子）：

| N | 直链 | ≤ 1 个支链 | ≤ 2 个支链 |
|---|---:|---:|---:|
| 8 | 648 | 2,031 | 2,801 |
| 10 | 1,827 | 9,420 | 18,932 |
| 12 | 4,416 | 34,931 | 103,431 |
| 13 | 6,555 | 62,721 | 223,088 |
| 19 | 44,027 | — | — |

* 外推到 C25：≤ 2 个支链约 5 × 10⁸ 个物种、约 10¹¹ 个反应（约 4000 核时、约 500 GB 压缩反应表），
  不可行；≤ 1 个支链约 3 × 10⁷；`config_c25_tiered.yaml`（≤ 12 个重原子 2 个支链，更大的只留直链）约 1.0 × 10⁶ 个物种、
  约 3 × 10⁸ 个反应、约 10 核时；全部直链约 8 × 10⁵。
* 网络结构（N = 9）：每个物种参与的反应数中位数 231，最多的是 C3 开壳层物种（如 `[CH]=C=C` 参与 14.6 万个），
  因为它们属于小分子池、要和所有物种配对。

## 4. 化学检查

* N = 12 网络的类别统计（`tools/classify_species.py`，结果 `results/classes_n12.json`）：
  开壳层 89%；含醚键 49%；非甲基支链 31%；烯醇醚 12%；季碳 11%；烯醇 5%；
  其余卡宾/双自由基、炔、两个 O 在同一个碳上或不同碳上、累积烯、乙烯酮、酯、羧酸、无 H 碳簇都在 ≤ 5 个重原子里。
* ≤ 5 个重原子：1,135 个物种中 945 个开壳层（中间体），190 个闭壳层；
  闭壳层中 132 个带醚、酯、烯醇、炔、累积烯/乙烯酮、两个 O 分在不同碳上、无 H 碳簇、季碳等结构，
  剩 58 个常规分子（H2O、CO2、CH4、C2–C5 的烷烃、烯烃、醇、醛、酮、羧酸等），
  其中还有 5 个偕二醇和若干二烯、不饱和醛酮属于少见副产物。
* 乙二醛 O=CC=O 合法（稳定闭壳层分子），最直接的来源是 HCO + HCO 结合；
  网络里 HCO + HCO 同时有歧化（→ CO + CH2O）等 4 个通道。它进来是因为 ≤ 5 档允许 2 个 O
  （本意是保留 CO2、HCOOH、HCOO、HOCO）。若要排除，可把这一档 `max_O` 改成 1，
  或加一条"多个 O 只能连在同一个碳上"的过滤规则（尚未实现）。
* 很多不是稳定产物的物种可能是表面中间体（自由价 ≈ 表面键数）：CHx*、HCO*、HCOH*、COH*、酰基、
  乙烯酮/乙烯酮基、烯基、亚乙烯基、次乙基等。曾启动一个多代理工作流逐类评审"公认中间体 / 被提出的中间体 /
  仅气相瞬态 / 不可能"，**该工作流卡住未完成，没有结论**，需要时重新做。
  已知的缺口：COH*（自由价 3）在气相网络 ≤ 5 档（自由价 ≤ 2）中被排除；大分子档自由价 ≤ 1 排除了 C6+ 的亚烷基、次烷基。
  表面网络（第 6 节）已包含这些。

## 5. 讨论过的方法问题

* **先枚举物种、再连线**：比逐层迭代简单。断键/成键反应可以对每个物种的每根键切一刀完整生成（线性），
  还能补全"大 + 大"结合；分子内 H 迁移也是逐物种处理。只有分子间 H 转移需要两两配对，
  可以存成"给 H 半反应 × 接 H 半反应"的隐式形式。真正的限制是数据量（≤ 2 个支链到 C25 约 10¹¹ 条反应，≤ 1 个支链约 10¹⁰）。
* **扩展 + 剪枝**：标准做法是 RMG 的基于速率的模型扩展（核心集 + 边缘集，按 MKM 通量把边缘物种提升进核心集），
  本质是以通量为评分的最佳优先搜索；也可以用束搜索、分支定界、MCTS（更适合找特定路线，不保证完备）。
  注意：按物种剪，正逆反应成对保留；边缘集不删除；降低阈值检验 α 和关键步骤是否收敛；不能按"是否稳定产物"剪。
* **回答 α 问题不需要显式算到 C25**：把 C1–C3/C4 的化学写完整，n ≥ 3 用链长无关的速率常数或模板延伸到 C20 左右即可；
  用选择性控制度 X = ∂ln α / ∂(−G_i/k_BT) 找控制步骤，并在不同温度、H2/CO 下检验。

## 6. 进行中：随机参数 MKM（`mkm/`）

详细设计见 `mkm/README.md`。已完成表面网络生成（n_max = 10：162 个物种、546 个可逆反应）
和随机能量模型（热力学一致，冒烟测试通过）；**未完成**：`kinetics.py`（稳态求解 + 隐式微分灵敏度）、
`run_alpha_study.py`（系综运行、类控制度统计、全局敏感性、作图）和测试。

## 7. 在本地继续：哪些可以用 GPU

注意：平均场 MKM 是常微分方程组（ODE）的稳态问题，不是偏微分方程（PDE）；
只有加入反应器或颗粒内的传质（轴向扩散、孔扩散）时才变成 PDE。

| 工作 | GPU 是否合适 | 建议 |
|---|---|---|
| 随机参数系综的 MKM 稳态（每个样本约 112 个变量） | **最合适** | JAX（`jax.vmap` 批量 + `diffrax` 刚性求解器，或直接批量牛顿迭代）或 PyTorch；雅可比用 `jax.jacfwd`；稳态灵敏度 dθ/dp = −J⁻¹ S diag(r) 是批量线性求解 |
| 选择性控制度的有限差分校验 | 合适 | 同上，批量扰动 |
| 物种枚举去重 | 可以 | 树上的 1 维 Weisfeiler–Lehman 颜色细化是完整同构判据，可批量矩阵化做哈希（用 128 位防碰撞）；但更大的收益来自"直接按规范形式生成 + 编译型代码" |
| 断键/成键、H 转移的反应生成 | 不太需要 | 主要是数据量和 I/O 问题，用半反应隐式存储 |
| 之后用 UMA 等 MLIP 算能量 | 合适 | MLIP 推理本身就跑在 GPU 上 |

**精度提醒**：链物种覆盖度跨很多个数量级、ODE 极度刚性，必须用 float64（JAX 需 `jax.config.update("jax_enable_x64", True)`），
并在对数覆盖度上做牛顿迭代。消费级 GPU 的双精度吞吐通常只有单精度的 1/32–1/64，A100/H100 这类数据中心卡才快；
样本数不大时 CPU 多进程可能同样快，建议先各测一次。

## 8. 本地复现

```bash
git clone -b claude/ftc-reaction-network-enum-t2ihw9 https://github.com/jinlijin6090/Usedcardata.git
cd Usedcardata/ftc_network
python -m venv venv && source venv/bin/activate
pip install rdkit numpy pyyaml pytest matplotlib       # 枚举、测试、作图
pip install -r requirements.txt                        # 另含 YARP 及其依赖（--yarp 需要）
python -m pytest -q tests                              # 15 个测试

python count_species.py --max-heavy 9                  # 只数物种
python run_enumeration.py -o output_n9 --max-heavy 9 --nproc 8 --compact     # 4 核约 1 分钟
python run_enumeration.py -o output_n12 --max-heavy 12 --nproc 8 --compact   # 4 核约 16 分钟
python subnet_stats.py output_n12 -o results/subnet_scale_n12.csv
python plot_scale.py results/subnet_scale_n12.csv -o results/scale_vs_max_heavy.png
# 作图需要中文字体：Linux 上可 `sudo apt install fonts-noto-cjk`（或 fonts-wqy-zenhei），装完删掉 ~/.cache/matplotlib
python tools/classify_species.py output_n12/species.csv results/classes_n12.json
python tools/compare_runs.py output_a output_b          # 比较两次运行

# 表面网络与随机能量（MKM 部分进行中）
python -c "from mkm.surface_network import SurfaceNetwork; print(SurfaceNetwork(n_max=10).build().summary())"
```

大的输出文件（N = 12 的压缩反应表 212 MB）没有放进仓库，按上面的命令重新生成即可。

## 9. 文件清单

| 路径 | 内容 |
|---|---|
| `ftc_enum/species.py` | 物种表示、键级分配、规范键、紧凑打包 |
| `ftc_enum/rules.py` | b1f0 / b0f1 / b1f1 规则 |
| `ftc_enum/network.py` | 闭包迭代、配对规则、多进程并行、瘦身 |
| `ftc_enum/yarp_annotate.py` | YARP Lewis 结构注释与去重校验 |
| `run_enumeration.py` | 气相网络枚举主程序 |
| `count_species.py` | 只数物种（规模估算） |
| `subnet_stats.py`、`plot_scale.py` | 按重原子上限统计子网络、作图 |
| `tools/classify_species.py`、`tools/compare_runs.py` | 物种化学分类、两次运行比对 |
| `config.yaml`、`config_c25_tiered.yaml` | 默认配置（≤ 2 个支链）、C25 可行版 |
| `results/` | N = 9、12 的 summary，子网络规模表与图，N = 12 物种分类统计 |
| `mkm/surface_network.py`、`mkm/energetics.py`、`mkm/README.md` | 随机参数 MKM（进行中）|
| `tests/test_rules.py` | 单元测试 |

## 10. 下一步

1. 实现 `mkm/kinetics.py`：吸附/脱附、速率、稳态求解（BDF + 对数牛顿）、隐式微分灵敏度、加和规则与有限差分检验。
2. 实现 `mkm/run_alpha_study.py`：100–1000 个随机样本，统计主控反应/类别、类控制度分布、全局敏感性并作图。
3. 根据结果决定：收窄能量范围（例如 CO 直接解离的反应能目前 −1.1 ~ +1.0 eV 偏宽）、加入烯烃再吸附、
   换不同条件（T、H2/CO）检验结论稳定性；之后把随机能量换成 UMA 等 MLIP 的能量。
4. 可选：气相网络加"多个 O 只能连在同一个碳上"的过滤；补做中间体评审。
