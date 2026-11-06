# ni_coke_md — Ni 团簇上的碳沉积 / 退火 / C+O 共沉积（ReaxFF）

一个自带运行环境的完整包：拷到服务器上解压就能跑，不需要 root、联网、conda 或 module。
包含三项任务的输入生成、断点续跑、单算例分析、汇总报告和 GIF 动画。

| 任务 | 内容 | 默认参数 |
|---|---|---|
| **task1_flux** | 纯 C 沉积，通量收敛测试 | 插入间隔 10/20/50 ps/C × 种子 1/2/3，每个算例插 120 个 C，沉积后再跑 50 ps |
| **task2_anneal** | 把 task1 的最终结构长时间退火，看碳环能否整理成六元环 | 从 50 ps/C 的 3 个种子出发，1200 K 退火 2 ns |
| **task3_codep** | C + O 共沉积，O 代表 CO₂ → CO + O 留下的氧 | O/C = 0, 0.5, 1, 2 × 种子 1/2/3，C 间隔 20 ps，沉积后再跑 100 ps |

体系：自由的 Ni₁₁₆ 截角八面体团簇（约 1.3 nm，不带载体），1200 K，NVT，时间步长 0.25 fs，
力场为 LAMMPS 自带的 `ffield.reax.FC`（只用其中的 C/Ni/O）。

---

## 1. 目录结构

```
ni_coke_md/
├── README.md                 本文件
├── run_case.sh               运行/续跑单个算例
├── run_local.sh              无调度系统时在工作站上批量运行
├── env/
│   ├── parts/                自带运行环境的 3 个分片（git 用）；或直接放 nicoke_env.tar.gz
│   ├── (nicoke_env.tar.gz)   自带运行环境（Python 3.12 + LAMMPS 2025.7.22 含 REAXFF/OPENMP + numpy/scipy/networkx/matplotlib）
│   ├── install_env.sh        离线安装（解压 + 修正路径 + 自检）
│   ├── activate.sh           激活环境（所有脚本内部会自动 source）
│   ├── check_env.sh          自检 + 测速
│   └── build_env.sh          在有网机器上重新构建 nicoke_env.tar.gz（一般用不到）
├── ff/ffield.reax.FC         ReaxFF 参数文件（LAMMPS potentials 目录原样拷贝）
├── structures/               ni79 / ni116 / ni201 团簇（data 文件，C=1 Ni=2 O=3）
├── md/
│   ├── driver.py             核心驱动：平衡 → 逐个插入原子 → 删除脱附碎片 → 后续弛豫/退火，带断点
│   ├── build_cluster.py      生成其他尺寸的 Ni 团簇
│   └── common.py
├── campaign/
│   ├── campaign.json         ★ 所有任务参数都在这里改
│   └── campaign.py           make / list / status / estimate
├── validate/
│   ├── ff_checks.py          力场合理性检查（键长、结合能、吸附能）
│   └── ff_report_reference.md  我已经跑出来的检查结果
├── analysis/
│   ├── analyze_case.py       单算例：轨迹 → metrics.csv
│   ├── report.py             ★ 汇总所有算例 → analysis_out/REPORT.md + 图 (+ GIF)
│   └── render_gif.py         纯 matplotlib 渲染 GIF（不需要显卡/OpenGL）
├── pbs/                      PBS Pro 提交：pbs.conf / submit.sh / job.pbs
└── examples/                 之前 demo（2 ps/C）的 GIF 和曲线图，作参考
```

---

## 2. 在服务器上的完整流程

### 2.1 安装（一次性，约 1 分钟）

**方式 A：从 GitHub**（环境包已切成 3 个 95 MB 分片放在 `env/parts/`，安装时自动合并并校验）

```bash
git clone --single-branch -b ni-coke-md https://github.com/jinlijin6090/Usedcardata.git
cd Usedcardata/ni_coke_md
bash env/install_env.sh
```

**方式 B：从压缩包**

```bash
tar -xzf ni_coke_md_code.tar.gz && mv nicoke_env.tar.gz ni_coke_md/env/
cd ni_coke_md && bash env/install_env.sh
```

自检最后会打印 `OK` 和 ReaxFF 的测速结果（每 ps 需要几秒）。

**系统要求**：x86_64 Linux，glibc ≥ 2.27（CentOS/RHEL/Rocky 8 及以上、Ubuntu 18.04 及以上都可以；
**CentOS 7 不行**，见第 8 节）。查看 glibc 版本：`ldd --version | head -1`。

### 2.2 力场检查（强烈建议，约 1 分钟）

```bash
source env/activate.sh
python validate/ff_checks.py --threads 4      # 结果写到 validate/ff_report.md
```

### 2.3 生成算例 & 估算机时

```bash
python campaign/campaign.py make              # 生成 runs/<task>/<case>/config.json
python campaign/campaign.py estimate --sec-per-ps 20
```

按 4 线程约 20 s/ps 估算（实际速度见各算例的 `progress.log`）：

| 任务 | 算例数 | 单个最长 | 最长单算例墙钟 | 总核时（4 线程） |
|---|---|---|---|---|
| task1_flux | 9 | 6060 ps | ~34 h | ~650 |
| task2_anneal | 3 | 2000 ps | ~11 h | ~130 |
| task3_codep | 12 | 2510 ps | ~14 h | ~670 |

### 2.4 提交

- **PBS（你们的集群）**：
  ```bash
  vi pbs/pbs.conf                  # 第一次先检查队列名、核数、software 标签
  bash pbs/submit.sh --dry-run     # 预览
  bash pbs/submit.sh               # 提交全部；之后到 walltime 会自动续交，task1 完成后自动提交 task2
  ```
  详见 `pbs/README.md`。
- **没有调度系统的工作站**：`bash run_local.sh -j 4 -t 4 task1_flux`（同时跑 4 个算例，每个 4 线程）。

**提交顺序**：task1 和 task3 互相独立，可以同时提交；
**task2 依赖 task1 中 50 ps/C 的 3 个算例的 `final.data`**，要等它们跑完再提交。
task2 提前启动会直接报错退出，不会产生错误结果。

**线程数**：每个算例 4 个 OpenMP 线程效率最好。体系只有两三百个原子，更多线程基本不会更快。
与其给一个算例 16 线程，不如同时跑 4 个算例。

### 2.5 查看进度

```bash
python campaign/campaign.py status
tail runs/task1_flux/i50_s1/progress.log
```

### 2.6 分析（随时可以跑，未完成的算例会用已有的帧）

```bash
python analysis/report.py --jobs 8           # 所有任务
python analysis/report.py task1_flux --gif   # 只看 task1，并为每组的种子 1 生成 GIF
```

结果在 `analysis_out/`：`REPORT.md`、`task1_flux.png`、`task2_anneal.png`、`task3_codep.png`，以及 `*.gif`。

单个算例：
```bash
python analysis/analyze_case.py runs/task3_codep/OC1_s2      # → metrics.csv, rings_final.csv
python analysis/render_gif.py runs/task3_codep/OC1_s2 --frames 200
```

---

## 3. 断点续跑

- 每个算例每隔 `checkpoint_every_ps`（默认 50 ps）写一次 `ckpt.data` + `state.json`。
- 作业被杀掉（walltime 到了、节点故障）以后，**原样重新运行同一条命令即可续跑**。
  轨迹和事件日志会自动截断到断点处，不会出现重复帧。
- `run_case.sh CASE THREADS MAX_HOURS`：跑满 `MAX_HOURS` 小时后，在下一个断点处干净退出（退出码 3）。
  建议设为作业 walltime 减 1 小时，然后重新提交。
- 已完成的算例（`state.json` 里 `phase = done`）再次运行会直接跳过。
- 续跑时速度从断点文件里读出，热浴（Nosé–Hoover）的内部变量会重新初始化，影响可以忽略。

---

## 4. 每个算例目录里的文件

| 文件 | 内容 |
|---|---|
| `config.json` | 参数（由 campaign.py 生成） |
| `traj.lammpstrj` | 轨迹（id type x y z q），可以直接用 OVITO / VMD 打开 |
| `events.csv` | 每次插入和每次删除脱附碎片的时间、种类（如 `CO`、`CO2`、`O2`、`C`） |
| `progress.log` | 进度和速度（ps/h） |
| `state.json`, `ckpt.data` | 断点 |
| `final.data` | 最终结构（含速度），task2 从这里接着跑 |
| `log.lammps.N` | LAMMPS 日志（每次续跑一个） |
| `metrics.csv` | 分析结果（运行分析后生成） |

---

## 5. 分析指标说明

`metrics.csv` 每行对应一帧：

| 列 | 含义 |
|---|---|
| `C_ins`, `O_ins` | 累计插入的 C、O 原子数 |
| `C_rem`, `C_rem_CO`, `C_rem_CO2`, `C_rem_other` | 累计以脱附碎片形式删除的 C（总数 / 以 CO / 以 CO₂ / 其他） |
| `C_interior` | 位于 Ni 原子凸包内部的 C（亚表层、体相溶解） |
| `C_surface` | 在凸包外、但与团簇相连的 C（表面吸附、表面碳链和碳岛） |
| `C_detached` | 不和 Ni 相连的 C（还没落到团簇上，或刚脱附） |
| `C1`, `C2`, `C3_5`, `C6p` | 所在 C–C 簇大小为 1 / 2 / 3–5 / ≥6 的 C 原子数（C–C 键判据 1.85 Å） |
| `largest` | 最大 C–C 簇的原子数 |
| `r5`, `r6`, `r7`, `rings` | 最小环基中的五元、六元、七元环数，以及环的总数 |
| `sp`, `sp2`, `sp3` | 有 2 / 3 / ≥4 个 C 邻居的 C 原子数（sp2 多说明石墨化） |
| `O_interior`, `O_surface` | O 在团簇内部 / 表面 |

**怎么看各项任务的结果**

- **task1**：横轴是"已沉积的 C 原子数"而不是时间，这样不同通量才能直接比较。
  报告里给出两个最慢通量（20 和 50 ps/C）之间的差值除以种子间标准差：
  |Δ|/σ < 1 说明通量已经收敛，50 ps/C 的结果可以当作"慢通量极限"。
  如果不收敛，就要把 `c_interval_ps` 再加一档（比如 100）。
- **task2**：关注 `r6` 是否上升、`r5`/`r7` 是否下降，以及六元环比例 `f6` 和 sp2 比例的变化。
  2 ns 内不一定能看到成片的石墨烯。如果 `f6` 一直在缓慢上升，可以加长 `anneal_ps`，
  或者在 `T` 里加一个 1400 K 做对比（注意小团簇在高温下会熔化）。
- **task3**：比较不同 O/C 下"留在颗粒上的 C"和"以 CO 形式移除的 C"。
  如果 O/C ≥ 1 时 `C6p` 和环数明显低于 O/C = 0，说明 O 的清碳作用在这个模型里起效。

---

## 6. 模型的物理假设和局限（写论文时要交代）

1. **沉积法本身就是加速**：直接插入 C，跳过了 CH₄ 解离。即使 50 ps/C，通量也比真实情况高几个数量级。
   所以只做定性和相对比较（谁先成核、通量或 O 的影响），不给绝对速率。
2. **O 原子代表 CO₂ 解离后留下的 O**。默认 CO₂ → CO + O 这一步很快，产生的 CO 直接离开，所以不模拟。
   如果需要显式模拟 CO₂，要另外写插入分子的逻辑。
3. **脱附碎片会被删除**：离开团簇超过 6 Å 的碎片（CO、CO₂、O₂、单个 C 等）记录到 `events.csv` 后删除，
   相当于这些产物被气流带走。
4. **没有载体**（Al₂O₃），**没有 H**。
5. **小团簇在 1200 K 下结构比较无序**。碳溶解量可能偏高，不同尺寸的结果不能直接外推。
   可以用 `structures/ni201.data` 对比尺寸效应：改 `campaign.json` 里的 `cluster` 字段。
6. **力场**：见下一节。

---

## 7. 力场的来源和可信度

- `ff/ffield.reax.FC` 来自 LAMMPS 官方 `potentials/` 目录，文件头只引用了 Singh et al., PRB 87, 104114 (2013)（C/F/H 部分）。
  根据 matsci.org 论坛上的说法，其中的 Ni 和 Ni–C–H 参数与
  Mueller, van Duin, Goddard, *J. Phys. Chem. C* 114, 4939 (2010) 一致，**但我没能核实**。
  **Ni–O 参数的来源不明。**
- 检查结果（`validate/ff_report_reference.md`，参考值为近似值）：
  - 可以接受：Ni 晶格常数 3.61 Å（实验 3.52）；Ni 结合能 4.50 eV（实验 4.44）；
    石墨烯 C–C 键长 1.433 Å；CO 键长 1.126 Å、解离能 11.5 eV（实验 11.1）；
    C 在 Ni(111) 上的吸附能 −6.7 eV；CO 在 Ni(111) 上的吸附能 −1.74 eV。
  - **有偏差**：CO₂ 的 C=O 键长 1.25 Å（实验 1.16）；O₂ 键长 1.26 Å（实验 1.21）；
    O 在 Ni(111) 上的吸附偏强（−5.9 eV，GGA 文献值约 −5.0 到 −5.6）。
  - 结论：**task1 和 task2（只有 C/Ni）比较可靠；task3 的 O 化学只能看趋势**。
    如果有更可靠的 Ni/C/O ReaxFF 参数文件，放进 `ff/`，改 `campaign.json` 里的 `ffield` 即可。
    元素顺序要求是 `C Ni O`，脚本里写死了 `pair_coeff * * <ffield> C Ni O`。

---

## 8. 常见问题

**glibc 太老（CentOS 7）**：预编译的 LAMMPS 需要 glibc ≥ 2.27。有两个办法：
(a) 如果集群有 Apptainer/Singularity，在 Rocky 8 容器里跑；
(b) 如果集群已经装了带 REAXFF、OPENMP 和 Python 接口的 LAMMPS，可以直接用它跑 `md/driver.py`。
driver 只依赖 `lammps` Python 模块，以及 numpy、scipy、networkx。

**`Not enough space for bonds`**：driver 里已经设置了 `safezone 3.0 mincap 200`。
如果还是报错，把 `md/driver.py` 里这两个数调大。

**换团簇尺寸**：`python md/build_cluster.py --length 7 --cutoff 2 -o structures/ni201.data`，
然后改 `campaign.json` 里的 `cluster`。

**改参数以后**：已经存在的 `config.json` 不会被覆盖。确实要覆盖，加 `--force`；但正在跑的算例不要覆盖。
推荐的做法是改参数后用新的种子或新的间隔值生成新算例。

**磁盘**：轨迹每帧约 15 KB。默认设置下整个 campaign 的轨迹加起来不到 1 GB。
