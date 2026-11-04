# EC 化学空间硬枚举（ec_space）

这个程序从碳酸乙烯酯（EC，O=C1OCCO1）出发，枚举整个 C/H/O(/Li) 化学空间：
- 闭壳层分子，以及所有自由基、离子、Li⁺ 配位物；
- 各物种之间的基元反应；
- 结果写成可以用 SQL 查询的数据库，并生成一个本地网页，可以直接浏览、查看每个物种的全部反应，以及从 EC 出发的最短路径。

程序自带完整的 Python 环境（RDKit、numpy、DuckDB），不依赖服务器上的 Python 或 conda，也不需要联网。它会用满分配到的 CPU，中断后可以从断点继续。

---

## 1. 三步运行

1. **放好文件**
   - 把 `ec_space/` 整个目录拷到服务器上；
   - 把环境包放进 `ec_space/env/`：可以是 `ecenv.tar.gz`，也可以是分卷 `ecenv.tar.gz.part-00 … part-03`，分卷会自动合并。
2. **（可选）先自检，约 1 分钟**
   ```bash
   bash ec_space/run.sh --selftest
   ```
   第一次运行会把环境解压到 `ec_space/env/ecenv/`，只做一次。
3. **提交计算**：在你自己的提交文件里只写这一行：
   ```bash
   bash /path/to/ec_space/run.sh /path/to/ec_space/config.ini
   ```
   `submit_example.slurm` 和 `submit_example.pbs` 是 SLURM 和 PBS 的示例。程序会自动读取作业分配到的核数（SLURM / PBS / LSF / SGE），否则使用本机全部核。

**Windows / macOS 或已有 Python 环境**：先按仓库根目录的 `environment.yml`（或 `pip install rdkit==2026.9.1 numpy duckdb`）装好依赖，再用跨平台启动器：
```bash
python ec_space/run.py --selftest
python ec_space/run.py ec_space/config.ini
```
多进程在 Linux 上用 fork，在 Windows / macOS 上自动改用 spawn，两种方式结果一致。

想改规模，可以改 `config.ini`，也可以直接在命令行覆盖：
```bash
bash run.sh config.ini --set space.max_heavy=10 --set write.derive_max_heavy=8
bash run.sh config.ini --stages estimate        # 只看规模估算，不真正计算
```

## 2. 运行流程

程序依次执行下面六个阶段，每个阶段都可以断点续跑：

| 阶段 | 做什么 | 输出 |
|---|---|---|
| estimate | 在本机精确枚举前 7 层，外推出各层的物种数、反应数、CPU 和磁盘需求 | `ESTIMATE.md` |
| enumerate | 逐层生成全部闭壳层分子（每次加一个 C 或 O），在磁盘上按哈希分片去重 | `closed/L{n}/s*.smi.gz` |
| derive | 对每个闭壳层分子生成它派生的全部自由基、离子、Li⁺ 配位物和基元反应 | `species/`、`reactions/`、`derive/L{n}.json` |
| solvent | 每个自由基 / 负离子与溶剂分子（`[solvent] solvents`）的双分子反应，见第 7 节 | `solvent/<TAG>/` |
| stats | 汇总统计 | `stats.json`、`SUMMARY_<TAG>.md`（中文表格） |
| viewer | 生成本地网页 | `viewer_<TAG>/index.html`，双击即可打开，不需要联网 |

- **断点续跑**：每完成一层就写一个 `DONE.json`。作业中断或超时后，用同一条命令重新提交，已完成的层会被跳过。
- **超大层只计数**：重原子数大于 `derive_max_heavy` 的层，派生物种和反应不逐条写出，只按 `count_sample` 抽样统计数量。闭壳层分子始终全部写出。

## 3. 枚举了什么

### 3.1 闭壳层分子（硬枚举）

- 元素：C、O、H，价态 C 4、O 2，H 隐含；所有键级（单、双、三键）；任意环、任意大小的环、芳香环。
- 每次加一个重原子，接到已有原子的任意子集上。这样能遍历所有连通分子，因为任何分子都能删掉一个"不切断分子"的原子退回上一层。结果用 RDKit 规范 SMILES 去重（经过往返校验，防止同一分子出现两种写法）。
- **只去掉以下情况：**
  - 8 元以下的环里出现三键或累积双键（几何上无法成环）；
  - O–O–O 链；`allow_peroxide = no` 时不允许任何 O–O。
  - 单个过氧键 O–O 默认保留。
  - 可选 `strain_filter = plausible`：再去掉同一原子处在两个 ≤4 元环里、桥环（含反 Bredt 烯烃）、丁三烯及更长的累积多烯。
- 另外手动加入的特殊物种：CO（`[C-]#[O+]`），因为它不能用常规价态写出。
- **正确性**：关闭过氧键时，1–9 个重原子的分子数（2、5、17、77、379、2 270、15 196、115 323、964 497）与另一套独立实现的枚举结果完全一致（第 8 层的独立结果是 115 325，多出的 2 个是同一分子被写成两种 SMILES，用往返校验修正后为 115 323）。

### 3.2 派生物种

每个物种都能把去掉的 H 或电子放回去，唯一地还原出它的"母体"分子，所以不同母体之间不会产生重复：

| 类别 | 构造方式 | 配置项 |
|---|---|---|
| 自由基 R• | 去掉一个 H•（C 中心或 O 中心） | `radicals` |
| 负离子 R⁻ | 去掉一个 H⁺（碳负离子、烷氧负离子） | `anions` |
| 正离子 R⁺ | 去掉一个 H⁻，默认只取碳正离子 | `cation_centers` |
| 质子化 MH⁺ | 在 O 上加 H⁺（氧鎓离子） | `protonation` |
| 自由基阴 / 阳离子 M•⁻ / M•⁺ | 整体加 / 失一个电子，不定域，写作 `母体SMILES\|-1` / `\|+1` | `radical_ions` |
| 双位点物种 | 远端自由基阴离子（如 EC 开环后的 •CH₂CH₂OCO₂⁻）、O 中心二价阴离子，以及可选的双自由基 / 卡宾、两性离子、远端自由基阳离子 | `pairs`、`carbenes` |
| Li⁺ 配位物 | 以上任一物种的 O 上配位 1–2 个 Li⁺（O→Li⁺ 配位键）。负离子 + Li⁺ 即锂盐，例如 LEDC 由二价阴离子和 2 个 Li⁺ 组成 | `li_max` |

只保留总电荷在 `charges`（默认 −2、−1、0、+1）范围内的物种。

### 3.3 基元反应

每条反应只从一侧生成一次，不需要全局去重。反应物或产物超出配置范围的反应不写出。

| 反应族 | 形式 |
|---|---|
| `h_loss` / `deprotonation` / `hydride_loss` / `protonation` | M ⇌ R• + H•；M ⇌ R⁻ + H⁺；M ⇌ R⁺ + H⁻；M + H⁺ ⇌ MH⁺ |
| `electron_capture` / `electron_loss` / `site_electron_transfer` | M + e⁻ ⇌ M•⁻；M ⇌ M•⁺ + e⁻；R• + e⁻ ⇌ R⁻；R• ⇌ R⁺ + e⁻ |
| `li_binding` | X + Li⁺ ⇌ X·Li⁺（每次一个 Li） |
| `homolysis` / `heterolysis` | 每根单键：A–B ⇌ A• + B•；A–B ⇌ A⁺ + B⁻（两种极性）。逆反应即自由基复合 / 离子复合 |
| `beta_scission` | 自由基 β-断裂（含开环）。逆反应即自由基加成 |
| `radical_anion_fragmentation` / `radical_cation_fragmentation` | M•⁻ ⇌ A• + B⁻，M•⁺ ⇌ A• + B⁺。环上的键断开即开环，例如 EC•⁻ → •CH₂CH₂OCO₂⁻ |
| `pair_steps` | 双位点物种由两个单位点物种各走一步形成 |

**不枚举的反应**：任意两个"大"物种之间的双分子反应（例如 R• + R′H 的 H 抽提）。这类反应数量约为物种数的平方，在本空间里是 10²⁰ 量级，没法列出来。自由基复合、离子复合、自由基加成已经作为上面断键反应的逆反应包含在内。

## 4. 规模与资源

下表是默认配置（`max_heavy = 12`，`derive_max_heavy = 9`）用 `--stages estimate` 得到的估算。前 7 层是精确值，其余按实测增长外推，误差约 ±50%。

| 重原子 | 闭壳层分子 | 全部物种 | 基元反应 | 逐条写出 |
|---|---|---|---|---|
| 7 | 17 851 | 1.36×10⁶ | 2.6×10⁶ | 是 |
| 8 | 1.3×10⁵ | 1.3×10⁷ | 2.5×10⁷ | 是 |
| 9 | 1.1×10⁶ | 1.2×10⁸ | 2.5×10⁸ | 是 |
| 10 | 9×10⁶ | 1.2×10⁹ | 2.6×10⁹ | 抽样计数 |
| 11 | 9×10⁷ | 1.3×10¹⁰ | 2.9×10¹⁰ | 抽样计数 |
| 12 | 9×10⁸ | 1.5×10¹¹ | 3.3×10¹¹ | 抽样计数 |

默认配置的资源需求：

| 项目 | 需求 |
|---|---|
| CPU | 约 1 400 核时（64 核约 22 小时） |
| 写出数据 | 约 8 GB |
| 临时文件峰值 | 约 26 GB |
| 内存 | 建议每核 2 GB |

如果要把派生物种和反应逐条写到更大的层，按下表调整 `derive_max_heavy`：

| `derive_max_heavy` | 新增写出 | 新增 CPU |
|---|---|---|
| 10 | 约 +35 GB | 约 +200 核时 |
| 11 | 约 +380 GB | 约 +2 500 核时 |
| 12 | 约 4 TB | 约 3 万核时 |

实测的压缩后大小：每个物种约 7 字节，每条反应约 8.4 字节，每个闭壳层分子约 4.5 字节。

## 5. 数据库与查询

所有文件都是 gzip 压缩的 TSV，没有表头行（列名在 `species/HEADER.tsv` 和 `reactions/HEADER.tsv` 里），可以直接拼接：

- `species/L{n}/part-*.tsv.gz`：`key  type  charge  multiplicity  C  H  O  Li  heavy  parent`
- `reactions/L{n}/part-*.tsv.gz`：`family  reactants  products  parent_heavy`
  - 多个反应物或产物之间用 ` + ` 分隔；`[e-]`、`[Li+]`、`[H]`、`[H+]`、`[H-]` 是池物种，定义在 `species/POOL.tsv`。
- `closed/L{n}/s*.smi.gz`：全部闭壳层分子，包括只计数的大层。

用自带的 DuckDB 直接写 SQL 查询（不需要先导入）：
```bash
P=ec_space/env/ecenv/bin/python; export PYTHONPATH=ec_space
$P -m ecspace.query ec_space_out "SELECT type, count(*) FROM species GROUP BY 1 ORDER BY 2 DESC"
$P -m ecspace.query ec_space_out --reactions-of "O=C1OCCO1"
$P -m ecspace.query ec_space_out --export ledc.parquet "SELECT * FROM species WHERE Li = 2 AND C = 4 AND O = 6"
```

## 6. 网页预览

`ec_space_out/viewer_<TAG>/index.html`（例如 `viewer_EC/`；不设溶剂时为 `viewer/`）双击就能打开，不需要联网，也可以拷到别的电脑上看。
- 列是重原子数，行是物种类别。放大到哪里，就读取哪里的标签和反应。
- 单击物种可以查看结构式、它参与的全部基元反应，以及从 EC 出发的最短路径（e⁻、Li⁺、H•、H⁺、H⁻ 视为始终可用）。
- 例如 EC + e⁻ → EC•⁻ → •CH₂CH₂OCO₂⁻ 是 2 步。
- 网页默认只包含重原子 ≤ 7 的物种（约 140 万个，浏览器里可以流畅加载），可在 `[viewer] max_heavy` 调整。每格物种过多时按 `max_per_cell` 抽样。

## 7. 其他溶剂与混合溶剂：EC、DMC、EA、DME、EC/DMC

`presets/` 里有 5 个预设，详细说明见 `presets/README.md`：

| 预设 | 溶剂 | 起点 |
|---|---|---|
| `presets/EC.ini` | EC | EC |
| `presets/DMC.ini` | DMC | DMC |
| `presets/EA.ini` | EA | EA |
| `presets/DME.ini` | DME（1,2-二甲氧基乙烷） | DME |
| `presets/EC_DMC.ini` | EC + DMC | EC 和 DMC |

要点：
- EC、DMC、EA、DME 都是 6 个重原子，"两个分子"都是 ≤ 12 个重原子。闭壳层分子和全部自由基、离子、Li⁺ 配位物的空间**与溶剂无关**，只算一次即可。
- 随溶剂变化的是新增的 `solvent` 阶段：每个定域、不含 Li 的自由基 / 负离子，与每种溶剂分子之间的双分子反应。共四类：
  - 从溶剂抽 H；
  - 夺溶剂质子；
  - 亲核酰基取代（酯交换、EC 开环、Claisen 型 C-酰化）；
  - SN2 烷基转移（DMC 甲基化、DME 醚断裂）。

  混合溶剂时 EC↔DMC 之间的交叉反应自动包括在内，互逆的反应只记一次。
- **操作**：
  1. 第一次 `bash run.sh presets/EC.ini`，完整计算；
  2. 之后每换一种溶剂，只跑 `bash run.sh presets/DMC.ini --stages solvent,stats,viewer`，共用同一个 `ec_space_out/`。
  3. 结果分别在 `solvent/<TAG>/`、`SUMMARY_<TAG>.md` 和 `viewer_<TAG>/`。

规模（重原子 ≤ 8 的测试中实测：溶剂反应数与基础反应数之比各层基本不变，下表按此比例外推到默认配置）：

| 体系 | 溶剂反应 / 基础反应 | ≤ 12 重原子的溶剂反应（外推） | 写出（≤ 9 层） | CPU |
|---|---|---|---|---|
| EC、DMC | 0.22 | 约 8×10¹⁰ | 约 0.5 GB | 每种溶剂约 50–150 核时（实测为派生阶段的 0.3–0.4 倍），远小于基础空间的约 1 400 核时 |
| DME、EC+DMC | 0.43 | 约 1.6×10¹¹ | 约 0.9 GB | 同上 |
| EA | 0.65 | 约 2.4×10¹¹ | 约 1.4 GB | 同上 |

## 8. 文件

```
ec_space/
  run.sh                 一键运行（首次自动解压环境）
  config.ini             配置（中文注释）
  submit_example.slurm   SLURM 提交示例
  submit_example.pbs     PBS 提交示例
  ecspace/               程序（chem.py 化学规则，pipeline.py 并行枚举与去重，solvent.py 溶剂双分子反应，viewer.py 网页数据，query.py 查询）
  viewer/                网页模板与 smiles-drawer（MIT 许可证）
  presets/               溶剂预设：EC、DMC、EA、DME、EC_DMC（说明见 presets/README.md）
  tests/selftest.sh      自检（含 EC/DMC 溶剂反应检查）
  env/                   放环境包 ecenv.tar.gz 或分卷 ecenv.tar.gz.part-*（Python 3.11、RDKit 2026.09.1、numpy 2.4、DuckDB 2.0；conda-pack 打包，Linux x86-64）
```
