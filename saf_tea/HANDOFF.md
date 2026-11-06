# 交接说明：CO₂ 制 SAF 技术经济模型（给本地 agent）

> 这份文件是写给接手的本地 agent 的。先读完本文，再读 `README.md`（结果和假设的中文说明）。

## 0. 一句话背景

用户在复现并扩展这篇论文的经济模型：
Sendlhofer 等，*Comparative techno-economic assessment of CO₂-based SAF pathways: FT, MtJ and EtJ across EU, USA and Middle East*，**Fuel 421 (2026) 139021**（预印本 SSRN 5894389）。

已完成三件事：
1. 用论文补充材料（SI）复现了论文的 4 条路线 × 3 个地区（12 个场景），误差 ±1%；
2. 加入中国地区，并新增 EtJ-2（CO₂ 一步热催化加氢制乙醇）；
3. 只对中国新增 HEFA 和 DRM-FT（甲烷干重整 + 费托），并做了中国 7 条路线的敏感性分析。

**用户有授权的 Aspen Plus**，接下来很可能要在本地用 Aspen 替换新路线里的简化计算（见第 6 节）。

## 1. 代码位置与运行

- 仓库：`github.com/jinlijin6090/Usedcardata`，分支 `claude/claude-md-mltvs6p2zb8snd6c-mCmjo`，目录 `saf_tea/`
- 依赖：Python 3.10+，`openpyxl`、`matplotlib`、`python-docx`（只有 `parse_si.py` 需要；`data/` 里已有解析结果）。中文图表字体用的是 WenQuanYi Zen Hei，本地没有的话，在 `run_all.py` 开头换成本机的中文字体（如 Microsoft YaHei、PingFang SC）。
- 运行：
  ```bash
  cd saf_tea
  python3 run_all.py      # 复现 + 5 路线 × 4 地区 + EtJ-2 / 中国 MtJ-2 敏感性 → results/
  python3 cn_routes.py    # 中国 7 路线 + 敏感性 + 交叉点（必须在 run_all.py 之后运行，会往同一个 xlsx 里追加工作表）
  ```
- 回归检查：`run_all.py` 会打印 `Replication |dev| max 1.0 %`。改代码后这个数必须保持 ≤ 1%。

## 2. 文件地图

| 文件 | 内容 |
|---|---|
| `model.py` | 核心 NPC 模型（论文式 12–16 + SI 表 31–34）；EU/US/ME 地区参数（论文表 1–3）；各路线 H₂ 流量；`npc()` 函数 |
| `extensions.py` | `CN_ASSUMPTIONS`（中国参数）、`china_region()`；`ETJ2_ASSUMPTIONS`、`etj2_balance()` / `etj2_equipment()` / `etj2_quantities()` |
| `cn_routes.py` | `HEFA_ASSUMPTIONS`、`DRM_ASSUMPTIONS`、`hefa()`、`drm_ft()`；`cn_scenario()`（中国 7 路线统一入口）；敏感性、扫描曲线、交叉点 |
| `run_all.py` | 复现对比、20 场景、美国→中国分解、EtJ-2 盈亏平衡、出图、写 Excel |
| `parse_si.py` | 把 `data/si_tables_raw.json`（SI docx 的 50 张表）转成 `data/equipment_costs.csv`、`data/varopex_si.csv` |
| `data/fig8_digitized.csv` | 论文图 8 的像素读数（复现目标） |
| `results/` | 全部输出（CSV、PNG、`SAF_TEA_results.xlsx`） |
| `README.md` | 给用户看的中文说明：假设表、结果、结论、局限 |

## 3. 复现中确认的关键设定（不要随意改）

- **运行小时与产量**：7884 h/a（由 SI 中 CO₂ 费用 ÷ CO₂ 进料量反推）；SAF 13.7 t/h → 10.8 万 t/a。正文写的"12 万吨"是名义产能，NPC 的分母用 10.8 万吨。
- **工艺段 FCI**：4.93 × EC，即 (1 + 2.51 + 0.78) × 1.15；EtJ 再乘 1.2。
- **电解 FCI**：1.70 × (800 欧元/kW × 55 kWh/kg × H₂)；只计电解槽本体，不再加地区系数。
- **资本回收**：WC = TCI 的 10%；ACC 按论文式 12；SR = floor(20 × 7884 / 60000 − 1) = 1；ARC = 0.77 × FCI_电解 / 20。
- **fixOpex 只按工艺段 FCI 计**；化验费取 **0.2 × 操作人工（OLC）**。SI 表 33 写的是 0.2 × EC，但按 EC 算对不上图 8，用户已被告知这一点。
- **EtJ 加氢补充氢**：按图 7a 电解投资反推为 267 kg/h。
- **varOpex**：EU/US/ME 的原 4 条路线直接用 SI 表 35–46 的金额。中国和新路线用"物理量 × 单价"：物理量由 EU 表 ÷ EU 单价反推（`quantities_from_eu`）。这个方法用 US/ME 的单价回代验证过，误差约 1–2%。
- **SI 里的笔误**：中东 MtJ-2 表标的是"€"，实际是 USD；美国 MtJ-2 的冷却水数值抄了欧洲。两处都不影响结果。
- **汇率**：1 USD = 0.8449 EUR = 7.1128 CNY（2025-09-19）= 3.75 SAR。

## 4. 用户已拍板的决定

- **中国的设备价格** = 美国设备清单 × 0.61（Towler & Sinnott，"China indigenous"）；电解槽单独按 450 USD/kW（设备口径）计。
- **EtJ-2 与 MtJ-2 对齐**（用户明确要求，分两次提出）：
  - 碳损失和 H₂ 吹扫损失与 MtJ-2 甲醇合成段相同（SI 表 9：碳收率 99.33%，H₂ 多耗 0.78%），不生成 CH₄：`carbon_yield = 1.0`；
  - 投资裕量 `capex_safety = 1.0`（MtJ-2 不加裕量）；
  - 反应热全部回收为中压蒸汽：`steam_recovery = 1.0`。依据是 MtJ-2 外输蒸汽约 74 MW，超过它约 48 MW 的反应热。
- **HEFA 和 DRM-FT 只算中国**。DRM-FT 默认用化石天然气价格，必须向用户注明：要算低碳 SAF，得用生物甲烷或沼气。

## 5. 当前结果（USD/kg SAF，化石 Jet A-1 ≈ 1.04）

| 地区 | FTS | MtJ-1 | MtJ-2 | EtJ | EtJ-2 | DRM-FT | HEFA |
|---|---|---|---|---|---|---|---|
| 欧洲 | 6.44 | 5.80 | 4.52 | 5.82 | 4.72 | — | — |
| 美国 | 3.49 | 3.45 | 2.68 | 3.95 | 2.89 | — | — |
| 中东 | 8.26 | 6.98 | 5.36 | 6.74 | 5.65 | — | — |
| 中国 | 2.51 | 2.45 | 1.84 | 2.60 | 1.96 | 1.75 | 2.10 |

**中国交叉点**（在这个价格以下，前一条路线比后一条便宜）：
- MtJ-2 比 HEFA 便宜：绿电 < 0.30 元/kWh
- HEFA 比 MtJ-2 便宜：地沟油 < 6,580 元/t
- DRM-FT 比 MtJ-2 便宜：甲烷 < 3.41 元/Nm³

**EtJ-2 与 MtJ-2 剩下的 0.12 差距**，几乎全部来自蒸汽和液化气的副产收入（下游乙醇脱水吸热，甲醇制丙烯放热）。

## 6. 下一步（推荐顺序）——Aspen 接入

模型里 **没有 Aspen 支撑** 的部分（风险由高到低）：

| 路线 | 目前的简化做法 | 用 Aspen 要做的 |
|---|---|---|
| EtJ-2 合成段 | 化学计量 + 等温压缩功 / 0.75；设备按 MtJ-2 甲醇合成段放缩（回路 0.6 次方，进料压缩 0.8 次方）；不建乙醇精馏 | 用用户的催化剂数据（单程转化率、选择性、副产物分布）建反应器 + 循环回路 + 闪蒸；做热集成，算真实外输蒸汽；用 APEA 出设备清单 |
| DRM-FT | 电加热重整炉，设备只是"FTS 清单 + 15 MUSD"；额外吸热 44 MW 全按用电计；碳效率沿用 FTS | 用 RGibbs 或动力学重整炉，评估 CH₄:CO₂ 比和积碳边界；FT 循环回路重新收敛；APEA 估价 |
| HEFA | 文献典型收率；FCI 按 600 USD/(t/a) 液体产品估 | 优先级较低：可以先用文献物料平衡，必要时再建 Aspen 模型 |

**接口约定**：Aspen 跑完后，把结果整理成与 SI 相同的格式，就能接入，不用改模型结构：
1. **设备清单 CSV**：列为 `route,item,EU_EUR,US_USD,ME_SAR`，追加到 `data/equipment_costs.csv`（中国用 `US_USD` × 地区系数）。
2. **年用量**（每年的物理量，不是金额）：电 MWh、净化水 m³、CO₂ t、MP/LP 蒸汽 t、轻烃 / 石脑油 / 柴油 / 蜡 kg、冷却水（每 1000 m³）、废水 m³。键名与 `model.varopex_from_quantities()` 一致。
3. **电解 H₂ 流量 kg/h**：写进 `model.H2_KGH`。
4. 然后在 `run_all.scenario()` 或 `cn_routes.cn_scenario()` 里给这条路线加一个分支，用上面三项替换现在的 `etj2_*` 或 `drm_*` 函数。

## 7. 待用户核实的估计值（交付时要提醒）

- **中国**：操作工时薪 81 元/h、CO₂ 原料 200 元/t、除盐水 15 元/m³、废水 5 元/m³、蒸汽 160–180 元/t、副产品价格、地区系数 0.61（Towler & Sinnott 原表未核实）。
- **HEFA**：收率、耗氢、投资、地沟油 7,500 元/t（没有找到权威的国内报价）。
- **DRM-FT**：重整炉造价、甲烷 3.0 元/Nm³（由 LNG 报价折算）。
- **论文本身**：不计氧气副产收入；只算投产第一年（不考虑货币时间价值）；全额贷款。

## 8. 其他产出

- 中文 16 页讲解幻灯片（论文故事版），在 claude.ai artifact：https://claude.ai/artifact/UmPepnqBNikRrHtJJnMEgR（可导出 PPTX）
- 源文件（用户上传的，不在仓库里，在打包的 `inputs/` 目录中）：
  - 正式版论文 `Fuel_2026_main.pdf`
  - 补充材料 `Fuel_2026_SI.docx`
  - 预印本 `SSRN_5894389_preprint.pdf`

## 9. 工作约定

- **语言**：用户用中文交流，交付物（README、图表、Excel）都用中文。
- **提交**：提交信息用英文；结尾加用户要求的 co-author 行（如有）。不要把 PDF/docx 原文提交进仓库。
- **新假设**：一律写进对应的 `*_ASSUMPTIONS` 字典，注明来源，标清是"估计"还是"有出处"，并同步到 README 的假设表。
