# 钴基费托：文献中的计算与实验数据

> **数据来源说明**：以下数值全部来自网络检索得到的论文摘要与搜索摘要。本环境的网络策略挡住了出版社与开放获取全文站点（pubs.acs.org、sciencedirect.com、pmc.ncbi.nlm.nih.gov、pure.tue.nl、mdpi.com、pubs.rsc.org 等），因此**没有一条经过原文核对**。标 ⚠ 的条目在摘要中的归属不够明确。引用前请逐条核对原文。

## 1. DFT 计算数据（活化能 Ea，eV）

### 1.1 CO 活化

| 步骤 | 表面 | Ea | 来源 |
|---|---|---|---|
| CO* → C* + O* | Co(0001) 平台 | 2.33（0.11 ML CO）→ 2.79（0.67 ML） | Zijlstra et al., ACS Catal. 2019, 10.1021/acscatal.9b01967 |
| CO* → C* + O* | Co(0001) | 2.40（1/4 ML）、2.26（1/9 ML），PW91 | Ge & Neurock, J. Phys. Chem. B 2006, 10.1021/jp060477i |
| CO* → C* + O* | Co(0001) | 2.82 | Inderwildi et al., J. Phys. Chem. C 2008, 10.1021/jp710674q |
| CO* → C* + O* | Co(0001) | 2.28 | Zhuo et al., J. Phys. Chem. C 2009, 10.1021/jp900281h |
| CO* → C* + O* | Co(11-21) 台阶 | ≈1.04，几乎不随 CO 覆盖度变化 ⚠ | Zijlstra et al., ACS Catal. 2020, 10.1021/acscatal.0c02420 |
| CO* → C* + O* | Co(11-20) | 2.02 | Ge & Neurock 2006 |
| CO* → C* + O* | Co(11-24) 扭折 | 0.92 | Ge & Neurock 2006 |
| CO* → C* + O* | fcc Co(321) 扭折 | 0.91 | Petersen et al., ACS Catal. 2017, 10.1021/acscatal.6b02843 |
| CO* → C* + O* | hcp Co(10-12) / (11-20) | 1.34 / 1.39 ⚠ | Qi et al.（综述），Catal. Lett. 2015 |
| CO* → C* + O* | Co(10-12) 两种终止面 | 1.55 / 1.01 | Nat. Commun. 2026, 10.1038/s41467-026-73224-8 |
| CO* + H* → HCO* | Co(0001) | 1.51 | Zhuo et al. 2009 |
| HCO* → CH* + O* | Co(0001) | 0.93（综述另给 0.73） | Zhuo et al. 2009 |
| H2CO* → CH2* + O* | Co(0001) | 0.70 | Zhuo et al. 2009 |
| HCO 路径 / HCOH 路径（总体） | Co(0001)，0.5 ML CO | 1.43 / 0.98 ⚠ | 综述归于 Ojeda / Iglesia 2010 |
| 氢助活化（总体） | Co 台阶 | 1.57 ⚠ | Zijlstra et al. 2020 |

### 1.2 CHx 加氢

| 步骤 | 平台 | 台阶 | 来源 |
|---|---|---|---|
| C* + H* → CH* | 0.41–0.85 | — | 综述（Molecules 2023, 10.3390/molecules28186525）⚠ |
| CH2* + H* → CH3* | 0.55–0.63 | 0.41–0.43 | 同上 ⚠ |
| CH3* + H* → CH4 | 0.96–1.09 | 0.88–0.96 | 同上；CH3+H 是 CHx 加氢中最慢的一步（Gong, Raval & Hu, J. Chem. Phys. 2004/05） |

### 1.3 除氧

| 步骤 | 表面 | Ea | 来源 |
|---|---|---|---|
| O* + H* → OH* | Co(0001) / 台阶 | 1.72 / 0.73 | Gong, Raval & Hu, Surf. Sci. 2004 |
| OH* + H* → H2O | Co(0001) | 1.42–1.61 ⚠ | 综述 |
| 2OH* → H2O + O* | — | 0.72（DFT）⚠；实验 0.41–0.52 | Weststrate 等 |
| O + 2H → H2O（总体） | Co(0001) | DFT 1.32，实验 1.34 ± 0.07 ⚠ | 未能确定出处 |
| O + H（实验） | Co(0001) | ≈1.24（速率控制） | J. Phys. Chem. C 2023 |

### 1.4 C–C 偶联与链终止

| 步骤 | 表面 | Ea | 来源 |
|---|---|---|---|
| CH* + CH2* | Co(0001) | 0.74 ⚠ | 可能为 Yao et al. 2019 |
| C1 + C1（有 CO 共吸附） | Co(0001) | 0.05–1.84：CH2+CH2 0.12，C+CH ≈0.70，CH3+CH3 1.84 | Catalysts 2019, 10.3390/catal9060551 |
| RCH* + CO*（CO 插入） | Co(0001) | 0.83 | Zhuo et al. 2009 |
| RCH2CHO* C–O 断裂 | Co(0001) | 0.52 | Zhuo et al. 2009 |
| 链终止（烷烃、β-H 消除）、烯烃脱附 | — | **未找到钴的数值**；Zijlstra 2020 仅称烯烃终止“略高于”链增长 | — |

### 1.5 吸附能

| 物种 | 表面 | 吸附能 | 来源 |
|---|---|---|---|
| CO | Co(0001) | 1.76 ± 0.12（PW91/PBE，高估）；实验 1.10 | 多篇；Zijlstra 2019 |
| CO | Co 台阶 | 1.24（RPBE） | Zijlstra 2019/2020 |
| H | Co(0001) | 0.97–1.1（TPD，与 DFT 一致；未注明每 H 还是每 H2） | van Helden et al., ACS Catal. 2012, 10.1021/cs2006586 |
| 乙烯 vs 更重的烯烃 | — | 乙烯强约 0.2 ⚠ | 综述片段 |

## 2. 实验数据

### 2.1 活性与选择性

| 催化剂 | Co 粒径 | T (°C) | P (bar) | H2/CO | TOF (s⁻¹) | 选择性 | 来源 |
|---|---|---|---|---|---|---|---|
| Co/CNF | 16 → 2.6 nm | 220 | 35 | 2 | 23×10⁻³ → 1.4×10⁻³ | C5+ 85 → 51 wt%；小颗粒 CH4 更高 | Bezemer et al., JACS 2006, 128, 3956 |
| Co/CNF | 2.6–27 nm | 220 | 1 | 2 | 6 nm（1 bar）/ 8 nm（35 bar）以上不变 | 6 nm 以下 CH4 增加 | 同上 |
| Co/ITQ-2 | 5.6 → 141 nm | 220 | 20 | — | 1.2×10⁻³ → 8.6×10⁻³，约 10 nm 以上不变 | — | Prieto et al., J. Catal. 2009, 266, 129 |
| 多种载体 Co | 分散度 0.01–0.12 | 200 | 20 | 2.1 | 1.6–3.0×10⁻²（与分散度无关） | C5+ > 80 % | Iglesia et al., J. Catal. 1992, 137, 212 |
| CoRu/Al2O3 | — | 200 | 20 | 2 | 2.5 ± 1.0×10⁻² | — | ACS Catal. 2021 |
| Co/TiO2 ± Ru | — | ≈200 | ≈20 | ≈2 | 加 Ru 约 ×3 | C5+ 84.5 → 91.1 % | Iglesia / Soled 1993 |
| 10 % Co/TiO2（CSTR） | — | — | — | — | — | CH4 15 → 6 %，C5+ 80 → 91 %（转化率 10 → 60 %） | Davis 组, Catal. Lett. |
| Co(-Re)/γ-Al2O3 | — | 210 | 20 | 2.1 / 1.5 / 1.0 | 随 H2/CO 降低而下降 | CH4 随 H2/CO 降低而下降 | Tristantini / Lögdberg, Fuel Process. Technol. 2007 |
| 36 种 Co 催化剂 | — | 210 | 20 | 2.1 | — | CH4 与 C5+ 选择性线性相关 | Lögdberg et al., J. Catal. 2010, 274, 84 |
| Oryx GTL（Co/Pt/Al2O3 浆态床） | — | ≈230 | 20–25 | — | — | α 0.90–0.92 | 学位论文转引 Shell / Sasol 数据 |

### 2.2 动力学

| 量 | 数值 | 条件 / 来源 |
|---|---|---|
| 速率方程 | r = a·P_CO·P_H2 / (1 + b·P_CO)² | 220–240 °C，5–15 bar — Yates & Satterfield, Energy Fuels 1991 |
| 表观活化能 | 93–95 kJ/mol | 同上 |
| 表观活化能 | 100.4 kJ/mol | Co–Re/Al2O3 — Todic et al., IECR 2013 |
| 表观活化能 | 169 → 107 kJ/mol（P_CO 2.7 → 6 bar） | Co-Pt/Al2O3 — Keyvanloo et al., Catal. Today 2016 |
| CH4 生成活化能 | ≈135.5 kJ/mol | Ma / Davis 2011 |
| H2 级数 / CO 级数 | 0.74 / −0.24 | Co/TiO2，200 °C，10 atm — Zennaro |
| H2 级数 / CO 级数 | 0.56–0.77 / +0.31（210 °C）→ −0.64（240 °C） | Keyvanloo 2016 |
| H2 级数 / CO 级数 | 0.5 / −0.25 | Co/SiO2 — Das |
| H2/D2 同位素效应 | 反向（r_H/r_D < 1） | Ojeda / Iglesia 2010 |
| SSITKA | τCHx ≈ 10 s（> 6 nm），≈ 22 s（2.6 nm） | Co/CNF，210 °C，1.85 bar，H2/CO = 10 — den Breejen et al., JACS 2009 |

### 2.3 偏离 ASF 的规律

- CH4 高于 ASF 预测；C1 这一步的 α 只有 0.52–0.57，更长链大于 0.8（ACS Catal. 2019，DTU）。
- C2 偏少，归因于乙烯快速再吸附。
- α 随链长增大，烯烃/烷烃比随链长下降到接近零（Iglesia 的烯烃再吸附模型；Kuipers 1995；van der Laan & Beenackers 1999）。
- 约 90 % 的 C3 终止生成丙烯（Rytter 2020）。

## 3. 已发表的钴基 MKM / kMC

| 研究 | 模型 | 主要结论 | 与实验的吻合 |
|---|---|---|---|
| Zijlstra et al., ACS Catal. 2020, 10.1021/acscatal.0c02420 | 台阶 Co(11-21) + 平台 Co(0001) 双位点，含侧向相互作用 | 直接 CO 解离与碳化物链增长几乎只在台阶上发生；平台通过 CH/CH3 迁移增加 CH4；速率由除氧和 CO 解离控制 | α ≈ 0.4–0.5，低于真实催化剂 |
| Rommens, Gunasooriya & Saeys, ACS Catal. 2024, 10.1021/acscatal.3c04844 | 600 多个基元步骤，DFT 在 CO 饱和的平台上计算 | 链增长为 CH 插入次烷基 | CH4 18 %，α 0.83，TOF 0.084 s⁻¹，Ea 107 kJ/mol；改用低覆盖度 DFT 时 CH4 升到 43 % |
| Yao, Guo, Mao & Hu, ACS Catal. 2019, 10.1021/acscatal.9b01150 | Co(0001)，含覆盖度效应 | 低 CHx 覆盖度下 CO 插入占优，高覆盖度下碳化物占优 | 计入覆盖度效应后 TOF 提高约 10⁶ 倍；Co(0001) 偏向生成烯烃 |
| Chen et al., ACS Catal. 2017；Pestman et al., ACS Catal. 2019 | 结合实验的动力学分析 | 速率由 C 加氢或加氢终止控制，而非 C–O 断裂；直接 CO 解离比氢助路径更符合数据 | — |
| Inderwildi 2008；Storsæter 2006；Yang / Holmen 2013 | DFT / UBI-QEP | 甲酰基（氢助）路径占优；CO 插入为链增长 | — |
| Su, Zhao, Liu & Li, Catal. Sci. Technol. 2017 | DFT | Co(0001) 和台阶上为 CO 插入，Co(10-11) 上为碳化物机理 | — |
