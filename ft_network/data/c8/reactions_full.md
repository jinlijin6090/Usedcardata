# 费托合成直链烃基元反应网络全表 (C1–C8)

由 `scripts/build_network.py --nmax 8` 自动生成；共 258 个物种、578 个可逆基元步骤，每一步均通过 C/H/O 原子守恒与表面位点守恒校验。`*` 为表面空位。


## 通用步骤（与碳数无关）

### 1 反应物吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0001 | `CO(g) + * ⇌ CO*` | CO 分子吸附 |
| R0002 | `H2(g) + 2 * ⇌ 2 H*` | H2 解离吸附 |

## C1 化学

### 2a CO 直接解离

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0003 | `CO* + * ⇌ C* + O*` | 碳化物机理的起点 |
| R0004 | `2 CO* ⇌ C* + CO2*` | Boudouard 歧化 |

### 2b CO 氢助活化(加氢)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0005 | `CO* + H* ⇌ HCO* + *` |  |
| R0006 | `CO* + H* ⇌ COH* + *` |  |
| R0007 | `HCO* + H* ⇌ CH2O* + *` |  |
| R0008 | `HCO* + H* ⇌ HCOH* + *` |  |
| R0009 | `COH* + H* ⇌ HCOH* + *` |  |
| R0010 | `CH2O* + H* ⇌ CH3O* + *` |  |
| R0011 | `CH2O* + H* ⇌ CH2OH* + *` |  |
| R0012 | `HCOH* + H* ⇌ CH2OH* + *` |  |

### 2c 氢助 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0013 | `HCO* + * ⇌ CH* + O*` |  |
| R0014 | `COH* + * ⇌ C* + OH*` |  |
| R0015 | `HCOH* + * ⇌ CH* + OH*` |  |
| R0016 | `CH2O* + * ⇌ CH2* + O*` |  |
| R0017 | `CH2OH* + * ⇌ CH2* + OH*` |  |
| R0018 | `CH3O* + * ⇌ CH3* + O*` |  |

## 通用步骤（与碳数无关）

### 3a 除氧生成 H2O

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0019 | `O* + H* ⇌ OH* + *` |  |
| R0020 | `OH* + H* ⇌ H2O* + *` |  |
| R0021 | `2 OH* ⇌ H2O* + O*` | 羟基歧化 |
| R0022 | `H2O* ⇌ H2O(g) + *` | H2O 脱附 |

## C1 化学

### 3b 除氧生成 CO2 / 水煤气变换

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0023 | `CO* + O* ⇌ CO2* + *` |  |
| R0024 | `CO* + OH* ⇌ COOH* + *` | 羧基路径 |
| R0025 | `COOH* + * ⇌ CO2* + H*` |  |
| R0026 | `COOH* + O* ⇌ CO2* + OH*` |  |
| R0027 | `COOH* + OH* ⇌ CO2* + H2O*` |  |
| R0028 | `HCO* + O* ⇌ HCOO* + *` | 甲酸根路径 |
| R0029 | `HCOO* + * ⇌ CO2* + H*` |  |
| R0030 | `CO2* ⇌ CO2(g) + *` | CO2 脱附 |

### 4 C1 加氢(链引发单体 CHx 生成)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0031 | `C* + H* ⇌ CH* + *` |  |
| R0032 | `CH* + H* ⇌ CH2* + *` |  |
| R0033 | `CH2* + H* ⇌ CH3* + *` | 链引发物种 CH3* |
| R0034 | `CH3* + H* ⇌ CH4* + *` | 甲烷化(链终止 n=1) |
| R0035 | `CH4* ⇌ CH4(g) + *` | CH4 脱附 |

## C1 → C2

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0036 | `2 C* ⇌ CC* + *` | 与表面 C 偶联 |
| R0037 | `C* + CH* ⇌ CHC* + *` |  |
| R0038 | `C* + CH2* ⇌ CH2C* + *` |  |
| R0039 | `C* + CH3* ⇌ CH3C* + *` |  |
| R0040 | `2 CH* ⇌ CHCH* + *` |  |
| R0041 | `CH* + CH2* ⇌ CH2CH* + *` |  |
| R0042 | `CH* + CH3* ⇌ CH3CH* + *` |  |
| R0043 | `2 CH2* ⇌ CH2CH2* + *` |  |
| R0044 | `CH2* + CH3* ⇌ CH3CH2* + *` |  |
| R0045 | `2 CH3* ⇌ CH3CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0046 | `C* + CO* ⇌ CCO* + *` |  |
| R0047 | `CH* + CO* ⇌ CHCO* + *` |  |
| R0048 | `CH2* + CO* ⇌ CH2CO* + *` |  |
| R0049 | `CH3* + CO* ⇌ CH3CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0050 | `CCO* + H* ⇌ CCHO* + *` | C 端加氢 |
| R0051 | `CCHO* + H* ⇌ CCH2O* + *` | C 端加氢 |
| R0052 | `CCOH* + H* ⇌ CCHOH* + *` | C 端加氢 |
| R0053 | `CCHOH* + H* ⇌ CCH2OH* + *` | C 端加氢 |
| R0054 | `CCO* + H* ⇌ CCOH* + *` | O 端加氢 |
| R0055 | `CCHO* + H* ⇌ CCHOH* + *` | O 端加氢 |
| R0056 | `CCH2O* + H* ⇌ CCH2OH* + *` | O 端加氢 |
| R0057 | `CHCO* + H* ⇌ CHCHO* + *` | C 端加氢 |
| R0058 | `CHCHO* + H* ⇌ CHCH2O* + *` | C 端加氢 |
| R0059 | `CHCOH* + H* ⇌ CHCHOH* + *` | C 端加氢 |
| R0060 | `CHCHOH* + H* ⇌ CHCH2OH* + *` | C 端加氢 |
| R0061 | `CHCO* + H* ⇌ CHCOH* + *` | O 端加氢 |
| R0062 | `CHCHO* + H* ⇌ CHCHOH* + *` | O 端加氢 |
| R0063 | `CHCH2O* + H* ⇌ CHCH2OH* + *` | O 端加氢 |
| R0064 | `CH2CO* + H* ⇌ CH2CHO* + *` | C 端加氢 |
| R0065 | `CH2CHO* + H* ⇌ CH2CH2O* + *` | C 端加氢 |
| R0066 | `CH2COH* + H* ⇌ CH2CHOH* + *` | C 端加氢 |
| R0067 | `CH2CHOH* + H* ⇌ CH2CH2OH* + *` | C 端加氢 |
| R0068 | `CH2CO* + H* ⇌ CH2COH* + *` | O 端加氢 |
| R0069 | `CH2CHO* + H* ⇌ CH2CHOH* + *` | O 端加氢 |
| R0070 | `CH2CH2O* + H* ⇌ CH2CH2OH* + *` | O 端加氢 |
| R0071 | `CH3CO* + H* ⇌ CH3CHO* + *` | C 端加氢 |
| R0072 | `CH3CHO* + H* ⇌ CH3CH2O* + *` | C 端加氢 |
| R0073 | `CH3COH* + H* ⇌ CH3CHOH* + *` | C 端加氢 |
| R0074 | `CH3CHOH* + H* ⇌ CH3CH2OH* + *` | C 端加氢 |
| R0075 | `CH3CO* + H* ⇌ CH3COH* + *` | O 端加氢 |
| R0076 | `CH3CHO* + H* ⇌ CH3CHOH* + *` | O 端加氢 |
| R0077 | `CH3CH2O* + H* ⇌ CH3CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0078 | `CCO* + * ⇌ CC* + O*` |  |
| R0079 | `CCOH* + * ⇌ CC* + OH*` |  |
| R0080 | `CCHO* + * ⇌ CHC* + O*` |  |
| R0081 | `CCHOH* + * ⇌ CHC* + OH*` |  |
| R0082 | `CCH2O* + * ⇌ CH2C* + O*` |  |
| R0083 | `CCH2OH* + * ⇌ CH2C* + OH*` |  |
| R0084 | `CHCO* + * ⇌ CHC* + O*` |  |
| R0085 | `CHCOH* + * ⇌ CHC* + OH*` |  |
| R0086 | `CHCHO* + * ⇌ CHCH* + O*` |  |
| R0087 | `CHCHOH* + * ⇌ CHCH* + OH*` |  |
| R0088 | `CHCH2O* + * ⇌ CH2CH* + O*` |  |
| R0089 | `CHCH2OH* + * ⇌ CH2CH* + OH*` |  |
| R0090 | `CH2CO* + * ⇌ CH2C* + O*` |  |
| R0091 | `CH2COH* + * ⇌ CH2C* + OH*` |  |
| R0092 | `CH2CHO* + * ⇌ CH2CH* + O*` |  |
| R0093 | `CH2CHOH* + * ⇌ CH2CH* + OH*` |  |
| R0094 | `CH2CH2O* + * ⇌ CH2CH2* + O*` |  |
| R0095 | `CH2CH2OH* + * ⇌ CH2CH2* + OH*` |  |
| R0096 | `CH3CO* + * ⇌ CH3C* + O*` |  |
| R0097 | `CH3COH* + * ⇌ CH3C* + OH*` |  |
| R0098 | `CH3CHO* + * ⇌ CH3CH* + O*` |  |
| R0099 | `CH3CHOH* + * ⇌ CH3CH* + OH*` |  |
| R0100 | `CH3CH2O* + * ⇌ CH3CH2* + O*` |  |
| R0101 | `CH3CH2OH* + * ⇌ CH3CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0102 | `2 HCOH* ⇌ CHCOH* + H2O*` | 两个羟基卡宾脱水缩合 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0103 | `CC* + H* ⇌ CHC* + *` |  |
| R0104 | `CHC* + H* ⇌ CH2C* + *` |  |
| R0105 | `CHC* + H* ⇌ CHCH* + *` |  |
| R0106 | `CHCH* + H* ⇌ CH2CH* + *` |  |
| R0107 | `CH2C* + H* ⇌ CH3C* + *` |  |
| R0108 | `CH2C* + H* ⇌ CH2CH* + *` |  |
| R0109 | `CH2CH* + H* ⇌ CH3CH* + *` |  |
| R0110 | `CH2CH* + H* ⇌ CH2CH2* + *` |  |
| R0111 | `CH2CH2* + H* ⇌ CH3CH2* + *` | 逆反应 = 乙基 β-H 消除 → 乙烯 |
| R0112 | `CH3C* + H* ⇌ CH3CH* + *` |  |
| R0113 | `CH3CH* + H* ⇌ CH3CH2* + *` |  |
| R0114 | `CH3CH2* + H* ⇌ CH3CH3* + *` | 链终止: 乙基加氢 → 乙烷 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0115 | `CH3CH3* ⇌ C2H6(g) + *` | 乙烷脱附 |
| R0116 | `CH2CH2* ⇌ C2H4(g) + *` | 乙烯脱附; 逆向 = 再吸附 |

## C2 → C3

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0117 | `CH3C* + C* ⇌ CH3CC* + *` | 与表面 C 偶联 |
| R0118 | `CH3C* + CH* ⇌ CH3CCH* + *` | 次烷基 + CH |
| R0119 | `CH3C* + CH2* ⇌ CH3CCH2* + *` |  |
| R0120 | `CH3C* + CH3* ⇌ CH3CCH3* + *` |  |
| R0121 | `CH3CH* + C* ⇌ CH3CHC* + *` | 与表面 C 偶联 |
| R0122 | `CH3CH* + CH* ⇌ CH3CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0123 | `CH3CH* + CH2* ⇌ CH3CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0124 | `CH3CH* + CH3* ⇌ CH3CHCH3* + *` |  |
| R0125 | `CH3CH2* + C* ⇌ CH3CH2C* + *` | 与表面 C 偶联 |
| R0126 | `CH3CH2* + CH* ⇌ CH3CH2CH* + *` |  |
| R0127 | `CH3CH2* + CH2* ⇌ CH3CH2CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0128 | `CH3CH2* + CH3* ⇌ CH3CH2CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0129 | `CH3C* + CO* ⇌ CH3CCO* + *` |  |
| R0130 | `CH3CH* + CO* ⇌ CH3CHCO* + *` |  |
| R0131 | `CH3CH2* + CO* ⇌ CH3CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0132 | `CH3CCO* + H* ⇌ CH3CCHO* + *` | C 端加氢 |
| R0133 | `CH3CCHO* + H* ⇌ CH3CCH2O* + *` | C 端加氢 |
| R0134 | `CH3CCOH* + H* ⇌ CH3CCHOH* + *` | C 端加氢 |
| R0135 | `CH3CCHOH* + H* ⇌ CH3CCH2OH* + *` | C 端加氢 |
| R0136 | `CH3CCO* + H* ⇌ CH3CCOH* + *` | O 端加氢 |
| R0137 | `CH3CCHO* + H* ⇌ CH3CCHOH* + *` | O 端加氢 |
| R0138 | `CH3CCH2O* + H* ⇌ CH3CCH2OH* + *` | O 端加氢 |
| R0139 | `CH3CHCO* + H* ⇌ CH3CHCHO* + *` | C 端加氢 |
| R0140 | `CH3CHCHO* + H* ⇌ CH3CHCH2O* + *` | C 端加氢 |
| R0141 | `CH3CHCOH* + H* ⇌ CH3CHCHOH* + *` | C 端加氢 |
| R0142 | `CH3CHCHOH* + H* ⇌ CH3CHCH2OH* + *` | C 端加氢 |
| R0143 | `CH3CHCO* + H* ⇌ CH3CHCOH* + *` | O 端加氢 |
| R0144 | `CH3CHCHO* + H* ⇌ CH3CHCHOH* + *` | O 端加氢 |
| R0145 | `CH3CHCH2O* + H* ⇌ CH3CHCH2OH* + *` | O 端加氢 |
| R0146 | `CH3CH2CO* + H* ⇌ CH3CH2CHO* + *` | C 端加氢 |
| R0147 | `CH3CH2CHO* + H* ⇌ CH3CH2CH2O* + *` | C 端加氢 |
| R0148 | `CH3CH2COH* + H* ⇌ CH3CH2CHOH* + *` | C 端加氢 |
| R0149 | `CH3CH2CHOH* + H* ⇌ CH3CH2CH2OH* + *` | C 端加氢 |
| R0150 | `CH3CH2CO* + H* ⇌ CH3CH2COH* + *` | O 端加氢 |
| R0151 | `CH3CH2CHO* + H* ⇌ CH3CH2CHOH* + *` | O 端加氢 |
| R0152 | `CH3CH2CH2O* + H* ⇌ CH3CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0153 | `CH3CCO* + * ⇌ CH3CC* + O*` |  |
| R0154 | `CH3CCOH* + * ⇌ CH3CC* + OH*` |  |
| R0155 | `CH3CCHO* + * ⇌ CH3CCH* + O*` |  |
| R0156 | `CH3CCHOH* + * ⇌ CH3CCH* + OH*` |  |
| R0157 | `CH3CCH2O* + * ⇌ CH3CCH2* + O*` |  |
| R0158 | `CH3CCH2OH* + * ⇌ CH3CCH2* + OH*` |  |
| R0159 | `CH3CHCO* + * ⇌ CH3CHC* + O*` |  |
| R0160 | `CH3CHCOH* + * ⇌ CH3CHC* + OH*` |  |
| R0161 | `CH3CHCHO* + * ⇌ CH3CHCH* + O*` |  |
| R0162 | `CH3CHCHOH* + * ⇌ CH3CHCH* + OH*` |  |
| R0163 | `CH3CHCH2O* + * ⇌ CH3CHCH2* + O*` |  |
| R0164 | `CH3CHCH2OH* + * ⇌ CH3CHCH2* + OH*` |  |
| R0165 | `CH3CH2CO* + * ⇌ CH3CH2C* + O*` |  |
| R0166 | `CH3CH2COH* + * ⇌ CH3CH2C* + OH*` |  |
| R0167 | `CH3CH2CHO* + * ⇌ CH3CH2CH* + O*` |  |
| R0168 | `CH3CH2CHOH* + * ⇌ CH3CH2CH* + OH*` |  |
| R0169 | `CH3CH2CH2O* + * ⇌ CH3CH2CH2* + O*` |  |
| R0170 | `CH3CH2CH2OH* + * ⇌ CH3CH2CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0171 | `CH3COH* + HCOH* ⇌ CH3CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0172 | `CH2CH* + CH2* ⇌ CH2CHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0173 | `CH2CHCH2* ⇌ CH3CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0174 | `CH2CHCH2* + H* ⇌ CH3CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0175 | `CH3CC* + H* ⇌ CH3CCH* + *` | α-C 加氢 |
| R0176 | `CH3CC* + H* ⇌ CH3CHC* + *` | β-C 加氢 |
| R0177 | `CH3CCH* + H* ⇌ CH3CCH2* + *` | α-C 加氢 |
| R0178 | `CH3CCH* + H* ⇌ CH3CHCH* + *` | β-C 加氢 |
| R0179 | `CH3CCH2* + H* ⇌ CH3CCH3* + *` | α-C 加氢 |
| R0180 | `CH3CCH2* + H* ⇌ CH3CHCH2* + *` | β-C 加氢 |
| R0181 | `CH3CCH3* + H* ⇌ CH3CHCH3* + *` | β-C 加氢 |
| R0182 | `CH3CHC* + H* ⇌ CH3CHCH* + *` | α-C 加氢 |
| R0183 | `CH3CHC* + H* ⇌ CH3CH2C* + *` | β-C 加氢 |
| R0184 | `CH3CHCH* + H* ⇌ CH3CHCH2* + *` | α-C 加氢 |
| R0185 | `CH3CHCH* + H* ⇌ CH3CH2CH* + *` | β-C 加氢 |
| R0186 | `CH3CHCH2* + H* ⇌ CH3CHCH3* + *` | α-C 加氢 |
| R0187 | `CH3CHCH2* + H* ⇌ CH3CH2CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0188 | `CH3CHCH3* + H* ⇌ CH3CH2CH3* + *` | β-C 加氢 |
| R0189 | `CH3CH2C* + H* ⇌ CH3CH2CH* + *` | α-C 加氢 |
| R0190 | `CH3CH2CH* + H* ⇌ CH3CH2CH2* + *` | α-C 加氢 |
| R0191 | `CH3CH2CH2* + H* ⇌ CH3CH2CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0192 | `CH3CH2CH3* ⇌ C3H8(g) + *` | 正构烷烃脱附 |
| R0193 | `CH3CHCH2* ⇌ C3H6(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C3 → C4

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0194 | `CH3CH2C* + C* ⇌ CH3CH2CC* + *` | 与表面 C 偶联 |
| R0195 | `CH3CH2C* + CH* ⇌ CH3CH2CCH* + *` | 次烷基 + CH |
| R0196 | `CH3CH2C* + CH2* ⇌ CH3CH2CCH2* + *` |  |
| R0197 | `CH3CH2C* + CH3* ⇌ CH3CH2CCH3* + *` |  |
| R0198 | `CH3CH2CH* + C* ⇌ CH3CH2CHC* + *` | 与表面 C 偶联 |
| R0199 | `CH3CH2CH* + CH* ⇌ CH3CH2CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0200 | `CH3CH2CH* + CH2* ⇌ CH3CH2CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0201 | `CH3CH2CH* + CH3* ⇌ CH3CH2CHCH3* + *` |  |
| R0202 | `CH3CH2CH2* + C* ⇌ CH3(CH2)2C* + *` | 与表面 C 偶联 |
| R0203 | `CH3CH2CH2* + CH* ⇌ CH3(CH2)2CH* + *` |  |
| R0204 | `CH3CH2CH2* + CH2* ⇌ CH3(CH2)2CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0205 | `CH3CH2CH2* + CH3* ⇌ CH3(CH2)2CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0206 | `CH3CH2C* + CO* ⇌ CH3CH2CCO* + *` |  |
| R0207 | `CH3CH2CH* + CO* ⇌ CH3CH2CHCO* + *` |  |
| R0208 | `CH3CH2CH2* + CO* ⇌ CH3CH2CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0209 | `CH3CH2CCO* + H* ⇌ CH3CH2CCHO* + *` | C 端加氢 |
| R0210 | `CH3CH2CCHO* + H* ⇌ CH3CH2CCH2O* + *` | C 端加氢 |
| R0211 | `CH3CH2CCOH* + H* ⇌ CH3CH2CCHOH* + *` | C 端加氢 |
| R0212 | `CH3CH2CCHOH* + H* ⇌ CH3CH2CCH2OH* + *` | C 端加氢 |
| R0213 | `CH3CH2CCO* + H* ⇌ CH3CH2CCOH* + *` | O 端加氢 |
| R0214 | `CH3CH2CCHO* + H* ⇌ CH3CH2CCHOH* + *` | O 端加氢 |
| R0215 | `CH3CH2CCH2O* + H* ⇌ CH3CH2CCH2OH* + *` | O 端加氢 |
| R0216 | `CH3CH2CHCO* + H* ⇌ CH3CH2CHCHO* + *` | C 端加氢 |
| R0217 | `CH3CH2CHCHO* + H* ⇌ CH3CH2CHCH2O* + *` | C 端加氢 |
| R0218 | `CH3CH2CHCOH* + H* ⇌ CH3CH2CHCHOH* + *` | C 端加氢 |
| R0219 | `CH3CH2CHCHOH* + H* ⇌ CH3CH2CHCH2OH* + *` | C 端加氢 |
| R0220 | `CH3CH2CHCO* + H* ⇌ CH3CH2CHCOH* + *` | O 端加氢 |
| R0221 | `CH3CH2CHCHO* + H* ⇌ CH3CH2CHCHOH* + *` | O 端加氢 |
| R0222 | `CH3CH2CHCH2O* + H* ⇌ CH3CH2CHCH2OH* + *` | O 端加氢 |
| R0223 | `CH3CH2CH2CO* + H* ⇌ CH3CH2CH2CHO* + *` | C 端加氢 |
| R0224 | `CH3CH2CH2CHO* + H* ⇌ CH3CH2CH2CH2O* + *` | C 端加氢 |
| R0225 | `CH3CH2CH2COH* + H* ⇌ CH3CH2CH2CHOH* + *` | C 端加氢 |
| R0226 | `CH3CH2CH2CHOH* + H* ⇌ CH3CH2CH2CH2OH* + *` | C 端加氢 |
| R0227 | `CH3CH2CH2CO* + H* ⇌ CH3CH2CH2COH* + *` | O 端加氢 |
| R0228 | `CH3CH2CH2CHO* + H* ⇌ CH3CH2CH2CHOH* + *` | O 端加氢 |
| R0229 | `CH3CH2CH2CH2O* + H* ⇌ CH3CH2CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0230 | `CH3CH2CCO* + * ⇌ CH3CH2CC* + O*` |  |
| R0231 | `CH3CH2CCOH* + * ⇌ CH3CH2CC* + OH*` |  |
| R0232 | `CH3CH2CCHO* + * ⇌ CH3CH2CCH* + O*` |  |
| R0233 | `CH3CH2CCHOH* + * ⇌ CH3CH2CCH* + OH*` |  |
| R0234 | `CH3CH2CCH2O* + * ⇌ CH3CH2CCH2* + O*` |  |
| R0235 | `CH3CH2CCH2OH* + * ⇌ CH3CH2CCH2* + OH*` |  |
| R0236 | `CH3CH2CHCO* + * ⇌ CH3CH2CHC* + O*` |  |
| R0237 | `CH3CH2CHCOH* + * ⇌ CH3CH2CHC* + OH*` |  |
| R0238 | `CH3CH2CHCHO* + * ⇌ CH3CH2CHCH* + O*` |  |
| R0239 | `CH3CH2CHCHOH* + * ⇌ CH3CH2CHCH* + OH*` |  |
| R0240 | `CH3CH2CHCH2O* + * ⇌ CH3CH2CHCH2* + O*` |  |
| R0241 | `CH3CH2CHCH2OH* + * ⇌ CH3CH2CHCH2* + OH*` |  |
| R0242 | `CH3CH2CH2CO* + * ⇌ CH3(CH2)2C* + O*` |  |
| R0243 | `CH3CH2CH2COH* + * ⇌ CH3(CH2)2C* + OH*` |  |
| R0244 | `CH3CH2CH2CHO* + * ⇌ CH3(CH2)2CH* + O*` |  |
| R0245 | `CH3CH2CH2CHOH* + * ⇌ CH3(CH2)2CH* + OH*` |  |
| R0246 | `CH3CH2CH2CH2O* + * ⇌ CH3(CH2)2CH2* + O*` |  |
| R0247 | `CH3CH2CH2CH2OH* + * ⇌ CH3(CH2)2CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0248 | `CH3CH2COH* + HCOH* ⇌ CH3CH2CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0249 | `CH3CHCH* + CH2* ⇌ CH3CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0250 | `CH3CHCHCH2* ⇌ CH3CH2CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0251 | `CH3CHCHCH2* + H* ⇌ CH3CH2CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0252 | `CH3CH2CC* + H* ⇌ CH3CH2CCH* + *` | α-C 加氢 |
| R0253 | `CH3CH2CC* + H* ⇌ CH3CH2CHC* + *` | β-C 加氢 |
| R0254 | `CH3CH2CCH* + H* ⇌ CH3CH2CCH2* + *` | α-C 加氢 |
| R0255 | `CH3CH2CCH* + H* ⇌ CH3CH2CHCH* + *` | β-C 加氢 |
| R0256 | `CH3CH2CCH2* + H* ⇌ CH3CH2CCH3* + *` | α-C 加氢 |
| R0257 | `CH3CH2CCH2* + H* ⇌ CH3CH2CHCH2* + *` | β-C 加氢 |
| R0258 | `CH3CH2CCH3* + H* ⇌ CH3CH2CHCH3* + *` | β-C 加氢 |
| R0259 | `CH3CH2CHC* + H* ⇌ CH3CH2CHCH* + *` | α-C 加氢 |
| R0260 | `CH3CH2CHC* + H* ⇌ CH3(CH2)2C* + *` | β-C 加氢 |
| R0261 | `CH3CH2CHCH* + H* ⇌ CH3CH2CHCH2* + *` | α-C 加氢 |
| R0262 | `CH3CH2CHCH* + H* ⇌ CH3(CH2)2CH* + *` | β-C 加氢 |
| R0263 | `CH3CH2CHCH2* + H* ⇌ CH3CH2CHCH3* + *` | α-C 加氢 |
| R0264 | `CH3CH2CHCH2* + H* ⇌ CH3(CH2)2CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0265 | `CH3CH2CHCH3* + H* ⇌ CH3(CH2)2CH3* + *` | β-C 加氢 |
| R0266 | `CH3(CH2)2C* + H* ⇌ CH3(CH2)2CH* + *` | α-C 加氢 |
| R0267 | `CH3(CH2)2CH* + H* ⇌ CH3(CH2)2CH2* + *` | α-C 加氢 |
| R0268 | `CH3(CH2)2CH2* + H* ⇌ CH3(CH2)2CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0269 | `CH3(CH2)2CH3* ⇌ C4H10(g) + *` | 正构烷烃脱附 |
| R0270 | `CH3CH2CHCH2* ⇌ C4H8(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C4 → C5

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0271 | `CH3(CH2)2C* + C* ⇌ CH3(CH2)2CC* + *` | 与表面 C 偶联 |
| R0272 | `CH3(CH2)2C* + CH* ⇌ CH3(CH2)2CCH* + *` | 次烷基 + CH |
| R0273 | `CH3(CH2)2C* + CH2* ⇌ CH3(CH2)2CCH2* + *` |  |
| R0274 | `CH3(CH2)2C* + CH3* ⇌ CH3(CH2)2CCH3* + *` |  |
| R0275 | `CH3(CH2)2CH* + C* ⇌ CH3(CH2)2CHC* + *` | 与表面 C 偶联 |
| R0276 | `CH3(CH2)2CH* + CH* ⇌ CH3(CH2)2CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0277 | `CH3(CH2)2CH* + CH2* ⇌ CH3(CH2)2CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0278 | `CH3(CH2)2CH* + CH3* ⇌ CH3(CH2)2CHCH3* + *` |  |
| R0279 | `CH3(CH2)2CH2* + C* ⇌ CH3(CH2)3C* + *` | 与表面 C 偶联 |
| R0280 | `CH3(CH2)2CH2* + CH* ⇌ CH3(CH2)3CH* + *` |  |
| R0281 | `CH3(CH2)2CH2* + CH2* ⇌ CH3(CH2)3CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0282 | `CH3(CH2)2CH2* + CH3* ⇌ CH3(CH2)3CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0283 | `CH3(CH2)2C* + CO* ⇌ CH3(CH2)2CCO* + *` |  |
| R0284 | `CH3(CH2)2CH* + CO* ⇌ CH3(CH2)2CHCO* + *` |  |
| R0285 | `CH3(CH2)2CH2* + CO* ⇌ CH3(CH2)2CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0286 | `CH3(CH2)2CCO* + H* ⇌ CH3(CH2)2CCHO* + *` | C 端加氢 |
| R0287 | `CH3(CH2)2CCHO* + H* ⇌ CH3(CH2)2CCH2O* + *` | C 端加氢 |
| R0288 | `CH3(CH2)2CCOH* + H* ⇌ CH3(CH2)2CCHOH* + *` | C 端加氢 |
| R0289 | `CH3(CH2)2CCHOH* + H* ⇌ CH3(CH2)2CCH2OH* + *` | C 端加氢 |
| R0290 | `CH3(CH2)2CCO* + H* ⇌ CH3(CH2)2CCOH* + *` | O 端加氢 |
| R0291 | `CH3(CH2)2CCHO* + H* ⇌ CH3(CH2)2CCHOH* + *` | O 端加氢 |
| R0292 | `CH3(CH2)2CCH2O* + H* ⇌ CH3(CH2)2CCH2OH* + *` | O 端加氢 |
| R0293 | `CH3(CH2)2CHCO* + H* ⇌ CH3(CH2)2CHCHO* + *` | C 端加氢 |
| R0294 | `CH3(CH2)2CHCHO* + H* ⇌ CH3(CH2)2CHCH2O* + *` | C 端加氢 |
| R0295 | `CH3(CH2)2CHCOH* + H* ⇌ CH3(CH2)2CHCHOH* + *` | C 端加氢 |
| R0296 | `CH3(CH2)2CHCHOH* + H* ⇌ CH3(CH2)2CHCH2OH* + *` | C 端加氢 |
| R0297 | `CH3(CH2)2CHCO* + H* ⇌ CH3(CH2)2CHCOH* + *` | O 端加氢 |
| R0298 | `CH3(CH2)2CHCHO* + H* ⇌ CH3(CH2)2CHCHOH* + *` | O 端加氢 |
| R0299 | `CH3(CH2)2CHCH2O* + H* ⇌ CH3(CH2)2CHCH2OH* + *` | O 端加氢 |
| R0300 | `CH3(CH2)2CH2CO* + H* ⇌ CH3(CH2)2CH2CHO* + *` | C 端加氢 |
| R0301 | `CH3(CH2)2CH2CHO* + H* ⇌ CH3(CH2)2CH2CH2O* + *` | C 端加氢 |
| R0302 | `CH3(CH2)2CH2COH* + H* ⇌ CH3(CH2)2CH2CHOH* + *` | C 端加氢 |
| R0303 | `CH3(CH2)2CH2CHOH* + H* ⇌ CH3(CH2)2CH2CH2OH* + *` | C 端加氢 |
| R0304 | `CH3(CH2)2CH2CO* + H* ⇌ CH3(CH2)2CH2COH* + *` | O 端加氢 |
| R0305 | `CH3(CH2)2CH2CHO* + H* ⇌ CH3(CH2)2CH2CHOH* + *` | O 端加氢 |
| R0306 | `CH3(CH2)2CH2CH2O* + H* ⇌ CH3(CH2)2CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0307 | `CH3(CH2)2CCO* + * ⇌ CH3(CH2)2CC* + O*` |  |
| R0308 | `CH3(CH2)2CCOH* + * ⇌ CH3(CH2)2CC* + OH*` |  |
| R0309 | `CH3(CH2)2CCHO* + * ⇌ CH3(CH2)2CCH* + O*` |  |
| R0310 | `CH3(CH2)2CCHOH* + * ⇌ CH3(CH2)2CCH* + OH*` |  |
| R0311 | `CH3(CH2)2CCH2O* + * ⇌ CH3(CH2)2CCH2* + O*` |  |
| R0312 | `CH3(CH2)2CCH2OH* + * ⇌ CH3(CH2)2CCH2* + OH*` |  |
| R0313 | `CH3(CH2)2CHCO* + * ⇌ CH3(CH2)2CHC* + O*` |  |
| R0314 | `CH3(CH2)2CHCOH* + * ⇌ CH3(CH2)2CHC* + OH*` |  |
| R0315 | `CH3(CH2)2CHCHO* + * ⇌ CH3(CH2)2CHCH* + O*` |  |
| R0316 | `CH3(CH2)2CHCHOH* + * ⇌ CH3(CH2)2CHCH* + OH*` |  |
| R0317 | `CH3(CH2)2CHCH2O* + * ⇌ CH3(CH2)2CHCH2* + O*` |  |
| R0318 | `CH3(CH2)2CHCH2OH* + * ⇌ CH3(CH2)2CHCH2* + OH*` |  |
| R0319 | `CH3(CH2)2CH2CO* + * ⇌ CH3(CH2)3C* + O*` |  |
| R0320 | `CH3(CH2)2CH2COH* + * ⇌ CH3(CH2)3C* + OH*` |  |
| R0321 | `CH3(CH2)2CH2CHO* + * ⇌ CH3(CH2)3CH* + O*` |  |
| R0322 | `CH3(CH2)2CH2CHOH* + * ⇌ CH3(CH2)3CH* + OH*` |  |
| R0323 | `CH3(CH2)2CH2CH2O* + * ⇌ CH3(CH2)3CH2* + O*` |  |
| R0324 | `CH3(CH2)2CH2CH2OH* + * ⇌ CH3(CH2)3CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0325 | `CH3CH2CH2COH* + HCOH* ⇌ CH3(CH2)2CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0326 | `CH3CH2CHCH* + CH2* ⇌ CH3CH2CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0327 | `CH3CH2CHCHCH2* ⇌ CH3(CH2)2CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0328 | `CH3CH2CHCHCH2* + H* ⇌ CH3(CH2)2CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0329 | `CH3(CH2)2CC* + H* ⇌ CH3(CH2)2CCH* + *` | α-C 加氢 |
| R0330 | `CH3(CH2)2CC* + H* ⇌ CH3(CH2)2CHC* + *` | β-C 加氢 |
| R0331 | `CH3(CH2)2CCH* + H* ⇌ CH3(CH2)2CCH2* + *` | α-C 加氢 |
| R0332 | `CH3(CH2)2CCH* + H* ⇌ CH3(CH2)2CHCH* + *` | β-C 加氢 |
| R0333 | `CH3(CH2)2CCH2* + H* ⇌ CH3(CH2)2CCH3* + *` | α-C 加氢 |
| R0334 | `CH3(CH2)2CCH2* + H* ⇌ CH3(CH2)2CHCH2* + *` | β-C 加氢 |
| R0335 | `CH3(CH2)2CCH3* + H* ⇌ CH3(CH2)2CHCH3* + *` | β-C 加氢 |
| R0336 | `CH3(CH2)2CHC* + H* ⇌ CH3(CH2)2CHCH* + *` | α-C 加氢 |
| R0337 | `CH3(CH2)2CHC* + H* ⇌ CH3(CH2)3C* + *` | β-C 加氢 |
| R0338 | `CH3(CH2)2CHCH* + H* ⇌ CH3(CH2)2CHCH2* + *` | α-C 加氢 |
| R0339 | `CH3(CH2)2CHCH* + H* ⇌ CH3(CH2)3CH* + *` | β-C 加氢 |
| R0340 | `CH3(CH2)2CHCH2* + H* ⇌ CH3(CH2)2CHCH3* + *` | α-C 加氢 |
| R0341 | `CH3(CH2)2CHCH2* + H* ⇌ CH3(CH2)3CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0342 | `CH3(CH2)2CHCH3* + H* ⇌ CH3(CH2)3CH3* + *` | β-C 加氢 |
| R0343 | `CH3(CH2)3C* + H* ⇌ CH3(CH2)3CH* + *` | α-C 加氢 |
| R0344 | `CH3(CH2)3CH* + H* ⇌ CH3(CH2)3CH2* + *` | α-C 加氢 |
| R0345 | `CH3(CH2)3CH2* + H* ⇌ CH3(CH2)3CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0346 | `CH3(CH2)3CH3* ⇌ C5H12(g) + *` | 正构烷烃脱附 |
| R0347 | `CH3(CH2)2CHCH2* ⇌ C5H10(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C5 → C6

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0348 | `CH3(CH2)3C* + C* ⇌ CH3(CH2)3CC* + *` | 与表面 C 偶联 |
| R0349 | `CH3(CH2)3C* + CH* ⇌ CH3(CH2)3CCH* + *` | 次烷基 + CH |
| R0350 | `CH3(CH2)3C* + CH2* ⇌ CH3(CH2)3CCH2* + *` |  |
| R0351 | `CH3(CH2)3C* + CH3* ⇌ CH3(CH2)3CCH3* + *` |  |
| R0352 | `CH3(CH2)3CH* + C* ⇌ CH3(CH2)3CHC* + *` | 与表面 C 偶联 |
| R0353 | `CH3(CH2)3CH* + CH* ⇌ CH3(CH2)3CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0354 | `CH3(CH2)3CH* + CH2* ⇌ CH3(CH2)3CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0355 | `CH3(CH2)3CH* + CH3* ⇌ CH3(CH2)3CHCH3* + *` |  |
| R0356 | `CH3(CH2)3CH2* + C* ⇌ CH3(CH2)4C* + *` | 与表面 C 偶联 |
| R0357 | `CH3(CH2)3CH2* + CH* ⇌ CH3(CH2)4CH* + *` |  |
| R0358 | `CH3(CH2)3CH2* + CH2* ⇌ CH3(CH2)4CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0359 | `CH3(CH2)3CH2* + CH3* ⇌ CH3(CH2)4CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0360 | `CH3(CH2)3C* + CO* ⇌ CH3(CH2)3CCO* + *` |  |
| R0361 | `CH3(CH2)3CH* + CO* ⇌ CH3(CH2)3CHCO* + *` |  |
| R0362 | `CH3(CH2)3CH2* + CO* ⇌ CH3(CH2)3CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0363 | `CH3(CH2)3CCO* + H* ⇌ CH3(CH2)3CCHO* + *` | C 端加氢 |
| R0364 | `CH3(CH2)3CCHO* + H* ⇌ CH3(CH2)3CCH2O* + *` | C 端加氢 |
| R0365 | `CH3(CH2)3CCOH* + H* ⇌ CH3(CH2)3CCHOH* + *` | C 端加氢 |
| R0366 | `CH3(CH2)3CCHOH* + H* ⇌ CH3(CH2)3CCH2OH* + *` | C 端加氢 |
| R0367 | `CH3(CH2)3CCO* + H* ⇌ CH3(CH2)3CCOH* + *` | O 端加氢 |
| R0368 | `CH3(CH2)3CCHO* + H* ⇌ CH3(CH2)3CCHOH* + *` | O 端加氢 |
| R0369 | `CH3(CH2)3CCH2O* + H* ⇌ CH3(CH2)3CCH2OH* + *` | O 端加氢 |
| R0370 | `CH3(CH2)3CHCO* + H* ⇌ CH3(CH2)3CHCHO* + *` | C 端加氢 |
| R0371 | `CH3(CH2)3CHCHO* + H* ⇌ CH3(CH2)3CHCH2O* + *` | C 端加氢 |
| R0372 | `CH3(CH2)3CHCOH* + H* ⇌ CH3(CH2)3CHCHOH* + *` | C 端加氢 |
| R0373 | `CH3(CH2)3CHCHOH* + H* ⇌ CH3(CH2)3CHCH2OH* + *` | C 端加氢 |
| R0374 | `CH3(CH2)3CHCO* + H* ⇌ CH3(CH2)3CHCOH* + *` | O 端加氢 |
| R0375 | `CH3(CH2)3CHCHO* + H* ⇌ CH3(CH2)3CHCHOH* + *` | O 端加氢 |
| R0376 | `CH3(CH2)3CHCH2O* + H* ⇌ CH3(CH2)3CHCH2OH* + *` | O 端加氢 |
| R0377 | `CH3(CH2)3CH2CO* + H* ⇌ CH3(CH2)3CH2CHO* + *` | C 端加氢 |
| R0378 | `CH3(CH2)3CH2CHO* + H* ⇌ CH3(CH2)3CH2CH2O* + *` | C 端加氢 |
| R0379 | `CH3(CH2)3CH2COH* + H* ⇌ CH3(CH2)3CH2CHOH* + *` | C 端加氢 |
| R0380 | `CH3(CH2)3CH2CHOH* + H* ⇌ CH3(CH2)3CH2CH2OH* + *` | C 端加氢 |
| R0381 | `CH3(CH2)3CH2CO* + H* ⇌ CH3(CH2)3CH2COH* + *` | O 端加氢 |
| R0382 | `CH3(CH2)3CH2CHO* + H* ⇌ CH3(CH2)3CH2CHOH* + *` | O 端加氢 |
| R0383 | `CH3(CH2)3CH2CH2O* + H* ⇌ CH3(CH2)3CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0384 | `CH3(CH2)3CCO* + * ⇌ CH3(CH2)3CC* + O*` |  |
| R0385 | `CH3(CH2)3CCOH* + * ⇌ CH3(CH2)3CC* + OH*` |  |
| R0386 | `CH3(CH2)3CCHO* + * ⇌ CH3(CH2)3CCH* + O*` |  |
| R0387 | `CH3(CH2)3CCHOH* + * ⇌ CH3(CH2)3CCH* + OH*` |  |
| R0388 | `CH3(CH2)3CCH2O* + * ⇌ CH3(CH2)3CCH2* + O*` |  |
| R0389 | `CH3(CH2)3CCH2OH* + * ⇌ CH3(CH2)3CCH2* + OH*` |  |
| R0390 | `CH3(CH2)3CHCO* + * ⇌ CH3(CH2)3CHC* + O*` |  |
| R0391 | `CH3(CH2)3CHCOH* + * ⇌ CH3(CH2)3CHC* + OH*` |  |
| R0392 | `CH3(CH2)3CHCHO* + * ⇌ CH3(CH2)3CHCH* + O*` |  |
| R0393 | `CH3(CH2)3CHCHOH* + * ⇌ CH3(CH2)3CHCH* + OH*` |  |
| R0394 | `CH3(CH2)3CHCH2O* + * ⇌ CH3(CH2)3CHCH2* + O*` |  |
| R0395 | `CH3(CH2)3CHCH2OH* + * ⇌ CH3(CH2)3CHCH2* + OH*` |  |
| R0396 | `CH3(CH2)3CH2CO* + * ⇌ CH3(CH2)4C* + O*` |  |
| R0397 | `CH3(CH2)3CH2COH* + * ⇌ CH3(CH2)4C* + OH*` |  |
| R0398 | `CH3(CH2)3CH2CHO* + * ⇌ CH3(CH2)4CH* + O*` |  |
| R0399 | `CH3(CH2)3CH2CHOH* + * ⇌ CH3(CH2)4CH* + OH*` |  |
| R0400 | `CH3(CH2)3CH2CH2O* + * ⇌ CH3(CH2)4CH2* + O*` |  |
| R0401 | `CH3(CH2)3CH2CH2OH* + * ⇌ CH3(CH2)4CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0402 | `CH3(CH2)2CH2COH* + HCOH* ⇌ CH3(CH2)3CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0403 | `CH3(CH2)2CHCH* + CH2* ⇌ CH3(CH2)2CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0404 | `CH3(CH2)2CHCHCH2* ⇌ CH3(CH2)3CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0405 | `CH3(CH2)2CHCHCH2* + H* ⇌ CH3(CH2)3CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0406 | `CH3(CH2)3CC* + H* ⇌ CH3(CH2)3CCH* + *` | α-C 加氢 |
| R0407 | `CH3(CH2)3CC* + H* ⇌ CH3(CH2)3CHC* + *` | β-C 加氢 |
| R0408 | `CH3(CH2)3CCH* + H* ⇌ CH3(CH2)3CCH2* + *` | α-C 加氢 |
| R0409 | `CH3(CH2)3CCH* + H* ⇌ CH3(CH2)3CHCH* + *` | β-C 加氢 |
| R0410 | `CH3(CH2)3CCH2* + H* ⇌ CH3(CH2)3CCH3* + *` | α-C 加氢 |
| R0411 | `CH3(CH2)3CCH2* + H* ⇌ CH3(CH2)3CHCH2* + *` | β-C 加氢 |
| R0412 | `CH3(CH2)3CCH3* + H* ⇌ CH3(CH2)3CHCH3* + *` | β-C 加氢 |
| R0413 | `CH3(CH2)3CHC* + H* ⇌ CH3(CH2)3CHCH* + *` | α-C 加氢 |
| R0414 | `CH3(CH2)3CHC* + H* ⇌ CH3(CH2)4C* + *` | β-C 加氢 |
| R0415 | `CH3(CH2)3CHCH* + H* ⇌ CH3(CH2)3CHCH2* + *` | α-C 加氢 |
| R0416 | `CH3(CH2)3CHCH* + H* ⇌ CH3(CH2)4CH* + *` | β-C 加氢 |
| R0417 | `CH3(CH2)3CHCH2* + H* ⇌ CH3(CH2)3CHCH3* + *` | α-C 加氢 |
| R0418 | `CH3(CH2)3CHCH2* + H* ⇌ CH3(CH2)4CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0419 | `CH3(CH2)3CHCH3* + H* ⇌ CH3(CH2)4CH3* + *` | β-C 加氢 |
| R0420 | `CH3(CH2)4C* + H* ⇌ CH3(CH2)4CH* + *` | α-C 加氢 |
| R0421 | `CH3(CH2)4CH* + H* ⇌ CH3(CH2)4CH2* + *` | α-C 加氢 |
| R0422 | `CH3(CH2)4CH2* + H* ⇌ CH3(CH2)4CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0423 | `CH3(CH2)4CH3* ⇌ C6H14(g) + *` | 正构烷烃脱附 |
| R0424 | `CH3(CH2)3CHCH2* ⇌ C6H12(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C6 → C7

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0425 | `CH3(CH2)4C* + C* ⇌ CH3(CH2)4CC* + *` | 与表面 C 偶联 |
| R0426 | `CH3(CH2)4C* + CH* ⇌ CH3(CH2)4CCH* + *` | 次烷基 + CH |
| R0427 | `CH3(CH2)4C* + CH2* ⇌ CH3(CH2)4CCH2* + *` |  |
| R0428 | `CH3(CH2)4C* + CH3* ⇌ CH3(CH2)4CCH3* + *` |  |
| R0429 | `CH3(CH2)4CH* + C* ⇌ CH3(CH2)4CHC* + *` | 与表面 C 偶联 |
| R0430 | `CH3(CH2)4CH* + CH* ⇌ CH3(CH2)4CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0431 | `CH3(CH2)4CH* + CH2* ⇌ CH3(CH2)4CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0432 | `CH3(CH2)4CH* + CH3* ⇌ CH3(CH2)4CHCH3* + *` |  |
| R0433 | `CH3(CH2)4CH2* + C* ⇌ CH3(CH2)5C* + *` | 与表面 C 偶联 |
| R0434 | `CH3(CH2)4CH2* + CH* ⇌ CH3(CH2)5CH* + *` |  |
| R0435 | `CH3(CH2)4CH2* + CH2* ⇌ CH3(CH2)5CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0436 | `CH3(CH2)4CH2* + CH3* ⇌ CH3(CH2)5CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0437 | `CH3(CH2)4C* + CO* ⇌ CH3(CH2)4CCO* + *` |  |
| R0438 | `CH3(CH2)4CH* + CO* ⇌ CH3(CH2)4CHCO* + *` |  |
| R0439 | `CH3(CH2)4CH2* + CO* ⇌ CH3(CH2)4CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0440 | `CH3(CH2)4CCO* + H* ⇌ CH3(CH2)4CCHO* + *` | C 端加氢 |
| R0441 | `CH3(CH2)4CCHO* + H* ⇌ CH3(CH2)4CCH2O* + *` | C 端加氢 |
| R0442 | `CH3(CH2)4CCOH* + H* ⇌ CH3(CH2)4CCHOH* + *` | C 端加氢 |
| R0443 | `CH3(CH2)4CCHOH* + H* ⇌ CH3(CH2)4CCH2OH* + *` | C 端加氢 |
| R0444 | `CH3(CH2)4CCO* + H* ⇌ CH3(CH2)4CCOH* + *` | O 端加氢 |
| R0445 | `CH3(CH2)4CCHO* + H* ⇌ CH3(CH2)4CCHOH* + *` | O 端加氢 |
| R0446 | `CH3(CH2)4CCH2O* + H* ⇌ CH3(CH2)4CCH2OH* + *` | O 端加氢 |
| R0447 | `CH3(CH2)4CHCO* + H* ⇌ CH3(CH2)4CHCHO* + *` | C 端加氢 |
| R0448 | `CH3(CH2)4CHCHO* + H* ⇌ CH3(CH2)4CHCH2O* + *` | C 端加氢 |
| R0449 | `CH3(CH2)4CHCOH* + H* ⇌ CH3(CH2)4CHCHOH* + *` | C 端加氢 |
| R0450 | `CH3(CH2)4CHCHOH* + H* ⇌ CH3(CH2)4CHCH2OH* + *` | C 端加氢 |
| R0451 | `CH3(CH2)4CHCO* + H* ⇌ CH3(CH2)4CHCOH* + *` | O 端加氢 |
| R0452 | `CH3(CH2)4CHCHO* + H* ⇌ CH3(CH2)4CHCHOH* + *` | O 端加氢 |
| R0453 | `CH3(CH2)4CHCH2O* + H* ⇌ CH3(CH2)4CHCH2OH* + *` | O 端加氢 |
| R0454 | `CH3(CH2)4CH2CO* + H* ⇌ CH3(CH2)4CH2CHO* + *` | C 端加氢 |
| R0455 | `CH3(CH2)4CH2CHO* + H* ⇌ CH3(CH2)4CH2CH2O* + *` | C 端加氢 |
| R0456 | `CH3(CH2)4CH2COH* + H* ⇌ CH3(CH2)4CH2CHOH* + *` | C 端加氢 |
| R0457 | `CH3(CH2)4CH2CHOH* + H* ⇌ CH3(CH2)4CH2CH2OH* + *` | C 端加氢 |
| R0458 | `CH3(CH2)4CH2CO* + H* ⇌ CH3(CH2)4CH2COH* + *` | O 端加氢 |
| R0459 | `CH3(CH2)4CH2CHO* + H* ⇌ CH3(CH2)4CH2CHOH* + *` | O 端加氢 |
| R0460 | `CH3(CH2)4CH2CH2O* + H* ⇌ CH3(CH2)4CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0461 | `CH3(CH2)4CCO* + * ⇌ CH3(CH2)4CC* + O*` |  |
| R0462 | `CH3(CH2)4CCOH* + * ⇌ CH3(CH2)4CC* + OH*` |  |
| R0463 | `CH3(CH2)4CCHO* + * ⇌ CH3(CH2)4CCH* + O*` |  |
| R0464 | `CH3(CH2)4CCHOH* + * ⇌ CH3(CH2)4CCH* + OH*` |  |
| R0465 | `CH3(CH2)4CCH2O* + * ⇌ CH3(CH2)4CCH2* + O*` |  |
| R0466 | `CH3(CH2)4CCH2OH* + * ⇌ CH3(CH2)4CCH2* + OH*` |  |
| R0467 | `CH3(CH2)4CHCO* + * ⇌ CH3(CH2)4CHC* + O*` |  |
| R0468 | `CH3(CH2)4CHCOH* + * ⇌ CH3(CH2)4CHC* + OH*` |  |
| R0469 | `CH3(CH2)4CHCHO* + * ⇌ CH3(CH2)4CHCH* + O*` |  |
| R0470 | `CH3(CH2)4CHCHOH* + * ⇌ CH3(CH2)4CHCH* + OH*` |  |
| R0471 | `CH3(CH2)4CHCH2O* + * ⇌ CH3(CH2)4CHCH2* + O*` |  |
| R0472 | `CH3(CH2)4CHCH2OH* + * ⇌ CH3(CH2)4CHCH2* + OH*` |  |
| R0473 | `CH3(CH2)4CH2CO* + * ⇌ CH3(CH2)5C* + O*` |  |
| R0474 | `CH3(CH2)4CH2COH* + * ⇌ CH3(CH2)5C* + OH*` |  |
| R0475 | `CH3(CH2)4CH2CHO* + * ⇌ CH3(CH2)5CH* + O*` |  |
| R0476 | `CH3(CH2)4CH2CHOH* + * ⇌ CH3(CH2)5CH* + OH*` |  |
| R0477 | `CH3(CH2)4CH2CH2O* + * ⇌ CH3(CH2)5CH2* + O*` |  |
| R0478 | `CH3(CH2)4CH2CH2OH* + * ⇌ CH3(CH2)5CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0479 | `CH3(CH2)3CH2COH* + HCOH* ⇌ CH3(CH2)4CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0480 | `CH3(CH2)3CHCH* + CH2* ⇌ CH3(CH2)3CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0481 | `CH3(CH2)3CHCHCH2* ⇌ CH3(CH2)4CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0482 | `CH3(CH2)3CHCHCH2* + H* ⇌ CH3(CH2)4CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0483 | `CH3(CH2)4CC* + H* ⇌ CH3(CH2)4CCH* + *` | α-C 加氢 |
| R0484 | `CH3(CH2)4CC* + H* ⇌ CH3(CH2)4CHC* + *` | β-C 加氢 |
| R0485 | `CH3(CH2)4CCH* + H* ⇌ CH3(CH2)4CCH2* + *` | α-C 加氢 |
| R0486 | `CH3(CH2)4CCH* + H* ⇌ CH3(CH2)4CHCH* + *` | β-C 加氢 |
| R0487 | `CH3(CH2)4CCH2* + H* ⇌ CH3(CH2)4CCH3* + *` | α-C 加氢 |
| R0488 | `CH3(CH2)4CCH2* + H* ⇌ CH3(CH2)4CHCH2* + *` | β-C 加氢 |
| R0489 | `CH3(CH2)4CCH3* + H* ⇌ CH3(CH2)4CHCH3* + *` | β-C 加氢 |
| R0490 | `CH3(CH2)4CHC* + H* ⇌ CH3(CH2)4CHCH* + *` | α-C 加氢 |
| R0491 | `CH3(CH2)4CHC* + H* ⇌ CH3(CH2)5C* + *` | β-C 加氢 |
| R0492 | `CH3(CH2)4CHCH* + H* ⇌ CH3(CH2)4CHCH2* + *` | α-C 加氢 |
| R0493 | `CH3(CH2)4CHCH* + H* ⇌ CH3(CH2)5CH* + *` | β-C 加氢 |
| R0494 | `CH3(CH2)4CHCH2* + H* ⇌ CH3(CH2)4CHCH3* + *` | α-C 加氢 |
| R0495 | `CH3(CH2)4CHCH2* + H* ⇌ CH3(CH2)5CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0496 | `CH3(CH2)4CHCH3* + H* ⇌ CH3(CH2)5CH3* + *` | β-C 加氢 |
| R0497 | `CH3(CH2)5C* + H* ⇌ CH3(CH2)5CH* + *` | α-C 加氢 |
| R0498 | `CH3(CH2)5CH* + H* ⇌ CH3(CH2)5CH2* + *` | α-C 加氢 |
| R0499 | `CH3(CH2)5CH2* + H* ⇌ CH3(CH2)5CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0500 | `CH3(CH2)5CH3* ⇌ C7H16(g) + *` | 正构烷烃脱附 |
| R0501 | `CH3(CH2)4CHCH2* ⇌ C7H14(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C7 → C8

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0502 | `CH3(CH2)5C* + C* ⇌ CH3(CH2)5CC* + *` | 与表面 C 偶联 |
| R0503 | `CH3(CH2)5C* + CH* ⇌ CH3(CH2)5CCH* + *` | 次烷基 + CH |
| R0504 | `CH3(CH2)5C* + CH2* ⇌ CH3(CH2)5CCH2* + *` |  |
| R0505 | `CH3(CH2)5C* + CH3* ⇌ CH3(CH2)5CCH3* + *` |  |
| R0506 | `CH3(CH2)5CH* + C* ⇌ CH3(CH2)5CHC* + *` | 与表面 C 偶联 |
| R0507 | `CH3(CH2)5CH* + CH* ⇌ CH3(CH2)5CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0508 | `CH3(CH2)5CH* + CH2* ⇌ CH3(CH2)5CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0509 | `CH3(CH2)5CH* + CH3* ⇌ CH3(CH2)5CHCH3* + *` |  |
| R0510 | `CH3(CH2)5CH2* + C* ⇌ CH3(CH2)6C* + *` | 与表面 C 偶联 |
| R0511 | `CH3(CH2)5CH2* + CH* ⇌ CH3(CH2)6CH* + *` |  |
| R0512 | `CH3(CH2)5CH2* + CH2* ⇌ CH3(CH2)6CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0513 | `CH3(CH2)5CH2* + CH3* ⇌ CH3(CH2)6CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0514 | `CH3(CH2)5C* + CO* ⇌ CH3(CH2)5CCO* + *` |  |
| R0515 | `CH3(CH2)5CH* + CO* ⇌ CH3(CH2)5CHCO* + *` |  |
| R0516 | `CH3(CH2)5CH2* + CO* ⇌ CH3(CH2)5CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0517 | `CH3(CH2)5CCO* + H* ⇌ CH3(CH2)5CCHO* + *` | C 端加氢 |
| R0518 | `CH3(CH2)5CCHO* + H* ⇌ CH3(CH2)5CCH2O* + *` | C 端加氢 |
| R0519 | `CH3(CH2)5CCOH* + H* ⇌ CH3(CH2)5CCHOH* + *` | C 端加氢 |
| R0520 | `CH3(CH2)5CCHOH* + H* ⇌ CH3(CH2)5CCH2OH* + *` | C 端加氢 |
| R0521 | `CH3(CH2)5CCO* + H* ⇌ CH3(CH2)5CCOH* + *` | O 端加氢 |
| R0522 | `CH3(CH2)5CCHO* + H* ⇌ CH3(CH2)5CCHOH* + *` | O 端加氢 |
| R0523 | `CH3(CH2)5CCH2O* + H* ⇌ CH3(CH2)5CCH2OH* + *` | O 端加氢 |
| R0524 | `CH3(CH2)5CHCO* + H* ⇌ CH3(CH2)5CHCHO* + *` | C 端加氢 |
| R0525 | `CH3(CH2)5CHCHO* + H* ⇌ CH3(CH2)5CHCH2O* + *` | C 端加氢 |
| R0526 | `CH3(CH2)5CHCOH* + H* ⇌ CH3(CH2)5CHCHOH* + *` | C 端加氢 |
| R0527 | `CH3(CH2)5CHCHOH* + H* ⇌ CH3(CH2)5CHCH2OH* + *` | C 端加氢 |
| R0528 | `CH3(CH2)5CHCO* + H* ⇌ CH3(CH2)5CHCOH* + *` | O 端加氢 |
| R0529 | `CH3(CH2)5CHCHO* + H* ⇌ CH3(CH2)5CHCHOH* + *` | O 端加氢 |
| R0530 | `CH3(CH2)5CHCH2O* + H* ⇌ CH3(CH2)5CHCH2OH* + *` | O 端加氢 |
| R0531 | `CH3(CH2)5CH2CO* + H* ⇌ CH3(CH2)5CH2CHO* + *` | C 端加氢 |
| R0532 | `CH3(CH2)5CH2CHO* + H* ⇌ CH3(CH2)5CH2CH2O* + *` | C 端加氢 |
| R0533 | `CH3(CH2)5CH2COH* + H* ⇌ CH3(CH2)5CH2CHOH* + *` | C 端加氢 |
| R0534 | `CH3(CH2)5CH2CHOH* + H* ⇌ CH3(CH2)5CH2CH2OH* + *` | C 端加氢 |
| R0535 | `CH3(CH2)5CH2CO* + H* ⇌ CH3(CH2)5CH2COH* + *` | O 端加氢 |
| R0536 | `CH3(CH2)5CH2CHO* + H* ⇌ CH3(CH2)5CH2CHOH* + *` | O 端加氢 |
| R0537 | `CH3(CH2)5CH2CH2O* + H* ⇌ CH3(CH2)5CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0538 | `CH3(CH2)5CCO* + * ⇌ CH3(CH2)5CC* + O*` |  |
| R0539 | `CH3(CH2)5CCOH* + * ⇌ CH3(CH2)5CC* + OH*` |  |
| R0540 | `CH3(CH2)5CCHO* + * ⇌ CH3(CH2)5CCH* + O*` |  |
| R0541 | `CH3(CH2)5CCHOH* + * ⇌ CH3(CH2)5CCH* + OH*` |  |
| R0542 | `CH3(CH2)5CCH2O* + * ⇌ CH3(CH2)5CCH2* + O*` |  |
| R0543 | `CH3(CH2)5CCH2OH* + * ⇌ CH3(CH2)5CCH2* + OH*` |  |
| R0544 | `CH3(CH2)5CHCO* + * ⇌ CH3(CH2)5CHC* + O*` |  |
| R0545 | `CH3(CH2)5CHCOH* + * ⇌ CH3(CH2)5CHC* + OH*` |  |
| R0546 | `CH3(CH2)5CHCHO* + * ⇌ CH3(CH2)5CHCH* + O*` |  |
| R0547 | `CH3(CH2)5CHCHOH* + * ⇌ CH3(CH2)5CHCH* + OH*` |  |
| R0548 | `CH3(CH2)5CHCH2O* + * ⇌ CH3(CH2)5CHCH2* + O*` |  |
| R0549 | `CH3(CH2)5CHCH2OH* + * ⇌ CH3(CH2)5CHCH2* + OH*` |  |
| R0550 | `CH3(CH2)5CH2CO* + * ⇌ CH3(CH2)6C* + O*` |  |
| R0551 | `CH3(CH2)5CH2COH* + * ⇌ CH3(CH2)6C* + OH*` |  |
| R0552 | `CH3(CH2)5CH2CHO* + * ⇌ CH3(CH2)6CH* + O*` |  |
| R0553 | `CH3(CH2)5CH2CHOH* + * ⇌ CH3(CH2)6CH* + OH*` |  |
| R0554 | `CH3(CH2)5CH2CH2O* + * ⇌ CH3(CH2)6CH2* + O*` |  |
| R0555 | `CH3(CH2)5CH2CH2OH* + * ⇌ CH3(CH2)6CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0556 | `CH3(CH2)4CH2COH* + HCOH* ⇌ CH3(CH2)5CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0557 | `CH3(CH2)4CHCH* + CH2* ⇌ CH3(CH2)4CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0558 | `CH3(CH2)4CHCHCH2* ⇌ CH3(CH2)5CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0559 | `CH3(CH2)4CHCHCH2* + H* ⇌ CH3(CH2)5CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0560 | `CH3(CH2)5CC* + H* ⇌ CH3(CH2)5CCH* + *` | α-C 加氢 |
| R0561 | `CH3(CH2)5CC* + H* ⇌ CH3(CH2)5CHC* + *` | β-C 加氢 |
| R0562 | `CH3(CH2)5CCH* + H* ⇌ CH3(CH2)5CCH2* + *` | α-C 加氢 |
| R0563 | `CH3(CH2)5CCH* + H* ⇌ CH3(CH2)5CHCH* + *` | β-C 加氢 |
| R0564 | `CH3(CH2)5CCH2* + H* ⇌ CH3(CH2)5CCH3* + *` | α-C 加氢 |
| R0565 | `CH3(CH2)5CCH2* + H* ⇌ CH3(CH2)5CHCH2* + *` | β-C 加氢 |
| R0566 | `CH3(CH2)5CCH3* + H* ⇌ CH3(CH2)5CHCH3* + *` | β-C 加氢 |
| R0567 | `CH3(CH2)5CHC* + H* ⇌ CH3(CH2)5CHCH* + *` | α-C 加氢 |
| R0568 | `CH3(CH2)5CHC* + H* ⇌ CH3(CH2)6C* + *` | β-C 加氢 |
| R0569 | `CH3(CH2)5CHCH* + H* ⇌ CH3(CH2)5CHCH2* + *` | α-C 加氢 |
| R0570 | `CH3(CH2)5CHCH* + H* ⇌ CH3(CH2)6CH* + *` | β-C 加氢 |
| R0571 | `CH3(CH2)5CHCH2* + H* ⇌ CH3(CH2)5CHCH3* + *` | α-C 加氢 |
| R0572 | `CH3(CH2)5CHCH2* + H* ⇌ CH3(CH2)6CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0573 | `CH3(CH2)5CHCH3* + H* ⇌ CH3(CH2)6CH3* + *` | β-C 加氢 |
| R0574 | `CH3(CH2)6C* + H* ⇌ CH3(CH2)6CH* + *` | α-C 加氢 |
| R0575 | `CH3(CH2)6CH* + H* ⇌ CH3(CH2)6CH2* + *` | α-C 加氢 |
| R0576 | `CH3(CH2)6CH2* + H* ⇌ CH3(CH2)6CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0577 | `CH3(CH2)6CH3* ⇌ C8H18(g) + *` | 正构烷烃脱附 |
| R0578 | `CH3(CH2)5CHCH2* ⇌ C8H16(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |
