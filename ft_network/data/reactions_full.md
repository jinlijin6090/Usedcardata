# 费托合成直链烃基元反应网络全表 (C1–C20)

由 `scripts/build_network.py --nmax 20` 自动生成；共 654 个物种、1502 个可逆基元步骤，每一步均通过 C/H/O 原子守恒与表面位点守恒校验。`*` 为表面空位。


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

## C8 → C9

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0579 | `CH3(CH2)6C* + C* ⇌ CH3(CH2)6CC* + *` | 与表面 C 偶联 |
| R0580 | `CH3(CH2)6C* + CH* ⇌ CH3(CH2)6CCH* + *` | 次烷基 + CH |
| R0581 | `CH3(CH2)6C* + CH2* ⇌ CH3(CH2)6CCH2* + *` |  |
| R0582 | `CH3(CH2)6C* + CH3* ⇌ CH3(CH2)6CCH3* + *` |  |
| R0583 | `CH3(CH2)6CH* + C* ⇌ CH3(CH2)6CHC* + *` | 与表面 C 偶联 |
| R0584 | `CH3(CH2)6CH* + CH* ⇌ CH3(CH2)6CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0585 | `CH3(CH2)6CH* + CH2* ⇌ CH3(CH2)6CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0586 | `CH3(CH2)6CH* + CH3* ⇌ CH3(CH2)6CHCH3* + *` |  |
| R0587 | `CH3(CH2)6CH2* + C* ⇌ CH3(CH2)7C* + *` | 与表面 C 偶联 |
| R0588 | `CH3(CH2)6CH2* + CH* ⇌ CH3(CH2)7CH* + *` |  |
| R0589 | `CH3(CH2)6CH2* + CH2* ⇌ CH3(CH2)7CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0590 | `CH3(CH2)6CH2* + CH3* ⇌ CH3(CH2)7CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0591 | `CH3(CH2)6C* + CO* ⇌ CH3(CH2)6CCO* + *` |  |
| R0592 | `CH3(CH2)6CH* + CO* ⇌ CH3(CH2)6CHCO* + *` |  |
| R0593 | `CH3(CH2)6CH2* + CO* ⇌ CH3(CH2)6CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0594 | `CH3(CH2)6CCO* + H* ⇌ CH3(CH2)6CCHO* + *` | C 端加氢 |
| R0595 | `CH3(CH2)6CCHO* + H* ⇌ CH3(CH2)6CCH2O* + *` | C 端加氢 |
| R0596 | `CH3(CH2)6CCOH* + H* ⇌ CH3(CH2)6CCHOH* + *` | C 端加氢 |
| R0597 | `CH3(CH2)6CCHOH* + H* ⇌ CH3(CH2)6CCH2OH* + *` | C 端加氢 |
| R0598 | `CH3(CH2)6CCO* + H* ⇌ CH3(CH2)6CCOH* + *` | O 端加氢 |
| R0599 | `CH3(CH2)6CCHO* + H* ⇌ CH3(CH2)6CCHOH* + *` | O 端加氢 |
| R0600 | `CH3(CH2)6CCH2O* + H* ⇌ CH3(CH2)6CCH2OH* + *` | O 端加氢 |
| R0601 | `CH3(CH2)6CHCO* + H* ⇌ CH3(CH2)6CHCHO* + *` | C 端加氢 |
| R0602 | `CH3(CH2)6CHCHO* + H* ⇌ CH3(CH2)6CHCH2O* + *` | C 端加氢 |
| R0603 | `CH3(CH2)6CHCOH* + H* ⇌ CH3(CH2)6CHCHOH* + *` | C 端加氢 |
| R0604 | `CH3(CH2)6CHCHOH* + H* ⇌ CH3(CH2)6CHCH2OH* + *` | C 端加氢 |
| R0605 | `CH3(CH2)6CHCO* + H* ⇌ CH3(CH2)6CHCOH* + *` | O 端加氢 |
| R0606 | `CH3(CH2)6CHCHO* + H* ⇌ CH3(CH2)6CHCHOH* + *` | O 端加氢 |
| R0607 | `CH3(CH2)6CHCH2O* + H* ⇌ CH3(CH2)6CHCH2OH* + *` | O 端加氢 |
| R0608 | `CH3(CH2)6CH2CO* + H* ⇌ CH3(CH2)6CH2CHO* + *` | C 端加氢 |
| R0609 | `CH3(CH2)6CH2CHO* + H* ⇌ CH3(CH2)6CH2CH2O* + *` | C 端加氢 |
| R0610 | `CH3(CH2)6CH2COH* + H* ⇌ CH3(CH2)6CH2CHOH* + *` | C 端加氢 |
| R0611 | `CH3(CH2)6CH2CHOH* + H* ⇌ CH3(CH2)6CH2CH2OH* + *` | C 端加氢 |
| R0612 | `CH3(CH2)6CH2CO* + H* ⇌ CH3(CH2)6CH2COH* + *` | O 端加氢 |
| R0613 | `CH3(CH2)6CH2CHO* + H* ⇌ CH3(CH2)6CH2CHOH* + *` | O 端加氢 |
| R0614 | `CH3(CH2)6CH2CH2O* + H* ⇌ CH3(CH2)6CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0615 | `CH3(CH2)6CCO* + * ⇌ CH3(CH2)6CC* + O*` |  |
| R0616 | `CH3(CH2)6CCOH* + * ⇌ CH3(CH2)6CC* + OH*` |  |
| R0617 | `CH3(CH2)6CCHO* + * ⇌ CH3(CH2)6CCH* + O*` |  |
| R0618 | `CH3(CH2)6CCHOH* + * ⇌ CH3(CH2)6CCH* + OH*` |  |
| R0619 | `CH3(CH2)6CCH2O* + * ⇌ CH3(CH2)6CCH2* + O*` |  |
| R0620 | `CH3(CH2)6CCH2OH* + * ⇌ CH3(CH2)6CCH2* + OH*` |  |
| R0621 | `CH3(CH2)6CHCO* + * ⇌ CH3(CH2)6CHC* + O*` |  |
| R0622 | `CH3(CH2)6CHCOH* + * ⇌ CH3(CH2)6CHC* + OH*` |  |
| R0623 | `CH3(CH2)6CHCHO* + * ⇌ CH3(CH2)6CHCH* + O*` |  |
| R0624 | `CH3(CH2)6CHCHOH* + * ⇌ CH3(CH2)6CHCH* + OH*` |  |
| R0625 | `CH3(CH2)6CHCH2O* + * ⇌ CH3(CH2)6CHCH2* + O*` |  |
| R0626 | `CH3(CH2)6CHCH2OH* + * ⇌ CH3(CH2)6CHCH2* + OH*` |  |
| R0627 | `CH3(CH2)6CH2CO* + * ⇌ CH3(CH2)7C* + O*` |  |
| R0628 | `CH3(CH2)6CH2COH* + * ⇌ CH3(CH2)7C* + OH*` |  |
| R0629 | `CH3(CH2)6CH2CHO* + * ⇌ CH3(CH2)7CH* + O*` |  |
| R0630 | `CH3(CH2)6CH2CHOH* + * ⇌ CH3(CH2)7CH* + OH*` |  |
| R0631 | `CH3(CH2)6CH2CH2O* + * ⇌ CH3(CH2)7CH2* + O*` |  |
| R0632 | `CH3(CH2)6CH2CH2OH* + * ⇌ CH3(CH2)7CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0633 | `CH3(CH2)5CH2COH* + HCOH* ⇌ CH3(CH2)6CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0634 | `CH3(CH2)5CHCH* + CH2* ⇌ CH3(CH2)5CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0635 | `CH3(CH2)5CHCHCH2* ⇌ CH3(CH2)6CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0636 | `CH3(CH2)5CHCHCH2* + H* ⇌ CH3(CH2)6CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0637 | `CH3(CH2)6CC* + H* ⇌ CH3(CH2)6CCH* + *` | α-C 加氢 |
| R0638 | `CH3(CH2)6CC* + H* ⇌ CH3(CH2)6CHC* + *` | β-C 加氢 |
| R0639 | `CH3(CH2)6CCH* + H* ⇌ CH3(CH2)6CCH2* + *` | α-C 加氢 |
| R0640 | `CH3(CH2)6CCH* + H* ⇌ CH3(CH2)6CHCH* + *` | β-C 加氢 |
| R0641 | `CH3(CH2)6CCH2* + H* ⇌ CH3(CH2)6CCH3* + *` | α-C 加氢 |
| R0642 | `CH3(CH2)6CCH2* + H* ⇌ CH3(CH2)6CHCH2* + *` | β-C 加氢 |
| R0643 | `CH3(CH2)6CCH3* + H* ⇌ CH3(CH2)6CHCH3* + *` | β-C 加氢 |
| R0644 | `CH3(CH2)6CHC* + H* ⇌ CH3(CH2)6CHCH* + *` | α-C 加氢 |
| R0645 | `CH3(CH2)6CHC* + H* ⇌ CH3(CH2)7C* + *` | β-C 加氢 |
| R0646 | `CH3(CH2)6CHCH* + H* ⇌ CH3(CH2)6CHCH2* + *` | α-C 加氢 |
| R0647 | `CH3(CH2)6CHCH* + H* ⇌ CH3(CH2)7CH* + *` | β-C 加氢 |
| R0648 | `CH3(CH2)6CHCH2* + H* ⇌ CH3(CH2)6CHCH3* + *` | α-C 加氢 |
| R0649 | `CH3(CH2)6CHCH2* + H* ⇌ CH3(CH2)7CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0650 | `CH3(CH2)6CHCH3* + H* ⇌ CH3(CH2)7CH3* + *` | β-C 加氢 |
| R0651 | `CH3(CH2)7C* + H* ⇌ CH3(CH2)7CH* + *` | α-C 加氢 |
| R0652 | `CH3(CH2)7CH* + H* ⇌ CH3(CH2)7CH2* + *` | α-C 加氢 |
| R0653 | `CH3(CH2)7CH2* + H* ⇌ CH3(CH2)7CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0654 | `CH3(CH2)7CH3* ⇌ C9H20(g) + *` | 正构烷烃脱附 |
| R0655 | `CH3(CH2)6CHCH2* ⇌ C9H18(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C9 → C10

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0656 | `CH3(CH2)7C* + C* ⇌ CH3(CH2)7CC* + *` | 与表面 C 偶联 |
| R0657 | `CH3(CH2)7C* + CH* ⇌ CH3(CH2)7CCH* + *` | 次烷基 + CH |
| R0658 | `CH3(CH2)7C* + CH2* ⇌ CH3(CH2)7CCH2* + *` |  |
| R0659 | `CH3(CH2)7C* + CH3* ⇌ CH3(CH2)7CCH3* + *` |  |
| R0660 | `CH3(CH2)7CH* + C* ⇌ CH3(CH2)7CHC* + *` | 与表面 C 偶联 |
| R0661 | `CH3(CH2)7CH* + CH* ⇌ CH3(CH2)7CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0662 | `CH3(CH2)7CH* + CH2* ⇌ CH3(CH2)7CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0663 | `CH3(CH2)7CH* + CH3* ⇌ CH3(CH2)7CHCH3* + *` |  |
| R0664 | `CH3(CH2)7CH2* + C* ⇌ CH3(CH2)8C* + *` | 与表面 C 偶联 |
| R0665 | `CH3(CH2)7CH2* + CH* ⇌ CH3(CH2)8CH* + *` |  |
| R0666 | `CH3(CH2)7CH2* + CH2* ⇌ CH3(CH2)8CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0667 | `CH3(CH2)7CH2* + CH3* ⇌ CH3(CH2)8CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0668 | `CH3(CH2)7C* + CO* ⇌ CH3(CH2)7CCO* + *` |  |
| R0669 | `CH3(CH2)7CH* + CO* ⇌ CH3(CH2)7CHCO* + *` |  |
| R0670 | `CH3(CH2)7CH2* + CO* ⇌ CH3(CH2)7CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0671 | `CH3(CH2)7CCO* + H* ⇌ CH3(CH2)7CCHO* + *` | C 端加氢 |
| R0672 | `CH3(CH2)7CCHO* + H* ⇌ CH3(CH2)7CCH2O* + *` | C 端加氢 |
| R0673 | `CH3(CH2)7CCOH* + H* ⇌ CH3(CH2)7CCHOH* + *` | C 端加氢 |
| R0674 | `CH3(CH2)7CCHOH* + H* ⇌ CH3(CH2)7CCH2OH* + *` | C 端加氢 |
| R0675 | `CH3(CH2)7CCO* + H* ⇌ CH3(CH2)7CCOH* + *` | O 端加氢 |
| R0676 | `CH3(CH2)7CCHO* + H* ⇌ CH3(CH2)7CCHOH* + *` | O 端加氢 |
| R0677 | `CH3(CH2)7CCH2O* + H* ⇌ CH3(CH2)7CCH2OH* + *` | O 端加氢 |
| R0678 | `CH3(CH2)7CHCO* + H* ⇌ CH3(CH2)7CHCHO* + *` | C 端加氢 |
| R0679 | `CH3(CH2)7CHCHO* + H* ⇌ CH3(CH2)7CHCH2O* + *` | C 端加氢 |
| R0680 | `CH3(CH2)7CHCOH* + H* ⇌ CH3(CH2)7CHCHOH* + *` | C 端加氢 |
| R0681 | `CH3(CH2)7CHCHOH* + H* ⇌ CH3(CH2)7CHCH2OH* + *` | C 端加氢 |
| R0682 | `CH3(CH2)7CHCO* + H* ⇌ CH3(CH2)7CHCOH* + *` | O 端加氢 |
| R0683 | `CH3(CH2)7CHCHO* + H* ⇌ CH3(CH2)7CHCHOH* + *` | O 端加氢 |
| R0684 | `CH3(CH2)7CHCH2O* + H* ⇌ CH3(CH2)7CHCH2OH* + *` | O 端加氢 |
| R0685 | `CH3(CH2)7CH2CO* + H* ⇌ CH3(CH2)7CH2CHO* + *` | C 端加氢 |
| R0686 | `CH3(CH2)7CH2CHO* + H* ⇌ CH3(CH2)7CH2CH2O* + *` | C 端加氢 |
| R0687 | `CH3(CH2)7CH2COH* + H* ⇌ CH3(CH2)7CH2CHOH* + *` | C 端加氢 |
| R0688 | `CH3(CH2)7CH2CHOH* + H* ⇌ CH3(CH2)7CH2CH2OH* + *` | C 端加氢 |
| R0689 | `CH3(CH2)7CH2CO* + H* ⇌ CH3(CH2)7CH2COH* + *` | O 端加氢 |
| R0690 | `CH3(CH2)7CH2CHO* + H* ⇌ CH3(CH2)7CH2CHOH* + *` | O 端加氢 |
| R0691 | `CH3(CH2)7CH2CH2O* + H* ⇌ CH3(CH2)7CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0692 | `CH3(CH2)7CCO* + * ⇌ CH3(CH2)7CC* + O*` |  |
| R0693 | `CH3(CH2)7CCOH* + * ⇌ CH3(CH2)7CC* + OH*` |  |
| R0694 | `CH3(CH2)7CCHO* + * ⇌ CH3(CH2)7CCH* + O*` |  |
| R0695 | `CH3(CH2)7CCHOH* + * ⇌ CH3(CH2)7CCH* + OH*` |  |
| R0696 | `CH3(CH2)7CCH2O* + * ⇌ CH3(CH2)7CCH2* + O*` |  |
| R0697 | `CH3(CH2)7CCH2OH* + * ⇌ CH3(CH2)7CCH2* + OH*` |  |
| R0698 | `CH3(CH2)7CHCO* + * ⇌ CH3(CH2)7CHC* + O*` |  |
| R0699 | `CH3(CH2)7CHCOH* + * ⇌ CH3(CH2)7CHC* + OH*` |  |
| R0700 | `CH3(CH2)7CHCHO* + * ⇌ CH3(CH2)7CHCH* + O*` |  |
| R0701 | `CH3(CH2)7CHCHOH* + * ⇌ CH3(CH2)7CHCH* + OH*` |  |
| R0702 | `CH3(CH2)7CHCH2O* + * ⇌ CH3(CH2)7CHCH2* + O*` |  |
| R0703 | `CH3(CH2)7CHCH2OH* + * ⇌ CH3(CH2)7CHCH2* + OH*` |  |
| R0704 | `CH3(CH2)7CH2CO* + * ⇌ CH3(CH2)8C* + O*` |  |
| R0705 | `CH3(CH2)7CH2COH* + * ⇌ CH3(CH2)8C* + OH*` |  |
| R0706 | `CH3(CH2)7CH2CHO* + * ⇌ CH3(CH2)8CH* + O*` |  |
| R0707 | `CH3(CH2)7CH2CHOH* + * ⇌ CH3(CH2)8CH* + OH*` |  |
| R0708 | `CH3(CH2)7CH2CH2O* + * ⇌ CH3(CH2)8CH2* + O*` |  |
| R0709 | `CH3(CH2)7CH2CH2OH* + * ⇌ CH3(CH2)8CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0710 | `CH3(CH2)6CH2COH* + HCOH* ⇌ CH3(CH2)7CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0711 | `CH3(CH2)6CHCH* + CH2* ⇌ CH3(CH2)6CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0712 | `CH3(CH2)6CHCHCH2* ⇌ CH3(CH2)7CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0713 | `CH3(CH2)6CHCHCH2* + H* ⇌ CH3(CH2)7CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0714 | `CH3(CH2)7CC* + H* ⇌ CH3(CH2)7CCH* + *` | α-C 加氢 |
| R0715 | `CH3(CH2)7CC* + H* ⇌ CH3(CH2)7CHC* + *` | β-C 加氢 |
| R0716 | `CH3(CH2)7CCH* + H* ⇌ CH3(CH2)7CCH2* + *` | α-C 加氢 |
| R0717 | `CH3(CH2)7CCH* + H* ⇌ CH3(CH2)7CHCH* + *` | β-C 加氢 |
| R0718 | `CH3(CH2)7CCH2* + H* ⇌ CH3(CH2)7CCH3* + *` | α-C 加氢 |
| R0719 | `CH3(CH2)7CCH2* + H* ⇌ CH3(CH2)7CHCH2* + *` | β-C 加氢 |
| R0720 | `CH3(CH2)7CCH3* + H* ⇌ CH3(CH2)7CHCH3* + *` | β-C 加氢 |
| R0721 | `CH3(CH2)7CHC* + H* ⇌ CH3(CH2)7CHCH* + *` | α-C 加氢 |
| R0722 | `CH3(CH2)7CHC* + H* ⇌ CH3(CH2)8C* + *` | β-C 加氢 |
| R0723 | `CH3(CH2)7CHCH* + H* ⇌ CH3(CH2)7CHCH2* + *` | α-C 加氢 |
| R0724 | `CH3(CH2)7CHCH* + H* ⇌ CH3(CH2)8CH* + *` | β-C 加氢 |
| R0725 | `CH3(CH2)7CHCH2* + H* ⇌ CH3(CH2)7CHCH3* + *` | α-C 加氢 |
| R0726 | `CH3(CH2)7CHCH2* + H* ⇌ CH3(CH2)8CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0727 | `CH3(CH2)7CHCH3* + H* ⇌ CH3(CH2)8CH3* + *` | β-C 加氢 |
| R0728 | `CH3(CH2)8C* + H* ⇌ CH3(CH2)8CH* + *` | α-C 加氢 |
| R0729 | `CH3(CH2)8CH* + H* ⇌ CH3(CH2)8CH2* + *` | α-C 加氢 |
| R0730 | `CH3(CH2)8CH2* + H* ⇌ CH3(CH2)8CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0731 | `CH3(CH2)8CH3* ⇌ C10H22(g) + *` | 正构烷烃脱附 |
| R0732 | `CH3(CH2)7CHCH2* ⇌ C10H20(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C10 → C11

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0733 | `CH3(CH2)8C* + C* ⇌ CH3(CH2)8CC* + *` | 与表面 C 偶联 |
| R0734 | `CH3(CH2)8C* + CH* ⇌ CH3(CH2)8CCH* + *` | 次烷基 + CH |
| R0735 | `CH3(CH2)8C* + CH2* ⇌ CH3(CH2)8CCH2* + *` |  |
| R0736 | `CH3(CH2)8C* + CH3* ⇌ CH3(CH2)8CCH3* + *` |  |
| R0737 | `CH3(CH2)8CH* + C* ⇌ CH3(CH2)8CHC* + *` | 与表面 C 偶联 |
| R0738 | `CH3(CH2)8CH* + CH* ⇌ CH3(CH2)8CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0739 | `CH3(CH2)8CH* + CH2* ⇌ CH3(CH2)8CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0740 | `CH3(CH2)8CH* + CH3* ⇌ CH3(CH2)8CHCH3* + *` |  |
| R0741 | `CH3(CH2)8CH2* + C* ⇌ CH3(CH2)9C* + *` | 与表面 C 偶联 |
| R0742 | `CH3(CH2)8CH2* + CH* ⇌ CH3(CH2)9CH* + *` |  |
| R0743 | `CH3(CH2)8CH2* + CH2* ⇌ CH3(CH2)9CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0744 | `CH3(CH2)8CH2* + CH3* ⇌ CH3(CH2)9CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0745 | `CH3(CH2)8C* + CO* ⇌ CH3(CH2)8CCO* + *` |  |
| R0746 | `CH3(CH2)8CH* + CO* ⇌ CH3(CH2)8CHCO* + *` |  |
| R0747 | `CH3(CH2)8CH2* + CO* ⇌ CH3(CH2)8CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0748 | `CH3(CH2)8CCO* + H* ⇌ CH3(CH2)8CCHO* + *` | C 端加氢 |
| R0749 | `CH3(CH2)8CCHO* + H* ⇌ CH3(CH2)8CCH2O* + *` | C 端加氢 |
| R0750 | `CH3(CH2)8CCOH* + H* ⇌ CH3(CH2)8CCHOH* + *` | C 端加氢 |
| R0751 | `CH3(CH2)8CCHOH* + H* ⇌ CH3(CH2)8CCH2OH* + *` | C 端加氢 |
| R0752 | `CH3(CH2)8CCO* + H* ⇌ CH3(CH2)8CCOH* + *` | O 端加氢 |
| R0753 | `CH3(CH2)8CCHO* + H* ⇌ CH3(CH2)8CCHOH* + *` | O 端加氢 |
| R0754 | `CH3(CH2)8CCH2O* + H* ⇌ CH3(CH2)8CCH2OH* + *` | O 端加氢 |
| R0755 | `CH3(CH2)8CHCO* + H* ⇌ CH3(CH2)8CHCHO* + *` | C 端加氢 |
| R0756 | `CH3(CH2)8CHCHO* + H* ⇌ CH3(CH2)8CHCH2O* + *` | C 端加氢 |
| R0757 | `CH3(CH2)8CHCOH* + H* ⇌ CH3(CH2)8CHCHOH* + *` | C 端加氢 |
| R0758 | `CH3(CH2)8CHCHOH* + H* ⇌ CH3(CH2)8CHCH2OH* + *` | C 端加氢 |
| R0759 | `CH3(CH2)8CHCO* + H* ⇌ CH3(CH2)8CHCOH* + *` | O 端加氢 |
| R0760 | `CH3(CH2)8CHCHO* + H* ⇌ CH3(CH2)8CHCHOH* + *` | O 端加氢 |
| R0761 | `CH3(CH2)8CHCH2O* + H* ⇌ CH3(CH2)8CHCH2OH* + *` | O 端加氢 |
| R0762 | `CH3(CH2)8CH2CO* + H* ⇌ CH3(CH2)8CH2CHO* + *` | C 端加氢 |
| R0763 | `CH3(CH2)8CH2CHO* + H* ⇌ CH3(CH2)8CH2CH2O* + *` | C 端加氢 |
| R0764 | `CH3(CH2)8CH2COH* + H* ⇌ CH3(CH2)8CH2CHOH* + *` | C 端加氢 |
| R0765 | `CH3(CH2)8CH2CHOH* + H* ⇌ CH3(CH2)8CH2CH2OH* + *` | C 端加氢 |
| R0766 | `CH3(CH2)8CH2CO* + H* ⇌ CH3(CH2)8CH2COH* + *` | O 端加氢 |
| R0767 | `CH3(CH2)8CH2CHO* + H* ⇌ CH3(CH2)8CH2CHOH* + *` | O 端加氢 |
| R0768 | `CH3(CH2)8CH2CH2O* + H* ⇌ CH3(CH2)8CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0769 | `CH3(CH2)8CCO* + * ⇌ CH3(CH2)8CC* + O*` |  |
| R0770 | `CH3(CH2)8CCOH* + * ⇌ CH3(CH2)8CC* + OH*` |  |
| R0771 | `CH3(CH2)8CCHO* + * ⇌ CH3(CH2)8CCH* + O*` |  |
| R0772 | `CH3(CH2)8CCHOH* + * ⇌ CH3(CH2)8CCH* + OH*` |  |
| R0773 | `CH3(CH2)8CCH2O* + * ⇌ CH3(CH2)8CCH2* + O*` |  |
| R0774 | `CH3(CH2)8CCH2OH* + * ⇌ CH3(CH2)8CCH2* + OH*` |  |
| R0775 | `CH3(CH2)8CHCO* + * ⇌ CH3(CH2)8CHC* + O*` |  |
| R0776 | `CH3(CH2)8CHCOH* + * ⇌ CH3(CH2)8CHC* + OH*` |  |
| R0777 | `CH3(CH2)8CHCHO* + * ⇌ CH3(CH2)8CHCH* + O*` |  |
| R0778 | `CH3(CH2)8CHCHOH* + * ⇌ CH3(CH2)8CHCH* + OH*` |  |
| R0779 | `CH3(CH2)8CHCH2O* + * ⇌ CH3(CH2)8CHCH2* + O*` |  |
| R0780 | `CH3(CH2)8CHCH2OH* + * ⇌ CH3(CH2)8CHCH2* + OH*` |  |
| R0781 | `CH3(CH2)8CH2CO* + * ⇌ CH3(CH2)9C* + O*` |  |
| R0782 | `CH3(CH2)8CH2COH* + * ⇌ CH3(CH2)9C* + OH*` |  |
| R0783 | `CH3(CH2)8CH2CHO* + * ⇌ CH3(CH2)9CH* + O*` |  |
| R0784 | `CH3(CH2)8CH2CHOH* + * ⇌ CH3(CH2)9CH* + OH*` |  |
| R0785 | `CH3(CH2)8CH2CH2O* + * ⇌ CH3(CH2)9CH2* + O*` |  |
| R0786 | `CH3(CH2)8CH2CH2OH* + * ⇌ CH3(CH2)9CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0787 | `CH3(CH2)7CH2COH* + HCOH* ⇌ CH3(CH2)8CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0788 | `CH3(CH2)7CHCH* + CH2* ⇌ CH3(CH2)7CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0789 | `CH3(CH2)7CHCHCH2* ⇌ CH3(CH2)8CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0790 | `CH3(CH2)7CHCHCH2* + H* ⇌ CH3(CH2)8CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0791 | `CH3(CH2)8CC* + H* ⇌ CH3(CH2)8CCH* + *` | α-C 加氢 |
| R0792 | `CH3(CH2)8CC* + H* ⇌ CH3(CH2)8CHC* + *` | β-C 加氢 |
| R0793 | `CH3(CH2)8CCH* + H* ⇌ CH3(CH2)8CCH2* + *` | α-C 加氢 |
| R0794 | `CH3(CH2)8CCH* + H* ⇌ CH3(CH2)8CHCH* + *` | β-C 加氢 |
| R0795 | `CH3(CH2)8CCH2* + H* ⇌ CH3(CH2)8CCH3* + *` | α-C 加氢 |
| R0796 | `CH3(CH2)8CCH2* + H* ⇌ CH3(CH2)8CHCH2* + *` | β-C 加氢 |
| R0797 | `CH3(CH2)8CCH3* + H* ⇌ CH3(CH2)8CHCH3* + *` | β-C 加氢 |
| R0798 | `CH3(CH2)8CHC* + H* ⇌ CH3(CH2)8CHCH* + *` | α-C 加氢 |
| R0799 | `CH3(CH2)8CHC* + H* ⇌ CH3(CH2)9C* + *` | β-C 加氢 |
| R0800 | `CH3(CH2)8CHCH* + H* ⇌ CH3(CH2)8CHCH2* + *` | α-C 加氢 |
| R0801 | `CH3(CH2)8CHCH* + H* ⇌ CH3(CH2)9CH* + *` | β-C 加氢 |
| R0802 | `CH3(CH2)8CHCH2* + H* ⇌ CH3(CH2)8CHCH3* + *` | α-C 加氢 |
| R0803 | `CH3(CH2)8CHCH2* + H* ⇌ CH3(CH2)9CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0804 | `CH3(CH2)8CHCH3* + H* ⇌ CH3(CH2)9CH3* + *` | β-C 加氢 |
| R0805 | `CH3(CH2)9C* + H* ⇌ CH3(CH2)9CH* + *` | α-C 加氢 |
| R0806 | `CH3(CH2)9CH* + H* ⇌ CH3(CH2)9CH2* + *` | α-C 加氢 |
| R0807 | `CH3(CH2)9CH2* + H* ⇌ CH3(CH2)9CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0808 | `CH3(CH2)9CH3* ⇌ C11H24(g) + *` | 正构烷烃脱附 |
| R0809 | `CH3(CH2)8CHCH2* ⇌ C11H22(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C11 → C12

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0810 | `CH3(CH2)9C* + C* ⇌ CH3(CH2)9CC* + *` | 与表面 C 偶联 |
| R0811 | `CH3(CH2)9C* + CH* ⇌ CH3(CH2)9CCH* + *` | 次烷基 + CH |
| R0812 | `CH3(CH2)9C* + CH2* ⇌ CH3(CH2)9CCH2* + *` |  |
| R0813 | `CH3(CH2)9C* + CH3* ⇌ CH3(CH2)9CCH3* + *` |  |
| R0814 | `CH3(CH2)9CH* + C* ⇌ CH3(CH2)9CHC* + *` | 与表面 C 偶联 |
| R0815 | `CH3(CH2)9CH* + CH* ⇌ CH3(CH2)9CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0816 | `CH3(CH2)9CH* + CH2* ⇌ CH3(CH2)9CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0817 | `CH3(CH2)9CH* + CH3* ⇌ CH3(CH2)9CHCH3* + *` |  |
| R0818 | `CH3(CH2)9CH2* + C* ⇌ CH3(CH2)10C* + *` | 与表面 C 偶联 |
| R0819 | `CH3(CH2)9CH2* + CH* ⇌ CH3(CH2)10CH* + *` |  |
| R0820 | `CH3(CH2)9CH2* + CH2* ⇌ CH3(CH2)10CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0821 | `CH3(CH2)9CH2* + CH3* ⇌ CH3(CH2)10CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0822 | `CH3(CH2)9C* + CO* ⇌ CH3(CH2)9CCO* + *` |  |
| R0823 | `CH3(CH2)9CH* + CO* ⇌ CH3(CH2)9CHCO* + *` |  |
| R0824 | `CH3(CH2)9CH2* + CO* ⇌ CH3(CH2)9CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0825 | `CH3(CH2)9CCO* + H* ⇌ CH3(CH2)9CCHO* + *` | C 端加氢 |
| R0826 | `CH3(CH2)9CCHO* + H* ⇌ CH3(CH2)9CCH2O* + *` | C 端加氢 |
| R0827 | `CH3(CH2)9CCOH* + H* ⇌ CH3(CH2)9CCHOH* + *` | C 端加氢 |
| R0828 | `CH3(CH2)9CCHOH* + H* ⇌ CH3(CH2)9CCH2OH* + *` | C 端加氢 |
| R0829 | `CH3(CH2)9CCO* + H* ⇌ CH3(CH2)9CCOH* + *` | O 端加氢 |
| R0830 | `CH3(CH2)9CCHO* + H* ⇌ CH3(CH2)9CCHOH* + *` | O 端加氢 |
| R0831 | `CH3(CH2)9CCH2O* + H* ⇌ CH3(CH2)9CCH2OH* + *` | O 端加氢 |
| R0832 | `CH3(CH2)9CHCO* + H* ⇌ CH3(CH2)9CHCHO* + *` | C 端加氢 |
| R0833 | `CH3(CH2)9CHCHO* + H* ⇌ CH3(CH2)9CHCH2O* + *` | C 端加氢 |
| R0834 | `CH3(CH2)9CHCOH* + H* ⇌ CH3(CH2)9CHCHOH* + *` | C 端加氢 |
| R0835 | `CH3(CH2)9CHCHOH* + H* ⇌ CH3(CH2)9CHCH2OH* + *` | C 端加氢 |
| R0836 | `CH3(CH2)9CHCO* + H* ⇌ CH3(CH2)9CHCOH* + *` | O 端加氢 |
| R0837 | `CH3(CH2)9CHCHO* + H* ⇌ CH3(CH2)9CHCHOH* + *` | O 端加氢 |
| R0838 | `CH3(CH2)9CHCH2O* + H* ⇌ CH3(CH2)9CHCH2OH* + *` | O 端加氢 |
| R0839 | `CH3(CH2)9CH2CO* + H* ⇌ CH3(CH2)9CH2CHO* + *` | C 端加氢 |
| R0840 | `CH3(CH2)9CH2CHO* + H* ⇌ CH3(CH2)9CH2CH2O* + *` | C 端加氢 |
| R0841 | `CH3(CH2)9CH2COH* + H* ⇌ CH3(CH2)9CH2CHOH* + *` | C 端加氢 |
| R0842 | `CH3(CH2)9CH2CHOH* + H* ⇌ CH3(CH2)9CH2CH2OH* + *` | C 端加氢 |
| R0843 | `CH3(CH2)9CH2CO* + H* ⇌ CH3(CH2)9CH2COH* + *` | O 端加氢 |
| R0844 | `CH3(CH2)9CH2CHO* + H* ⇌ CH3(CH2)9CH2CHOH* + *` | O 端加氢 |
| R0845 | `CH3(CH2)9CH2CH2O* + H* ⇌ CH3(CH2)9CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0846 | `CH3(CH2)9CCO* + * ⇌ CH3(CH2)9CC* + O*` |  |
| R0847 | `CH3(CH2)9CCOH* + * ⇌ CH3(CH2)9CC* + OH*` |  |
| R0848 | `CH3(CH2)9CCHO* + * ⇌ CH3(CH2)9CCH* + O*` |  |
| R0849 | `CH3(CH2)9CCHOH* + * ⇌ CH3(CH2)9CCH* + OH*` |  |
| R0850 | `CH3(CH2)9CCH2O* + * ⇌ CH3(CH2)9CCH2* + O*` |  |
| R0851 | `CH3(CH2)9CCH2OH* + * ⇌ CH3(CH2)9CCH2* + OH*` |  |
| R0852 | `CH3(CH2)9CHCO* + * ⇌ CH3(CH2)9CHC* + O*` |  |
| R0853 | `CH3(CH2)9CHCOH* + * ⇌ CH3(CH2)9CHC* + OH*` |  |
| R0854 | `CH3(CH2)9CHCHO* + * ⇌ CH3(CH2)9CHCH* + O*` |  |
| R0855 | `CH3(CH2)9CHCHOH* + * ⇌ CH3(CH2)9CHCH* + OH*` |  |
| R0856 | `CH3(CH2)9CHCH2O* + * ⇌ CH3(CH2)9CHCH2* + O*` |  |
| R0857 | `CH3(CH2)9CHCH2OH* + * ⇌ CH3(CH2)9CHCH2* + OH*` |  |
| R0858 | `CH3(CH2)9CH2CO* + * ⇌ CH3(CH2)10C* + O*` |  |
| R0859 | `CH3(CH2)9CH2COH* + * ⇌ CH3(CH2)10C* + OH*` |  |
| R0860 | `CH3(CH2)9CH2CHO* + * ⇌ CH3(CH2)10CH* + O*` |  |
| R0861 | `CH3(CH2)9CH2CHOH* + * ⇌ CH3(CH2)10CH* + OH*` |  |
| R0862 | `CH3(CH2)9CH2CH2O* + * ⇌ CH3(CH2)10CH2* + O*` |  |
| R0863 | `CH3(CH2)9CH2CH2OH* + * ⇌ CH3(CH2)10CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0864 | `CH3(CH2)8CH2COH* + HCOH* ⇌ CH3(CH2)9CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0865 | `CH3(CH2)8CHCH* + CH2* ⇌ CH3(CH2)8CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0866 | `CH3(CH2)8CHCHCH2* ⇌ CH3(CH2)9CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0867 | `CH3(CH2)8CHCHCH2* + H* ⇌ CH3(CH2)9CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0868 | `CH3(CH2)9CC* + H* ⇌ CH3(CH2)9CCH* + *` | α-C 加氢 |
| R0869 | `CH3(CH2)9CC* + H* ⇌ CH3(CH2)9CHC* + *` | β-C 加氢 |
| R0870 | `CH3(CH2)9CCH* + H* ⇌ CH3(CH2)9CCH2* + *` | α-C 加氢 |
| R0871 | `CH3(CH2)9CCH* + H* ⇌ CH3(CH2)9CHCH* + *` | β-C 加氢 |
| R0872 | `CH3(CH2)9CCH2* + H* ⇌ CH3(CH2)9CCH3* + *` | α-C 加氢 |
| R0873 | `CH3(CH2)9CCH2* + H* ⇌ CH3(CH2)9CHCH2* + *` | β-C 加氢 |
| R0874 | `CH3(CH2)9CCH3* + H* ⇌ CH3(CH2)9CHCH3* + *` | β-C 加氢 |
| R0875 | `CH3(CH2)9CHC* + H* ⇌ CH3(CH2)9CHCH* + *` | α-C 加氢 |
| R0876 | `CH3(CH2)9CHC* + H* ⇌ CH3(CH2)10C* + *` | β-C 加氢 |
| R0877 | `CH3(CH2)9CHCH* + H* ⇌ CH3(CH2)9CHCH2* + *` | α-C 加氢 |
| R0878 | `CH3(CH2)9CHCH* + H* ⇌ CH3(CH2)10CH* + *` | β-C 加氢 |
| R0879 | `CH3(CH2)9CHCH2* + H* ⇌ CH3(CH2)9CHCH3* + *` | α-C 加氢 |
| R0880 | `CH3(CH2)9CHCH2* + H* ⇌ CH3(CH2)10CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0881 | `CH3(CH2)9CHCH3* + H* ⇌ CH3(CH2)10CH3* + *` | β-C 加氢 |
| R0882 | `CH3(CH2)10C* + H* ⇌ CH3(CH2)10CH* + *` | α-C 加氢 |
| R0883 | `CH3(CH2)10CH* + H* ⇌ CH3(CH2)10CH2* + *` | α-C 加氢 |
| R0884 | `CH3(CH2)10CH2* + H* ⇌ CH3(CH2)10CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0885 | `CH3(CH2)10CH3* ⇌ C12H26(g) + *` | 正构烷烃脱附 |
| R0886 | `CH3(CH2)9CHCH2* ⇌ C12H24(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C12 → C13

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0887 | `CH3(CH2)10C* + C* ⇌ CH3(CH2)10CC* + *` | 与表面 C 偶联 |
| R0888 | `CH3(CH2)10C* + CH* ⇌ CH3(CH2)10CCH* + *` | 次烷基 + CH |
| R0889 | `CH3(CH2)10C* + CH2* ⇌ CH3(CH2)10CCH2* + *` |  |
| R0890 | `CH3(CH2)10C* + CH3* ⇌ CH3(CH2)10CCH3* + *` |  |
| R0891 | `CH3(CH2)10CH* + C* ⇌ CH3(CH2)10CHC* + *` | 与表面 C 偶联 |
| R0892 | `CH3(CH2)10CH* + CH* ⇌ CH3(CH2)10CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0893 | `CH3(CH2)10CH* + CH2* ⇌ CH3(CH2)10CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0894 | `CH3(CH2)10CH* + CH3* ⇌ CH3(CH2)10CHCH3* + *` |  |
| R0895 | `CH3(CH2)10CH2* + C* ⇌ CH3(CH2)11C* + *` | 与表面 C 偶联 |
| R0896 | `CH3(CH2)10CH2* + CH* ⇌ CH3(CH2)11CH* + *` |  |
| R0897 | `CH3(CH2)10CH2* + CH2* ⇌ CH3(CH2)11CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0898 | `CH3(CH2)10CH2* + CH3* ⇌ CH3(CH2)11CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0899 | `CH3(CH2)10C* + CO* ⇌ CH3(CH2)10CCO* + *` |  |
| R0900 | `CH3(CH2)10CH* + CO* ⇌ CH3(CH2)10CHCO* + *` |  |
| R0901 | `CH3(CH2)10CH2* + CO* ⇌ CH3(CH2)10CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0902 | `CH3(CH2)10CCO* + H* ⇌ CH3(CH2)10CCHO* + *` | C 端加氢 |
| R0903 | `CH3(CH2)10CCHO* + H* ⇌ CH3(CH2)10CCH2O* + *` | C 端加氢 |
| R0904 | `CH3(CH2)10CCOH* + H* ⇌ CH3(CH2)10CCHOH* + *` | C 端加氢 |
| R0905 | `CH3(CH2)10CCHOH* + H* ⇌ CH3(CH2)10CCH2OH* + *` | C 端加氢 |
| R0906 | `CH3(CH2)10CCO* + H* ⇌ CH3(CH2)10CCOH* + *` | O 端加氢 |
| R0907 | `CH3(CH2)10CCHO* + H* ⇌ CH3(CH2)10CCHOH* + *` | O 端加氢 |
| R0908 | `CH3(CH2)10CCH2O* + H* ⇌ CH3(CH2)10CCH2OH* + *` | O 端加氢 |
| R0909 | `CH3(CH2)10CHCO* + H* ⇌ CH3(CH2)10CHCHO* + *` | C 端加氢 |
| R0910 | `CH3(CH2)10CHCHO* + H* ⇌ CH3(CH2)10CHCH2O* + *` | C 端加氢 |
| R0911 | `CH3(CH2)10CHCOH* + H* ⇌ CH3(CH2)10CHCHOH* + *` | C 端加氢 |
| R0912 | `CH3(CH2)10CHCHOH* + H* ⇌ CH3(CH2)10CHCH2OH* + *` | C 端加氢 |
| R0913 | `CH3(CH2)10CHCO* + H* ⇌ CH3(CH2)10CHCOH* + *` | O 端加氢 |
| R0914 | `CH3(CH2)10CHCHO* + H* ⇌ CH3(CH2)10CHCHOH* + *` | O 端加氢 |
| R0915 | `CH3(CH2)10CHCH2O* + H* ⇌ CH3(CH2)10CHCH2OH* + *` | O 端加氢 |
| R0916 | `CH3(CH2)10CH2CO* + H* ⇌ CH3(CH2)10CH2CHO* + *` | C 端加氢 |
| R0917 | `CH3(CH2)10CH2CHO* + H* ⇌ CH3(CH2)10CH2CH2O* + *` | C 端加氢 |
| R0918 | `CH3(CH2)10CH2COH* + H* ⇌ CH3(CH2)10CH2CHOH* + *` | C 端加氢 |
| R0919 | `CH3(CH2)10CH2CHOH* + H* ⇌ CH3(CH2)10CH2CH2OH* + *` | C 端加氢 |
| R0920 | `CH3(CH2)10CH2CO* + H* ⇌ CH3(CH2)10CH2COH* + *` | O 端加氢 |
| R0921 | `CH3(CH2)10CH2CHO* + H* ⇌ CH3(CH2)10CH2CHOH* + *` | O 端加氢 |
| R0922 | `CH3(CH2)10CH2CH2O* + H* ⇌ CH3(CH2)10CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0923 | `CH3(CH2)10CCO* + * ⇌ CH3(CH2)10CC* + O*` |  |
| R0924 | `CH3(CH2)10CCOH* + * ⇌ CH3(CH2)10CC* + OH*` |  |
| R0925 | `CH3(CH2)10CCHO* + * ⇌ CH3(CH2)10CCH* + O*` |  |
| R0926 | `CH3(CH2)10CCHOH* + * ⇌ CH3(CH2)10CCH* + OH*` |  |
| R0927 | `CH3(CH2)10CCH2O* + * ⇌ CH3(CH2)10CCH2* + O*` |  |
| R0928 | `CH3(CH2)10CCH2OH* + * ⇌ CH3(CH2)10CCH2* + OH*` |  |
| R0929 | `CH3(CH2)10CHCO* + * ⇌ CH3(CH2)10CHC* + O*` |  |
| R0930 | `CH3(CH2)10CHCOH* + * ⇌ CH3(CH2)10CHC* + OH*` |  |
| R0931 | `CH3(CH2)10CHCHO* + * ⇌ CH3(CH2)10CHCH* + O*` |  |
| R0932 | `CH3(CH2)10CHCHOH* + * ⇌ CH3(CH2)10CHCH* + OH*` |  |
| R0933 | `CH3(CH2)10CHCH2O* + * ⇌ CH3(CH2)10CHCH2* + O*` |  |
| R0934 | `CH3(CH2)10CHCH2OH* + * ⇌ CH3(CH2)10CHCH2* + OH*` |  |
| R0935 | `CH3(CH2)10CH2CO* + * ⇌ CH3(CH2)11C* + O*` |  |
| R0936 | `CH3(CH2)10CH2COH* + * ⇌ CH3(CH2)11C* + OH*` |  |
| R0937 | `CH3(CH2)10CH2CHO* + * ⇌ CH3(CH2)11CH* + O*` |  |
| R0938 | `CH3(CH2)10CH2CHOH* + * ⇌ CH3(CH2)11CH* + OH*` |  |
| R0939 | `CH3(CH2)10CH2CH2O* + * ⇌ CH3(CH2)11CH2* + O*` |  |
| R0940 | `CH3(CH2)10CH2CH2OH* + * ⇌ CH3(CH2)11CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0941 | `CH3(CH2)9CH2COH* + HCOH* ⇌ CH3(CH2)10CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0942 | `CH3(CH2)9CHCH* + CH2* ⇌ CH3(CH2)9CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R0943 | `CH3(CH2)9CHCHCH2* ⇌ CH3(CH2)10CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R0944 | `CH3(CH2)9CHCHCH2* + H* ⇌ CH3(CH2)10CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0945 | `CH3(CH2)10CC* + H* ⇌ CH3(CH2)10CCH* + *` | α-C 加氢 |
| R0946 | `CH3(CH2)10CC* + H* ⇌ CH3(CH2)10CHC* + *` | β-C 加氢 |
| R0947 | `CH3(CH2)10CCH* + H* ⇌ CH3(CH2)10CCH2* + *` | α-C 加氢 |
| R0948 | `CH3(CH2)10CCH* + H* ⇌ CH3(CH2)10CHCH* + *` | β-C 加氢 |
| R0949 | `CH3(CH2)10CCH2* + H* ⇌ CH3(CH2)10CCH3* + *` | α-C 加氢 |
| R0950 | `CH3(CH2)10CCH2* + H* ⇌ CH3(CH2)10CHCH2* + *` | β-C 加氢 |
| R0951 | `CH3(CH2)10CCH3* + H* ⇌ CH3(CH2)10CHCH3* + *` | β-C 加氢 |
| R0952 | `CH3(CH2)10CHC* + H* ⇌ CH3(CH2)10CHCH* + *` | α-C 加氢 |
| R0953 | `CH3(CH2)10CHC* + H* ⇌ CH3(CH2)11C* + *` | β-C 加氢 |
| R0954 | `CH3(CH2)10CHCH* + H* ⇌ CH3(CH2)10CHCH2* + *` | α-C 加氢 |
| R0955 | `CH3(CH2)10CHCH* + H* ⇌ CH3(CH2)11CH* + *` | β-C 加氢 |
| R0956 | `CH3(CH2)10CHCH2* + H* ⇌ CH3(CH2)10CHCH3* + *` | α-C 加氢 |
| R0957 | `CH3(CH2)10CHCH2* + H* ⇌ CH3(CH2)11CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R0958 | `CH3(CH2)10CHCH3* + H* ⇌ CH3(CH2)11CH3* + *` | β-C 加氢 |
| R0959 | `CH3(CH2)11C* + H* ⇌ CH3(CH2)11CH* + *` | α-C 加氢 |
| R0960 | `CH3(CH2)11CH* + H* ⇌ CH3(CH2)11CH2* + *` | α-C 加氢 |
| R0961 | `CH3(CH2)11CH2* + H* ⇌ CH3(CH2)11CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0962 | `CH3(CH2)11CH3* ⇌ C13H28(g) + *` | 正构烷烃脱附 |
| R0963 | `CH3(CH2)10CHCH2* ⇌ C13H26(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C13 → C14

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0964 | `CH3(CH2)11C* + C* ⇌ CH3(CH2)11CC* + *` | 与表面 C 偶联 |
| R0965 | `CH3(CH2)11C* + CH* ⇌ CH3(CH2)11CCH* + *` | 次烷基 + CH |
| R0966 | `CH3(CH2)11C* + CH2* ⇌ CH3(CH2)11CCH2* + *` |  |
| R0967 | `CH3(CH2)11C* + CH3* ⇌ CH3(CH2)11CCH3* + *` |  |
| R0968 | `CH3(CH2)11CH* + C* ⇌ CH3(CH2)11CHC* + *` | 与表面 C 偶联 |
| R0969 | `CH3(CH2)11CH* + CH* ⇌ CH3(CH2)11CHCH* + *` | 亚烷基 + CH → 烯基 |
| R0970 | `CH3(CH2)11CH* + CH2* ⇌ CH3(CH2)11CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R0971 | `CH3(CH2)11CH* + CH3* ⇌ CH3(CH2)11CHCH3* + *` |  |
| R0972 | `CH3(CH2)11CH2* + C* ⇌ CH3(CH2)12C* + *` | 与表面 C 偶联 |
| R0973 | `CH3(CH2)11CH2* + CH* ⇌ CH3(CH2)12CH* + *` |  |
| R0974 | `CH3(CH2)11CH2* + CH2* ⇌ CH3(CH2)12CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R0975 | `CH3(CH2)11CH2* + CH3* ⇌ CH3(CH2)12CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0976 | `CH3(CH2)11C* + CO* ⇌ CH3(CH2)11CCO* + *` |  |
| R0977 | `CH3(CH2)11CH* + CO* ⇌ CH3(CH2)11CHCO* + *` |  |
| R0978 | `CH3(CH2)11CH2* + CO* ⇌ CH3(CH2)11CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R0979 | `CH3(CH2)11CCO* + H* ⇌ CH3(CH2)11CCHO* + *` | C 端加氢 |
| R0980 | `CH3(CH2)11CCHO* + H* ⇌ CH3(CH2)11CCH2O* + *` | C 端加氢 |
| R0981 | `CH3(CH2)11CCOH* + H* ⇌ CH3(CH2)11CCHOH* + *` | C 端加氢 |
| R0982 | `CH3(CH2)11CCHOH* + H* ⇌ CH3(CH2)11CCH2OH* + *` | C 端加氢 |
| R0983 | `CH3(CH2)11CCO* + H* ⇌ CH3(CH2)11CCOH* + *` | O 端加氢 |
| R0984 | `CH3(CH2)11CCHO* + H* ⇌ CH3(CH2)11CCHOH* + *` | O 端加氢 |
| R0985 | `CH3(CH2)11CCH2O* + H* ⇌ CH3(CH2)11CCH2OH* + *` | O 端加氢 |
| R0986 | `CH3(CH2)11CHCO* + H* ⇌ CH3(CH2)11CHCHO* + *` | C 端加氢 |
| R0987 | `CH3(CH2)11CHCHO* + H* ⇌ CH3(CH2)11CHCH2O* + *` | C 端加氢 |
| R0988 | `CH3(CH2)11CHCOH* + H* ⇌ CH3(CH2)11CHCHOH* + *` | C 端加氢 |
| R0989 | `CH3(CH2)11CHCHOH* + H* ⇌ CH3(CH2)11CHCH2OH* + *` | C 端加氢 |
| R0990 | `CH3(CH2)11CHCO* + H* ⇌ CH3(CH2)11CHCOH* + *` | O 端加氢 |
| R0991 | `CH3(CH2)11CHCHO* + H* ⇌ CH3(CH2)11CHCHOH* + *` | O 端加氢 |
| R0992 | `CH3(CH2)11CHCH2O* + H* ⇌ CH3(CH2)11CHCH2OH* + *` | O 端加氢 |
| R0993 | `CH3(CH2)11CH2CO* + H* ⇌ CH3(CH2)11CH2CHO* + *` | C 端加氢 |
| R0994 | `CH3(CH2)11CH2CHO* + H* ⇌ CH3(CH2)11CH2CH2O* + *` | C 端加氢 |
| R0995 | `CH3(CH2)11CH2COH* + H* ⇌ CH3(CH2)11CH2CHOH* + *` | C 端加氢 |
| R0996 | `CH3(CH2)11CH2CHOH* + H* ⇌ CH3(CH2)11CH2CH2OH* + *` | C 端加氢 |
| R0997 | `CH3(CH2)11CH2CO* + H* ⇌ CH3(CH2)11CH2COH* + *` | O 端加氢 |
| R0998 | `CH3(CH2)11CH2CHO* + H* ⇌ CH3(CH2)11CH2CHOH* + *` | O 端加氢 |
| R0999 | `CH3(CH2)11CH2CH2O* + H* ⇌ CH3(CH2)11CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1000 | `CH3(CH2)11CCO* + * ⇌ CH3(CH2)11CC* + O*` |  |
| R1001 | `CH3(CH2)11CCOH* + * ⇌ CH3(CH2)11CC* + OH*` |  |
| R1002 | `CH3(CH2)11CCHO* + * ⇌ CH3(CH2)11CCH* + O*` |  |
| R1003 | `CH3(CH2)11CCHOH* + * ⇌ CH3(CH2)11CCH* + OH*` |  |
| R1004 | `CH3(CH2)11CCH2O* + * ⇌ CH3(CH2)11CCH2* + O*` |  |
| R1005 | `CH3(CH2)11CCH2OH* + * ⇌ CH3(CH2)11CCH2* + OH*` |  |
| R1006 | `CH3(CH2)11CHCO* + * ⇌ CH3(CH2)11CHC* + O*` |  |
| R1007 | `CH3(CH2)11CHCOH* + * ⇌ CH3(CH2)11CHC* + OH*` |  |
| R1008 | `CH3(CH2)11CHCHO* + * ⇌ CH3(CH2)11CHCH* + O*` |  |
| R1009 | `CH3(CH2)11CHCHOH* + * ⇌ CH3(CH2)11CHCH* + OH*` |  |
| R1010 | `CH3(CH2)11CHCH2O* + * ⇌ CH3(CH2)11CHCH2* + O*` |  |
| R1011 | `CH3(CH2)11CHCH2OH* + * ⇌ CH3(CH2)11CHCH2* + OH*` |  |
| R1012 | `CH3(CH2)11CH2CO* + * ⇌ CH3(CH2)12C* + O*` |  |
| R1013 | `CH3(CH2)11CH2COH* + * ⇌ CH3(CH2)12C* + OH*` |  |
| R1014 | `CH3(CH2)11CH2CHO* + * ⇌ CH3(CH2)12CH* + O*` |  |
| R1015 | `CH3(CH2)11CH2CHOH* + * ⇌ CH3(CH2)12CH* + OH*` |  |
| R1016 | `CH3(CH2)11CH2CH2O* + * ⇌ CH3(CH2)12CH2* + O*` |  |
| R1017 | `CH3(CH2)11CH2CH2OH* + * ⇌ CH3(CH2)12CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1018 | `CH3(CH2)10CH2COH* + HCOH* ⇌ CH3(CH2)11CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1019 | `CH3(CH2)10CHCH* + CH2* ⇌ CH3(CH2)10CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R1020 | `CH3(CH2)10CHCHCH2* ⇌ CH3(CH2)11CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R1021 | `CH3(CH2)10CHCHCH2* + H* ⇌ CH3(CH2)11CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1022 | `CH3(CH2)11CC* + H* ⇌ CH3(CH2)11CCH* + *` | α-C 加氢 |
| R1023 | `CH3(CH2)11CC* + H* ⇌ CH3(CH2)11CHC* + *` | β-C 加氢 |
| R1024 | `CH3(CH2)11CCH* + H* ⇌ CH3(CH2)11CCH2* + *` | α-C 加氢 |
| R1025 | `CH3(CH2)11CCH* + H* ⇌ CH3(CH2)11CHCH* + *` | β-C 加氢 |
| R1026 | `CH3(CH2)11CCH2* + H* ⇌ CH3(CH2)11CCH3* + *` | α-C 加氢 |
| R1027 | `CH3(CH2)11CCH2* + H* ⇌ CH3(CH2)11CHCH2* + *` | β-C 加氢 |
| R1028 | `CH3(CH2)11CCH3* + H* ⇌ CH3(CH2)11CHCH3* + *` | β-C 加氢 |
| R1029 | `CH3(CH2)11CHC* + H* ⇌ CH3(CH2)11CHCH* + *` | α-C 加氢 |
| R1030 | `CH3(CH2)11CHC* + H* ⇌ CH3(CH2)12C* + *` | β-C 加氢 |
| R1031 | `CH3(CH2)11CHCH* + H* ⇌ CH3(CH2)11CHCH2* + *` | α-C 加氢 |
| R1032 | `CH3(CH2)11CHCH* + H* ⇌ CH3(CH2)12CH* + *` | β-C 加氢 |
| R1033 | `CH3(CH2)11CHCH2* + H* ⇌ CH3(CH2)11CHCH3* + *` | α-C 加氢 |
| R1034 | `CH3(CH2)11CHCH2* + H* ⇌ CH3(CH2)12CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R1035 | `CH3(CH2)11CHCH3* + H* ⇌ CH3(CH2)12CH3* + *` | β-C 加氢 |
| R1036 | `CH3(CH2)12C* + H* ⇌ CH3(CH2)12CH* + *` | α-C 加氢 |
| R1037 | `CH3(CH2)12CH* + H* ⇌ CH3(CH2)12CH2* + *` | α-C 加氢 |
| R1038 | `CH3(CH2)12CH2* + H* ⇌ CH3(CH2)12CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1039 | `CH3(CH2)12CH3* ⇌ C14H30(g) + *` | 正构烷烃脱附 |
| R1040 | `CH3(CH2)11CHCH2* ⇌ C14H28(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C14 → C15

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1041 | `CH3(CH2)12C* + C* ⇌ CH3(CH2)12CC* + *` | 与表面 C 偶联 |
| R1042 | `CH3(CH2)12C* + CH* ⇌ CH3(CH2)12CCH* + *` | 次烷基 + CH |
| R1043 | `CH3(CH2)12C* + CH2* ⇌ CH3(CH2)12CCH2* + *` |  |
| R1044 | `CH3(CH2)12C* + CH3* ⇌ CH3(CH2)12CCH3* + *` |  |
| R1045 | `CH3(CH2)12CH* + C* ⇌ CH3(CH2)12CHC* + *` | 与表面 C 偶联 |
| R1046 | `CH3(CH2)12CH* + CH* ⇌ CH3(CH2)12CHCH* + *` | 亚烷基 + CH → 烯基 |
| R1047 | `CH3(CH2)12CH* + CH2* ⇌ CH3(CH2)12CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R1048 | `CH3(CH2)12CH* + CH3* ⇌ CH3(CH2)12CHCH3* + *` |  |
| R1049 | `CH3(CH2)12CH2* + C* ⇌ CH3(CH2)13C* + *` | 与表面 C 偶联 |
| R1050 | `CH3(CH2)12CH2* + CH* ⇌ CH3(CH2)13CH* + *` |  |
| R1051 | `CH3(CH2)12CH2* + CH2* ⇌ CH3(CH2)13CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R1052 | `CH3(CH2)12CH2* + CH3* ⇌ CH3(CH2)13CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1053 | `CH3(CH2)12C* + CO* ⇌ CH3(CH2)12CCO* + *` |  |
| R1054 | `CH3(CH2)12CH* + CO* ⇌ CH3(CH2)12CHCO* + *` |  |
| R1055 | `CH3(CH2)12CH2* + CO* ⇌ CH3(CH2)12CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1056 | `CH3(CH2)12CCO* + H* ⇌ CH3(CH2)12CCHO* + *` | C 端加氢 |
| R1057 | `CH3(CH2)12CCHO* + H* ⇌ CH3(CH2)12CCH2O* + *` | C 端加氢 |
| R1058 | `CH3(CH2)12CCOH* + H* ⇌ CH3(CH2)12CCHOH* + *` | C 端加氢 |
| R1059 | `CH3(CH2)12CCHOH* + H* ⇌ CH3(CH2)12CCH2OH* + *` | C 端加氢 |
| R1060 | `CH3(CH2)12CCO* + H* ⇌ CH3(CH2)12CCOH* + *` | O 端加氢 |
| R1061 | `CH3(CH2)12CCHO* + H* ⇌ CH3(CH2)12CCHOH* + *` | O 端加氢 |
| R1062 | `CH3(CH2)12CCH2O* + H* ⇌ CH3(CH2)12CCH2OH* + *` | O 端加氢 |
| R1063 | `CH3(CH2)12CHCO* + H* ⇌ CH3(CH2)12CHCHO* + *` | C 端加氢 |
| R1064 | `CH3(CH2)12CHCHO* + H* ⇌ CH3(CH2)12CHCH2O* + *` | C 端加氢 |
| R1065 | `CH3(CH2)12CHCOH* + H* ⇌ CH3(CH2)12CHCHOH* + *` | C 端加氢 |
| R1066 | `CH3(CH2)12CHCHOH* + H* ⇌ CH3(CH2)12CHCH2OH* + *` | C 端加氢 |
| R1067 | `CH3(CH2)12CHCO* + H* ⇌ CH3(CH2)12CHCOH* + *` | O 端加氢 |
| R1068 | `CH3(CH2)12CHCHO* + H* ⇌ CH3(CH2)12CHCHOH* + *` | O 端加氢 |
| R1069 | `CH3(CH2)12CHCH2O* + H* ⇌ CH3(CH2)12CHCH2OH* + *` | O 端加氢 |
| R1070 | `CH3(CH2)12CH2CO* + H* ⇌ CH3(CH2)12CH2CHO* + *` | C 端加氢 |
| R1071 | `CH3(CH2)12CH2CHO* + H* ⇌ CH3(CH2)12CH2CH2O* + *` | C 端加氢 |
| R1072 | `CH3(CH2)12CH2COH* + H* ⇌ CH3(CH2)12CH2CHOH* + *` | C 端加氢 |
| R1073 | `CH3(CH2)12CH2CHOH* + H* ⇌ CH3(CH2)12CH2CH2OH* + *` | C 端加氢 |
| R1074 | `CH3(CH2)12CH2CO* + H* ⇌ CH3(CH2)12CH2COH* + *` | O 端加氢 |
| R1075 | `CH3(CH2)12CH2CHO* + H* ⇌ CH3(CH2)12CH2CHOH* + *` | O 端加氢 |
| R1076 | `CH3(CH2)12CH2CH2O* + H* ⇌ CH3(CH2)12CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1077 | `CH3(CH2)12CCO* + * ⇌ CH3(CH2)12CC* + O*` |  |
| R1078 | `CH3(CH2)12CCOH* + * ⇌ CH3(CH2)12CC* + OH*` |  |
| R1079 | `CH3(CH2)12CCHO* + * ⇌ CH3(CH2)12CCH* + O*` |  |
| R1080 | `CH3(CH2)12CCHOH* + * ⇌ CH3(CH2)12CCH* + OH*` |  |
| R1081 | `CH3(CH2)12CCH2O* + * ⇌ CH3(CH2)12CCH2* + O*` |  |
| R1082 | `CH3(CH2)12CCH2OH* + * ⇌ CH3(CH2)12CCH2* + OH*` |  |
| R1083 | `CH3(CH2)12CHCO* + * ⇌ CH3(CH2)12CHC* + O*` |  |
| R1084 | `CH3(CH2)12CHCOH* + * ⇌ CH3(CH2)12CHC* + OH*` |  |
| R1085 | `CH3(CH2)12CHCHO* + * ⇌ CH3(CH2)12CHCH* + O*` |  |
| R1086 | `CH3(CH2)12CHCHOH* + * ⇌ CH3(CH2)12CHCH* + OH*` |  |
| R1087 | `CH3(CH2)12CHCH2O* + * ⇌ CH3(CH2)12CHCH2* + O*` |  |
| R1088 | `CH3(CH2)12CHCH2OH* + * ⇌ CH3(CH2)12CHCH2* + OH*` |  |
| R1089 | `CH3(CH2)12CH2CO* + * ⇌ CH3(CH2)13C* + O*` |  |
| R1090 | `CH3(CH2)12CH2COH* + * ⇌ CH3(CH2)13C* + OH*` |  |
| R1091 | `CH3(CH2)12CH2CHO* + * ⇌ CH3(CH2)13CH* + O*` |  |
| R1092 | `CH3(CH2)12CH2CHOH* + * ⇌ CH3(CH2)13CH* + OH*` |  |
| R1093 | `CH3(CH2)12CH2CH2O* + * ⇌ CH3(CH2)13CH2* + O*` |  |
| R1094 | `CH3(CH2)12CH2CH2OH* + * ⇌ CH3(CH2)13CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1095 | `CH3(CH2)11CH2COH* + HCOH* ⇌ CH3(CH2)12CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1096 | `CH3(CH2)11CHCH* + CH2* ⇌ CH3(CH2)11CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R1097 | `CH3(CH2)11CHCHCH2* ⇌ CH3(CH2)12CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R1098 | `CH3(CH2)11CHCHCH2* + H* ⇌ CH3(CH2)12CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1099 | `CH3(CH2)12CC* + H* ⇌ CH3(CH2)12CCH* + *` | α-C 加氢 |
| R1100 | `CH3(CH2)12CC* + H* ⇌ CH3(CH2)12CHC* + *` | β-C 加氢 |
| R1101 | `CH3(CH2)12CCH* + H* ⇌ CH3(CH2)12CCH2* + *` | α-C 加氢 |
| R1102 | `CH3(CH2)12CCH* + H* ⇌ CH3(CH2)12CHCH* + *` | β-C 加氢 |
| R1103 | `CH3(CH2)12CCH2* + H* ⇌ CH3(CH2)12CCH3* + *` | α-C 加氢 |
| R1104 | `CH3(CH2)12CCH2* + H* ⇌ CH3(CH2)12CHCH2* + *` | β-C 加氢 |
| R1105 | `CH3(CH2)12CCH3* + H* ⇌ CH3(CH2)12CHCH3* + *` | β-C 加氢 |
| R1106 | `CH3(CH2)12CHC* + H* ⇌ CH3(CH2)12CHCH* + *` | α-C 加氢 |
| R1107 | `CH3(CH2)12CHC* + H* ⇌ CH3(CH2)13C* + *` | β-C 加氢 |
| R1108 | `CH3(CH2)12CHCH* + H* ⇌ CH3(CH2)12CHCH2* + *` | α-C 加氢 |
| R1109 | `CH3(CH2)12CHCH* + H* ⇌ CH3(CH2)13CH* + *` | β-C 加氢 |
| R1110 | `CH3(CH2)12CHCH2* + H* ⇌ CH3(CH2)12CHCH3* + *` | α-C 加氢 |
| R1111 | `CH3(CH2)12CHCH2* + H* ⇌ CH3(CH2)13CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R1112 | `CH3(CH2)12CHCH3* + H* ⇌ CH3(CH2)13CH3* + *` | β-C 加氢 |
| R1113 | `CH3(CH2)13C* + H* ⇌ CH3(CH2)13CH* + *` | α-C 加氢 |
| R1114 | `CH3(CH2)13CH* + H* ⇌ CH3(CH2)13CH2* + *` | α-C 加氢 |
| R1115 | `CH3(CH2)13CH2* + H* ⇌ CH3(CH2)13CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1116 | `CH3(CH2)13CH3* ⇌ C15H32(g) + *` | 正构烷烃脱附 |
| R1117 | `CH3(CH2)12CHCH2* ⇌ C15H30(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C15 → C16

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1118 | `CH3(CH2)13C* + C* ⇌ CH3(CH2)13CC* + *` | 与表面 C 偶联 |
| R1119 | `CH3(CH2)13C* + CH* ⇌ CH3(CH2)13CCH* + *` | 次烷基 + CH |
| R1120 | `CH3(CH2)13C* + CH2* ⇌ CH3(CH2)13CCH2* + *` |  |
| R1121 | `CH3(CH2)13C* + CH3* ⇌ CH3(CH2)13CCH3* + *` |  |
| R1122 | `CH3(CH2)13CH* + C* ⇌ CH3(CH2)13CHC* + *` | 与表面 C 偶联 |
| R1123 | `CH3(CH2)13CH* + CH* ⇌ CH3(CH2)13CHCH* + *` | 亚烷基 + CH → 烯基 |
| R1124 | `CH3(CH2)13CH* + CH2* ⇌ CH3(CH2)13CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R1125 | `CH3(CH2)13CH* + CH3* ⇌ CH3(CH2)13CHCH3* + *` |  |
| R1126 | `CH3(CH2)13CH2* + C* ⇌ CH3(CH2)14C* + *` | 与表面 C 偶联 |
| R1127 | `CH3(CH2)13CH2* + CH* ⇌ CH3(CH2)14CH* + *` |  |
| R1128 | `CH3(CH2)13CH2* + CH2* ⇌ CH3(CH2)14CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R1129 | `CH3(CH2)13CH2* + CH3* ⇌ CH3(CH2)14CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1130 | `CH3(CH2)13C* + CO* ⇌ CH3(CH2)13CCO* + *` |  |
| R1131 | `CH3(CH2)13CH* + CO* ⇌ CH3(CH2)13CHCO* + *` |  |
| R1132 | `CH3(CH2)13CH2* + CO* ⇌ CH3(CH2)13CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1133 | `CH3(CH2)13CCO* + H* ⇌ CH3(CH2)13CCHO* + *` | C 端加氢 |
| R1134 | `CH3(CH2)13CCHO* + H* ⇌ CH3(CH2)13CCH2O* + *` | C 端加氢 |
| R1135 | `CH3(CH2)13CCOH* + H* ⇌ CH3(CH2)13CCHOH* + *` | C 端加氢 |
| R1136 | `CH3(CH2)13CCHOH* + H* ⇌ CH3(CH2)13CCH2OH* + *` | C 端加氢 |
| R1137 | `CH3(CH2)13CCO* + H* ⇌ CH3(CH2)13CCOH* + *` | O 端加氢 |
| R1138 | `CH3(CH2)13CCHO* + H* ⇌ CH3(CH2)13CCHOH* + *` | O 端加氢 |
| R1139 | `CH3(CH2)13CCH2O* + H* ⇌ CH3(CH2)13CCH2OH* + *` | O 端加氢 |
| R1140 | `CH3(CH2)13CHCO* + H* ⇌ CH3(CH2)13CHCHO* + *` | C 端加氢 |
| R1141 | `CH3(CH2)13CHCHO* + H* ⇌ CH3(CH2)13CHCH2O* + *` | C 端加氢 |
| R1142 | `CH3(CH2)13CHCOH* + H* ⇌ CH3(CH2)13CHCHOH* + *` | C 端加氢 |
| R1143 | `CH3(CH2)13CHCHOH* + H* ⇌ CH3(CH2)13CHCH2OH* + *` | C 端加氢 |
| R1144 | `CH3(CH2)13CHCO* + H* ⇌ CH3(CH2)13CHCOH* + *` | O 端加氢 |
| R1145 | `CH3(CH2)13CHCHO* + H* ⇌ CH3(CH2)13CHCHOH* + *` | O 端加氢 |
| R1146 | `CH3(CH2)13CHCH2O* + H* ⇌ CH3(CH2)13CHCH2OH* + *` | O 端加氢 |
| R1147 | `CH3(CH2)13CH2CO* + H* ⇌ CH3(CH2)13CH2CHO* + *` | C 端加氢 |
| R1148 | `CH3(CH2)13CH2CHO* + H* ⇌ CH3(CH2)13CH2CH2O* + *` | C 端加氢 |
| R1149 | `CH3(CH2)13CH2COH* + H* ⇌ CH3(CH2)13CH2CHOH* + *` | C 端加氢 |
| R1150 | `CH3(CH2)13CH2CHOH* + H* ⇌ CH3(CH2)13CH2CH2OH* + *` | C 端加氢 |
| R1151 | `CH3(CH2)13CH2CO* + H* ⇌ CH3(CH2)13CH2COH* + *` | O 端加氢 |
| R1152 | `CH3(CH2)13CH2CHO* + H* ⇌ CH3(CH2)13CH2CHOH* + *` | O 端加氢 |
| R1153 | `CH3(CH2)13CH2CH2O* + H* ⇌ CH3(CH2)13CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1154 | `CH3(CH2)13CCO* + * ⇌ CH3(CH2)13CC* + O*` |  |
| R1155 | `CH3(CH2)13CCOH* + * ⇌ CH3(CH2)13CC* + OH*` |  |
| R1156 | `CH3(CH2)13CCHO* + * ⇌ CH3(CH2)13CCH* + O*` |  |
| R1157 | `CH3(CH2)13CCHOH* + * ⇌ CH3(CH2)13CCH* + OH*` |  |
| R1158 | `CH3(CH2)13CCH2O* + * ⇌ CH3(CH2)13CCH2* + O*` |  |
| R1159 | `CH3(CH2)13CCH2OH* + * ⇌ CH3(CH2)13CCH2* + OH*` |  |
| R1160 | `CH3(CH2)13CHCO* + * ⇌ CH3(CH2)13CHC* + O*` |  |
| R1161 | `CH3(CH2)13CHCOH* + * ⇌ CH3(CH2)13CHC* + OH*` |  |
| R1162 | `CH3(CH2)13CHCHO* + * ⇌ CH3(CH2)13CHCH* + O*` |  |
| R1163 | `CH3(CH2)13CHCHOH* + * ⇌ CH3(CH2)13CHCH* + OH*` |  |
| R1164 | `CH3(CH2)13CHCH2O* + * ⇌ CH3(CH2)13CHCH2* + O*` |  |
| R1165 | `CH3(CH2)13CHCH2OH* + * ⇌ CH3(CH2)13CHCH2* + OH*` |  |
| R1166 | `CH3(CH2)13CH2CO* + * ⇌ CH3(CH2)14C* + O*` |  |
| R1167 | `CH3(CH2)13CH2COH* + * ⇌ CH3(CH2)14C* + OH*` |  |
| R1168 | `CH3(CH2)13CH2CHO* + * ⇌ CH3(CH2)14CH* + O*` |  |
| R1169 | `CH3(CH2)13CH2CHOH* + * ⇌ CH3(CH2)14CH* + OH*` |  |
| R1170 | `CH3(CH2)13CH2CH2O* + * ⇌ CH3(CH2)14CH2* + O*` |  |
| R1171 | `CH3(CH2)13CH2CH2OH* + * ⇌ CH3(CH2)14CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1172 | `CH3(CH2)12CH2COH* + HCOH* ⇌ CH3(CH2)13CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1173 | `CH3(CH2)12CHCH* + CH2* ⇌ CH3(CH2)12CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R1174 | `CH3(CH2)12CHCHCH2* ⇌ CH3(CH2)13CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R1175 | `CH3(CH2)12CHCHCH2* + H* ⇌ CH3(CH2)13CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1176 | `CH3(CH2)13CC* + H* ⇌ CH3(CH2)13CCH* + *` | α-C 加氢 |
| R1177 | `CH3(CH2)13CC* + H* ⇌ CH3(CH2)13CHC* + *` | β-C 加氢 |
| R1178 | `CH3(CH2)13CCH* + H* ⇌ CH3(CH2)13CCH2* + *` | α-C 加氢 |
| R1179 | `CH3(CH2)13CCH* + H* ⇌ CH3(CH2)13CHCH* + *` | β-C 加氢 |
| R1180 | `CH3(CH2)13CCH2* + H* ⇌ CH3(CH2)13CCH3* + *` | α-C 加氢 |
| R1181 | `CH3(CH2)13CCH2* + H* ⇌ CH3(CH2)13CHCH2* + *` | β-C 加氢 |
| R1182 | `CH3(CH2)13CCH3* + H* ⇌ CH3(CH2)13CHCH3* + *` | β-C 加氢 |
| R1183 | `CH3(CH2)13CHC* + H* ⇌ CH3(CH2)13CHCH* + *` | α-C 加氢 |
| R1184 | `CH3(CH2)13CHC* + H* ⇌ CH3(CH2)14C* + *` | β-C 加氢 |
| R1185 | `CH3(CH2)13CHCH* + H* ⇌ CH3(CH2)13CHCH2* + *` | α-C 加氢 |
| R1186 | `CH3(CH2)13CHCH* + H* ⇌ CH3(CH2)14CH* + *` | β-C 加氢 |
| R1187 | `CH3(CH2)13CHCH2* + H* ⇌ CH3(CH2)13CHCH3* + *` | α-C 加氢 |
| R1188 | `CH3(CH2)13CHCH2* + H* ⇌ CH3(CH2)14CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R1189 | `CH3(CH2)13CHCH3* + H* ⇌ CH3(CH2)14CH3* + *` | β-C 加氢 |
| R1190 | `CH3(CH2)14C* + H* ⇌ CH3(CH2)14CH* + *` | α-C 加氢 |
| R1191 | `CH3(CH2)14CH* + H* ⇌ CH3(CH2)14CH2* + *` | α-C 加氢 |
| R1192 | `CH3(CH2)14CH2* + H* ⇌ CH3(CH2)14CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1193 | `CH3(CH2)14CH3* ⇌ C16H34(g) + *` | 正构烷烃脱附 |
| R1194 | `CH3(CH2)13CHCH2* ⇌ C16H32(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C16 → C17

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1195 | `CH3(CH2)14C* + C* ⇌ CH3(CH2)14CC* + *` | 与表面 C 偶联 |
| R1196 | `CH3(CH2)14C* + CH* ⇌ CH3(CH2)14CCH* + *` | 次烷基 + CH |
| R1197 | `CH3(CH2)14C* + CH2* ⇌ CH3(CH2)14CCH2* + *` |  |
| R1198 | `CH3(CH2)14C* + CH3* ⇌ CH3(CH2)14CCH3* + *` |  |
| R1199 | `CH3(CH2)14CH* + C* ⇌ CH3(CH2)14CHC* + *` | 与表面 C 偶联 |
| R1200 | `CH3(CH2)14CH* + CH* ⇌ CH3(CH2)14CHCH* + *` | 亚烷基 + CH → 烯基 |
| R1201 | `CH3(CH2)14CH* + CH2* ⇌ CH3(CH2)14CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R1202 | `CH3(CH2)14CH* + CH3* ⇌ CH3(CH2)14CHCH3* + *` |  |
| R1203 | `CH3(CH2)14CH2* + C* ⇌ CH3(CH2)15C* + *` | 与表面 C 偶联 |
| R1204 | `CH3(CH2)14CH2* + CH* ⇌ CH3(CH2)15CH* + *` |  |
| R1205 | `CH3(CH2)14CH2* + CH2* ⇌ CH3(CH2)15CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R1206 | `CH3(CH2)14CH2* + CH3* ⇌ CH3(CH2)15CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1207 | `CH3(CH2)14C* + CO* ⇌ CH3(CH2)14CCO* + *` |  |
| R1208 | `CH3(CH2)14CH* + CO* ⇌ CH3(CH2)14CHCO* + *` |  |
| R1209 | `CH3(CH2)14CH2* + CO* ⇌ CH3(CH2)14CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1210 | `CH3(CH2)14CCO* + H* ⇌ CH3(CH2)14CCHO* + *` | C 端加氢 |
| R1211 | `CH3(CH2)14CCHO* + H* ⇌ CH3(CH2)14CCH2O* + *` | C 端加氢 |
| R1212 | `CH3(CH2)14CCOH* + H* ⇌ CH3(CH2)14CCHOH* + *` | C 端加氢 |
| R1213 | `CH3(CH2)14CCHOH* + H* ⇌ CH3(CH2)14CCH2OH* + *` | C 端加氢 |
| R1214 | `CH3(CH2)14CCO* + H* ⇌ CH3(CH2)14CCOH* + *` | O 端加氢 |
| R1215 | `CH3(CH2)14CCHO* + H* ⇌ CH3(CH2)14CCHOH* + *` | O 端加氢 |
| R1216 | `CH3(CH2)14CCH2O* + H* ⇌ CH3(CH2)14CCH2OH* + *` | O 端加氢 |
| R1217 | `CH3(CH2)14CHCO* + H* ⇌ CH3(CH2)14CHCHO* + *` | C 端加氢 |
| R1218 | `CH3(CH2)14CHCHO* + H* ⇌ CH3(CH2)14CHCH2O* + *` | C 端加氢 |
| R1219 | `CH3(CH2)14CHCOH* + H* ⇌ CH3(CH2)14CHCHOH* + *` | C 端加氢 |
| R1220 | `CH3(CH2)14CHCHOH* + H* ⇌ CH3(CH2)14CHCH2OH* + *` | C 端加氢 |
| R1221 | `CH3(CH2)14CHCO* + H* ⇌ CH3(CH2)14CHCOH* + *` | O 端加氢 |
| R1222 | `CH3(CH2)14CHCHO* + H* ⇌ CH3(CH2)14CHCHOH* + *` | O 端加氢 |
| R1223 | `CH3(CH2)14CHCH2O* + H* ⇌ CH3(CH2)14CHCH2OH* + *` | O 端加氢 |
| R1224 | `CH3(CH2)14CH2CO* + H* ⇌ CH3(CH2)14CH2CHO* + *` | C 端加氢 |
| R1225 | `CH3(CH2)14CH2CHO* + H* ⇌ CH3(CH2)14CH2CH2O* + *` | C 端加氢 |
| R1226 | `CH3(CH2)14CH2COH* + H* ⇌ CH3(CH2)14CH2CHOH* + *` | C 端加氢 |
| R1227 | `CH3(CH2)14CH2CHOH* + H* ⇌ CH3(CH2)14CH2CH2OH* + *` | C 端加氢 |
| R1228 | `CH3(CH2)14CH2CO* + H* ⇌ CH3(CH2)14CH2COH* + *` | O 端加氢 |
| R1229 | `CH3(CH2)14CH2CHO* + H* ⇌ CH3(CH2)14CH2CHOH* + *` | O 端加氢 |
| R1230 | `CH3(CH2)14CH2CH2O* + H* ⇌ CH3(CH2)14CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1231 | `CH3(CH2)14CCO* + * ⇌ CH3(CH2)14CC* + O*` |  |
| R1232 | `CH3(CH2)14CCOH* + * ⇌ CH3(CH2)14CC* + OH*` |  |
| R1233 | `CH3(CH2)14CCHO* + * ⇌ CH3(CH2)14CCH* + O*` |  |
| R1234 | `CH3(CH2)14CCHOH* + * ⇌ CH3(CH2)14CCH* + OH*` |  |
| R1235 | `CH3(CH2)14CCH2O* + * ⇌ CH3(CH2)14CCH2* + O*` |  |
| R1236 | `CH3(CH2)14CCH2OH* + * ⇌ CH3(CH2)14CCH2* + OH*` |  |
| R1237 | `CH3(CH2)14CHCO* + * ⇌ CH3(CH2)14CHC* + O*` |  |
| R1238 | `CH3(CH2)14CHCOH* + * ⇌ CH3(CH2)14CHC* + OH*` |  |
| R1239 | `CH3(CH2)14CHCHO* + * ⇌ CH3(CH2)14CHCH* + O*` |  |
| R1240 | `CH3(CH2)14CHCHOH* + * ⇌ CH3(CH2)14CHCH* + OH*` |  |
| R1241 | `CH3(CH2)14CHCH2O* + * ⇌ CH3(CH2)14CHCH2* + O*` |  |
| R1242 | `CH3(CH2)14CHCH2OH* + * ⇌ CH3(CH2)14CHCH2* + OH*` |  |
| R1243 | `CH3(CH2)14CH2CO* + * ⇌ CH3(CH2)15C* + O*` |  |
| R1244 | `CH3(CH2)14CH2COH* + * ⇌ CH3(CH2)15C* + OH*` |  |
| R1245 | `CH3(CH2)14CH2CHO* + * ⇌ CH3(CH2)15CH* + O*` |  |
| R1246 | `CH3(CH2)14CH2CHOH* + * ⇌ CH3(CH2)15CH* + OH*` |  |
| R1247 | `CH3(CH2)14CH2CH2O* + * ⇌ CH3(CH2)15CH2* + O*` |  |
| R1248 | `CH3(CH2)14CH2CH2OH* + * ⇌ CH3(CH2)15CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1249 | `CH3(CH2)13CH2COH* + HCOH* ⇌ CH3(CH2)14CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1250 | `CH3(CH2)13CHCH* + CH2* ⇌ CH3(CH2)13CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R1251 | `CH3(CH2)13CHCHCH2* ⇌ CH3(CH2)14CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R1252 | `CH3(CH2)13CHCHCH2* + H* ⇌ CH3(CH2)14CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1253 | `CH3(CH2)14CC* + H* ⇌ CH3(CH2)14CCH* + *` | α-C 加氢 |
| R1254 | `CH3(CH2)14CC* + H* ⇌ CH3(CH2)14CHC* + *` | β-C 加氢 |
| R1255 | `CH3(CH2)14CCH* + H* ⇌ CH3(CH2)14CCH2* + *` | α-C 加氢 |
| R1256 | `CH3(CH2)14CCH* + H* ⇌ CH3(CH2)14CHCH* + *` | β-C 加氢 |
| R1257 | `CH3(CH2)14CCH2* + H* ⇌ CH3(CH2)14CCH3* + *` | α-C 加氢 |
| R1258 | `CH3(CH2)14CCH2* + H* ⇌ CH3(CH2)14CHCH2* + *` | β-C 加氢 |
| R1259 | `CH3(CH2)14CCH3* + H* ⇌ CH3(CH2)14CHCH3* + *` | β-C 加氢 |
| R1260 | `CH3(CH2)14CHC* + H* ⇌ CH3(CH2)14CHCH* + *` | α-C 加氢 |
| R1261 | `CH3(CH2)14CHC* + H* ⇌ CH3(CH2)15C* + *` | β-C 加氢 |
| R1262 | `CH3(CH2)14CHCH* + H* ⇌ CH3(CH2)14CHCH2* + *` | α-C 加氢 |
| R1263 | `CH3(CH2)14CHCH* + H* ⇌ CH3(CH2)15CH* + *` | β-C 加氢 |
| R1264 | `CH3(CH2)14CHCH2* + H* ⇌ CH3(CH2)14CHCH3* + *` | α-C 加氢 |
| R1265 | `CH3(CH2)14CHCH2* + H* ⇌ CH3(CH2)15CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R1266 | `CH3(CH2)14CHCH3* + H* ⇌ CH3(CH2)15CH3* + *` | β-C 加氢 |
| R1267 | `CH3(CH2)15C* + H* ⇌ CH3(CH2)15CH* + *` | α-C 加氢 |
| R1268 | `CH3(CH2)15CH* + H* ⇌ CH3(CH2)15CH2* + *` | α-C 加氢 |
| R1269 | `CH3(CH2)15CH2* + H* ⇌ CH3(CH2)15CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1270 | `CH3(CH2)15CH3* ⇌ C17H36(g) + *` | 正构烷烃脱附 |
| R1271 | `CH3(CH2)14CHCH2* ⇌ C17H34(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C17 → C18

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1272 | `CH3(CH2)15C* + C* ⇌ CH3(CH2)15CC* + *` | 与表面 C 偶联 |
| R1273 | `CH3(CH2)15C* + CH* ⇌ CH3(CH2)15CCH* + *` | 次烷基 + CH |
| R1274 | `CH3(CH2)15C* + CH2* ⇌ CH3(CH2)15CCH2* + *` |  |
| R1275 | `CH3(CH2)15C* + CH3* ⇌ CH3(CH2)15CCH3* + *` |  |
| R1276 | `CH3(CH2)15CH* + C* ⇌ CH3(CH2)15CHC* + *` | 与表面 C 偶联 |
| R1277 | `CH3(CH2)15CH* + CH* ⇌ CH3(CH2)15CHCH* + *` | 亚烷基 + CH → 烯基 |
| R1278 | `CH3(CH2)15CH* + CH2* ⇌ CH3(CH2)15CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R1279 | `CH3(CH2)15CH* + CH3* ⇌ CH3(CH2)15CHCH3* + *` |  |
| R1280 | `CH3(CH2)15CH2* + C* ⇌ CH3(CH2)16C* + *` | 与表面 C 偶联 |
| R1281 | `CH3(CH2)15CH2* + CH* ⇌ CH3(CH2)16CH* + *` |  |
| R1282 | `CH3(CH2)15CH2* + CH2* ⇌ CH3(CH2)16CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R1283 | `CH3(CH2)15CH2* + CH3* ⇌ CH3(CH2)16CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1284 | `CH3(CH2)15C* + CO* ⇌ CH3(CH2)15CCO* + *` |  |
| R1285 | `CH3(CH2)15CH* + CO* ⇌ CH3(CH2)15CHCO* + *` |  |
| R1286 | `CH3(CH2)15CH2* + CO* ⇌ CH3(CH2)15CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1287 | `CH3(CH2)15CCO* + H* ⇌ CH3(CH2)15CCHO* + *` | C 端加氢 |
| R1288 | `CH3(CH2)15CCHO* + H* ⇌ CH3(CH2)15CCH2O* + *` | C 端加氢 |
| R1289 | `CH3(CH2)15CCOH* + H* ⇌ CH3(CH2)15CCHOH* + *` | C 端加氢 |
| R1290 | `CH3(CH2)15CCHOH* + H* ⇌ CH3(CH2)15CCH2OH* + *` | C 端加氢 |
| R1291 | `CH3(CH2)15CCO* + H* ⇌ CH3(CH2)15CCOH* + *` | O 端加氢 |
| R1292 | `CH3(CH2)15CCHO* + H* ⇌ CH3(CH2)15CCHOH* + *` | O 端加氢 |
| R1293 | `CH3(CH2)15CCH2O* + H* ⇌ CH3(CH2)15CCH2OH* + *` | O 端加氢 |
| R1294 | `CH3(CH2)15CHCO* + H* ⇌ CH3(CH2)15CHCHO* + *` | C 端加氢 |
| R1295 | `CH3(CH2)15CHCHO* + H* ⇌ CH3(CH2)15CHCH2O* + *` | C 端加氢 |
| R1296 | `CH3(CH2)15CHCOH* + H* ⇌ CH3(CH2)15CHCHOH* + *` | C 端加氢 |
| R1297 | `CH3(CH2)15CHCHOH* + H* ⇌ CH3(CH2)15CHCH2OH* + *` | C 端加氢 |
| R1298 | `CH3(CH2)15CHCO* + H* ⇌ CH3(CH2)15CHCOH* + *` | O 端加氢 |
| R1299 | `CH3(CH2)15CHCHO* + H* ⇌ CH3(CH2)15CHCHOH* + *` | O 端加氢 |
| R1300 | `CH3(CH2)15CHCH2O* + H* ⇌ CH3(CH2)15CHCH2OH* + *` | O 端加氢 |
| R1301 | `CH3(CH2)15CH2CO* + H* ⇌ CH3(CH2)15CH2CHO* + *` | C 端加氢 |
| R1302 | `CH3(CH2)15CH2CHO* + H* ⇌ CH3(CH2)15CH2CH2O* + *` | C 端加氢 |
| R1303 | `CH3(CH2)15CH2COH* + H* ⇌ CH3(CH2)15CH2CHOH* + *` | C 端加氢 |
| R1304 | `CH3(CH2)15CH2CHOH* + H* ⇌ CH3(CH2)15CH2CH2OH* + *` | C 端加氢 |
| R1305 | `CH3(CH2)15CH2CO* + H* ⇌ CH3(CH2)15CH2COH* + *` | O 端加氢 |
| R1306 | `CH3(CH2)15CH2CHO* + H* ⇌ CH3(CH2)15CH2CHOH* + *` | O 端加氢 |
| R1307 | `CH3(CH2)15CH2CH2O* + H* ⇌ CH3(CH2)15CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1308 | `CH3(CH2)15CCO* + * ⇌ CH3(CH2)15CC* + O*` |  |
| R1309 | `CH3(CH2)15CCOH* + * ⇌ CH3(CH2)15CC* + OH*` |  |
| R1310 | `CH3(CH2)15CCHO* + * ⇌ CH3(CH2)15CCH* + O*` |  |
| R1311 | `CH3(CH2)15CCHOH* + * ⇌ CH3(CH2)15CCH* + OH*` |  |
| R1312 | `CH3(CH2)15CCH2O* + * ⇌ CH3(CH2)15CCH2* + O*` |  |
| R1313 | `CH3(CH2)15CCH2OH* + * ⇌ CH3(CH2)15CCH2* + OH*` |  |
| R1314 | `CH3(CH2)15CHCO* + * ⇌ CH3(CH2)15CHC* + O*` |  |
| R1315 | `CH3(CH2)15CHCOH* + * ⇌ CH3(CH2)15CHC* + OH*` |  |
| R1316 | `CH3(CH2)15CHCHO* + * ⇌ CH3(CH2)15CHCH* + O*` |  |
| R1317 | `CH3(CH2)15CHCHOH* + * ⇌ CH3(CH2)15CHCH* + OH*` |  |
| R1318 | `CH3(CH2)15CHCH2O* + * ⇌ CH3(CH2)15CHCH2* + O*` |  |
| R1319 | `CH3(CH2)15CHCH2OH* + * ⇌ CH3(CH2)15CHCH2* + OH*` |  |
| R1320 | `CH3(CH2)15CH2CO* + * ⇌ CH3(CH2)16C* + O*` |  |
| R1321 | `CH3(CH2)15CH2COH* + * ⇌ CH3(CH2)16C* + OH*` |  |
| R1322 | `CH3(CH2)15CH2CHO* + * ⇌ CH3(CH2)16CH* + O*` |  |
| R1323 | `CH3(CH2)15CH2CHOH* + * ⇌ CH3(CH2)16CH* + OH*` |  |
| R1324 | `CH3(CH2)15CH2CH2O* + * ⇌ CH3(CH2)16CH2* + O*` |  |
| R1325 | `CH3(CH2)15CH2CH2OH* + * ⇌ CH3(CH2)16CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1326 | `CH3(CH2)14CH2COH* + HCOH* ⇌ CH3(CH2)15CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1327 | `CH3(CH2)14CHCH* + CH2* ⇌ CH3(CH2)14CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R1328 | `CH3(CH2)14CHCHCH2* ⇌ CH3(CH2)15CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R1329 | `CH3(CH2)14CHCHCH2* + H* ⇌ CH3(CH2)15CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1330 | `CH3(CH2)15CC* + H* ⇌ CH3(CH2)15CCH* + *` | α-C 加氢 |
| R1331 | `CH3(CH2)15CC* + H* ⇌ CH3(CH2)15CHC* + *` | β-C 加氢 |
| R1332 | `CH3(CH2)15CCH* + H* ⇌ CH3(CH2)15CCH2* + *` | α-C 加氢 |
| R1333 | `CH3(CH2)15CCH* + H* ⇌ CH3(CH2)15CHCH* + *` | β-C 加氢 |
| R1334 | `CH3(CH2)15CCH2* + H* ⇌ CH3(CH2)15CCH3* + *` | α-C 加氢 |
| R1335 | `CH3(CH2)15CCH2* + H* ⇌ CH3(CH2)15CHCH2* + *` | β-C 加氢 |
| R1336 | `CH3(CH2)15CCH3* + H* ⇌ CH3(CH2)15CHCH3* + *` | β-C 加氢 |
| R1337 | `CH3(CH2)15CHC* + H* ⇌ CH3(CH2)15CHCH* + *` | α-C 加氢 |
| R1338 | `CH3(CH2)15CHC* + H* ⇌ CH3(CH2)16C* + *` | β-C 加氢 |
| R1339 | `CH3(CH2)15CHCH* + H* ⇌ CH3(CH2)15CHCH2* + *` | α-C 加氢 |
| R1340 | `CH3(CH2)15CHCH* + H* ⇌ CH3(CH2)16CH* + *` | β-C 加氢 |
| R1341 | `CH3(CH2)15CHCH2* + H* ⇌ CH3(CH2)15CHCH3* + *` | α-C 加氢 |
| R1342 | `CH3(CH2)15CHCH2* + H* ⇌ CH3(CH2)16CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R1343 | `CH3(CH2)15CHCH3* + H* ⇌ CH3(CH2)16CH3* + *` | β-C 加氢 |
| R1344 | `CH3(CH2)16C* + H* ⇌ CH3(CH2)16CH* + *` | α-C 加氢 |
| R1345 | `CH3(CH2)16CH* + H* ⇌ CH3(CH2)16CH2* + *` | α-C 加氢 |
| R1346 | `CH3(CH2)16CH2* + H* ⇌ CH3(CH2)16CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1347 | `CH3(CH2)16CH3* ⇌ C18H38(g) + *` | 正构烷烃脱附 |
| R1348 | `CH3(CH2)15CHCH2* ⇌ C18H36(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C18 → C19

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1349 | `CH3(CH2)16C* + C* ⇌ CH3(CH2)16CC* + *` | 与表面 C 偶联 |
| R1350 | `CH3(CH2)16C* + CH* ⇌ CH3(CH2)16CCH* + *` | 次烷基 + CH |
| R1351 | `CH3(CH2)16C* + CH2* ⇌ CH3(CH2)16CCH2* + *` |  |
| R1352 | `CH3(CH2)16C* + CH3* ⇌ CH3(CH2)16CCH3* + *` |  |
| R1353 | `CH3(CH2)16CH* + C* ⇌ CH3(CH2)16CHC* + *` | 与表面 C 偶联 |
| R1354 | `CH3(CH2)16CH* + CH* ⇌ CH3(CH2)16CHCH* + *` | 亚烷基 + CH → 烯基 |
| R1355 | `CH3(CH2)16CH* + CH2* ⇌ CH3(CH2)16CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R1356 | `CH3(CH2)16CH* + CH3* ⇌ CH3(CH2)16CHCH3* + *` |  |
| R1357 | `CH3(CH2)16CH2* + C* ⇌ CH3(CH2)17C* + *` | 与表面 C 偶联 |
| R1358 | `CH3(CH2)16CH2* + CH* ⇌ CH3(CH2)17CH* + *` |  |
| R1359 | `CH3(CH2)16CH2* + CH2* ⇌ CH3(CH2)17CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R1360 | `CH3(CH2)16CH2* + CH3* ⇌ CH3(CH2)17CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1361 | `CH3(CH2)16C* + CO* ⇌ CH3(CH2)16CCO* + *` |  |
| R1362 | `CH3(CH2)16CH* + CO* ⇌ CH3(CH2)16CHCO* + *` |  |
| R1363 | `CH3(CH2)16CH2* + CO* ⇌ CH3(CH2)16CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1364 | `CH3(CH2)16CCO* + H* ⇌ CH3(CH2)16CCHO* + *` | C 端加氢 |
| R1365 | `CH3(CH2)16CCHO* + H* ⇌ CH3(CH2)16CCH2O* + *` | C 端加氢 |
| R1366 | `CH3(CH2)16CCOH* + H* ⇌ CH3(CH2)16CCHOH* + *` | C 端加氢 |
| R1367 | `CH3(CH2)16CCHOH* + H* ⇌ CH3(CH2)16CCH2OH* + *` | C 端加氢 |
| R1368 | `CH3(CH2)16CCO* + H* ⇌ CH3(CH2)16CCOH* + *` | O 端加氢 |
| R1369 | `CH3(CH2)16CCHO* + H* ⇌ CH3(CH2)16CCHOH* + *` | O 端加氢 |
| R1370 | `CH3(CH2)16CCH2O* + H* ⇌ CH3(CH2)16CCH2OH* + *` | O 端加氢 |
| R1371 | `CH3(CH2)16CHCO* + H* ⇌ CH3(CH2)16CHCHO* + *` | C 端加氢 |
| R1372 | `CH3(CH2)16CHCHO* + H* ⇌ CH3(CH2)16CHCH2O* + *` | C 端加氢 |
| R1373 | `CH3(CH2)16CHCOH* + H* ⇌ CH3(CH2)16CHCHOH* + *` | C 端加氢 |
| R1374 | `CH3(CH2)16CHCHOH* + H* ⇌ CH3(CH2)16CHCH2OH* + *` | C 端加氢 |
| R1375 | `CH3(CH2)16CHCO* + H* ⇌ CH3(CH2)16CHCOH* + *` | O 端加氢 |
| R1376 | `CH3(CH2)16CHCHO* + H* ⇌ CH3(CH2)16CHCHOH* + *` | O 端加氢 |
| R1377 | `CH3(CH2)16CHCH2O* + H* ⇌ CH3(CH2)16CHCH2OH* + *` | O 端加氢 |
| R1378 | `CH3(CH2)16CH2CO* + H* ⇌ CH3(CH2)16CH2CHO* + *` | C 端加氢 |
| R1379 | `CH3(CH2)16CH2CHO* + H* ⇌ CH3(CH2)16CH2CH2O* + *` | C 端加氢 |
| R1380 | `CH3(CH2)16CH2COH* + H* ⇌ CH3(CH2)16CH2CHOH* + *` | C 端加氢 |
| R1381 | `CH3(CH2)16CH2CHOH* + H* ⇌ CH3(CH2)16CH2CH2OH* + *` | C 端加氢 |
| R1382 | `CH3(CH2)16CH2CO* + H* ⇌ CH3(CH2)16CH2COH* + *` | O 端加氢 |
| R1383 | `CH3(CH2)16CH2CHO* + H* ⇌ CH3(CH2)16CH2CHOH* + *` | O 端加氢 |
| R1384 | `CH3(CH2)16CH2CH2O* + H* ⇌ CH3(CH2)16CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1385 | `CH3(CH2)16CCO* + * ⇌ CH3(CH2)16CC* + O*` |  |
| R1386 | `CH3(CH2)16CCOH* + * ⇌ CH3(CH2)16CC* + OH*` |  |
| R1387 | `CH3(CH2)16CCHO* + * ⇌ CH3(CH2)16CCH* + O*` |  |
| R1388 | `CH3(CH2)16CCHOH* + * ⇌ CH3(CH2)16CCH* + OH*` |  |
| R1389 | `CH3(CH2)16CCH2O* + * ⇌ CH3(CH2)16CCH2* + O*` |  |
| R1390 | `CH3(CH2)16CCH2OH* + * ⇌ CH3(CH2)16CCH2* + OH*` |  |
| R1391 | `CH3(CH2)16CHCO* + * ⇌ CH3(CH2)16CHC* + O*` |  |
| R1392 | `CH3(CH2)16CHCOH* + * ⇌ CH3(CH2)16CHC* + OH*` |  |
| R1393 | `CH3(CH2)16CHCHO* + * ⇌ CH3(CH2)16CHCH* + O*` |  |
| R1394 | `CH3(CH2)16CHCHOH* + * ⇌ CH3(CH2)16CHCH* + OH*` |  |
| R1395 | `CH3(CH2)16CHCH2O* + * ⇌ CH3(CH2)16CHCH2* + O*` |  |
| R1396 | `CH3(CH2)16CHCH2OH* + * ⇌ CH3(CH2)16CHCH2* + OH*` |  |
| R1397 | `CH3(CH2)16CH2CO* + * ⇌ CH3(CH2)17C* + O*` |  |
| R1398 | `CH3(CH2)16CH2COH* + * ⇌ CH3(CH2)17C* + OH*` |  |
| R1399 | `CH3(CH2)16CH2CHO* + * ⇌ CH3(CH2)17CH* + O*` |  |
| R1400 | `CH3(CH2)16CH2CHOH* + * ⇌ CH3(CH2)17CH* + OH*` |  |
| R1401 | `CH3(CH2)16CH2CH2O* + * ⇌ CH3(CH2)17CH2* + O*` |  |
| R1402 | `CH3(CH2)16CH2CH2OH* + * ⇌ CH3(CH2)17CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1403 | `CH3(CH2)15CH2COH* + HCOH* ⇌ CH3(CH2)16CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1404 | `CH3(CH2)15CHCH* + CH2* ⇌ CH3(CH2)15CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R1405 | `CH3(CH2)15CHCHCH2* ⇌ CH3(CH2)16CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R1406 | `CH3(CH2)15CHCHCH2* + H* ⇌ CH3(CH2)16CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1407 | `CH3(CH2)16CC* + H* ⇌ CH3(CH2)16CCH* + *` | α-C 加氢 |
| R1408 | `CH3(CH2)16CC* + H* ⇌ CH3(CH2)16CHC* + *` | β-C 加氢 |
| R1409 | `CH3(CH2)16CCH* + H* ⇌ CH3(CH2)16CCH2* + *` | α-C 加氢 |
| R1410 | `CH3(CH2)16CCH* + H* ⇌ CH3(CH2)16CHCH* + *` | β-C 加氢 |
| R1411 | `CH3(CH2)16CCH2* + H* ⇌ CH3(CH2)16CCH3* + *` | α-C 加氢 |
| R1412 | `CH3(CH2)16CCH2* + H* ⇌ CH3(CH2)16CHCH2* + *` | β-C 加氢 |
| R1413 | `CH3(CH2)16CCH3* + H* ⇌ CH3(CH2)16CHCH3* + *` | β-C 加氢 |
| R1414 | `CH3(CH2)16CHC* + H* ⇌ CH3(CH2)16CHCH* + *` | α-C 加氢 |
| R1415 | `CH3(CH2)16CHC* + H* ⇌ CH3(CH2)17C* + *` | β-C 加氢 |
| R1416 | `CH3(CH2)16CHCH* + H* ⇌ CH3(CH2)16CHCH2* + *` | α-C 加氢 |
| R1417 | `CH3(CH2)16CHCH* + H* ⇌ CH3(CH2)17CH* + *` | β-C 加氢 |
| R1418 | `CH3(CH2)16CHCH2* + H* ⇌ CH3(CH2)16CHCH3* + *` | α-C 加氢 |
| R1419 | `CH3(CH2)16CHCH2* + H* ⇌ CH3(CH2)17CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R1420 | `CH3(CH2)16CHCH3* + H* ⇌ CH3(CH2)17CH3* + *` | β-C 加氢 |
| R1421 | `CH3(CH2)17C* + H* ⇌ CH3(CH2)17CH* + *` | α-C 加氢 |
| R1422 | `CH3(CH2)17CH* + H* ⇌ CH3(CH2)17CH2* + *` | α-C 加氢 |
| R1423 | `CH3(CH2)17CH2* + H* ⇌ CH3(CH2)17CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1424 | `CH3(CH2)17CH3* ⇌ C19H40(g) + *` | 正构烷烃脱附 |
| R1425 | `CH3(CH2)16CHCH2* ⇌ C19H38(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |

## C19 → C20

### 5a 链增长: C–C 偶联(碳化物机理)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1426 | `CH3(CH2)17C* + C* ⇌ CH3(CH2)17CC* + *` | 与表面 C 偶联 |
| R1427 | `CH3(CH2)17C* + CH* ⇌ CH3(CH2)17CCH* + *` | 次烷基 + CH |
| R1428 | `CH3(CH2)17C* + CH2* ⇌ CH3(CH2)17CCH2* + *` |  |
| R1429 | `CH3(CH2)17C* + CH3* ⇌ CH3(CH2)17CCH3* + *` |  |
| R1430 | `CH3(CH2)17CH* + C* ⇌ CH3(CH2)17CHC* + *` | 与表面 C 偶联 |
| R1431 | `CH3(CH2)17CH* + CH* ⇌ CH3(CH2)17CHCH* + *` | 亚烷基 + CH → 烯基 |
| R1432 | `CH3(CH2)17CH* + CH2* ⇌ CH3(CH2)17CHCH2* + *` | 亚烷基 + CH2 → 吸附 α-烯烃 |
| R1433 | `CH3(CH2)17CH* + CH3* ⇌ CH3(CH2)17CHCH3* + *` |  |
| R1434 | `CH3(CH2)17CH2* + C* ⇌ CH3(CH2)18C* + *` | 与表面 C 偶联 |
| R1435 | `CH3(CH2)17CH2* + CH* ⇌ CH3(CH2)18CH* + *` |  |
| R1436 | `CH3(CH2)17CH2* + CH2* ⇌ CH3(CH2)18CH2* + *` | 烷基机理(Brady–Pettit): CH2 插入金属–烷基键 |
| R1437 | `CH3(CH2)17CH2* + CH3* ⇌ CH3(CH2)18CH3* + *` |  |

### 5b 链增长: CO 插入

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1438 | `CH3(CH2)17C* + CO* ⇌ CH3(CH2)17CCO* + *` |  |
| R1439 | `CH3(CH2)17CH* + CO* ⇌ CH3(CH2)17CHCO* + *` |  |
| R1440 | `CH3(CH2)17CH2* + CO* ⇌ CH3(CH2)17CH2CO* + *` | CO 插入金属–烷基键(Pichler–Schulz) |

### 5c CO 插入中间体加氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1441 | `CH3(CH2)17CCO* + H* ⇌ CH3(CH2)17CCHO* + *` | C 端加氢 |
| R1442 | `CH3(CH2)17CCHO* + H* ⇌ CH3(CH2)17CCH2O* + *` | C 端加氢 |
| R1443 | `CH3(CH2)17CCOH* + H* ⇌ CH3(CH2)17CCHOH* + *` | C 端加氢 |
| R1444 | `CH3(CH2)17CCHOH* + H* ⇌ CH3(CH2)17CCH2OH* + *` | C 端加氢 |
| R1445 | `CH3(CH2)17CCO* + H* ⇌ CH3(CH2)17CCOH* + *` | O 端加氢 |
| R1446 | `CH3(CH2)17CCHO* + H* ⇌ CH3(CH2)17CCHOH* + *` | O 端加氢 |
| R1447 | `CH3(CH2)17CCH2O* + H* ⇌ CH3(CH2)17CCH2OH* + *` | O 端加氢 |
| R1448 | `CH3(CH2)17CHCO* + H* ⇌ CH3(CH2)17CHCHO* + *` | C 端加氢 |
| R1449 | `CH3(CH2)17CHCHO* + H* ⇌ CH3(CH2)17CHCH2O* + *` | C 端加氢 |
| R1450 | `CH3(CH2)17CHCOH* + H* ⇌ CH3(CH2)17CHCHOH* + *` | C 端加氢 |
| R1451 | `CH3(CH2)17CHCHOH* + H* ⇌ CH3(CH2)17CHCH2OH* + *` | C 端加氢 |
| R1452 | `CH3(CH2)17CHCO* + H* ⇌ CH3(CH2)17CHCOH* + *` | O 端加氢 |
| R1453 | `CH3(CH2)17CHCHO* + H* ⇌ CH3(CH2)17CHCHOH* + *` | O 端加氢 |
| R1454 | `CH3(CH2)17CHCH2O* + H* ⇌ CH3(CH2)17CHCH2OH* + *` | O 端加氢 |
| R1455 | `CH3(CH2)17CH2CO* + H* ⇌ CH3(CH2)17CH2CHO* + *` | C 端加氢 |
| R1456 | `CH3(CH2)17CH2CHO* + H* ⇌ CH3(CH2)17CH2CH2O* + *` | C 端加氢 |
| R1457 | `CH3(CH2)17CH2COH* + H* ⇌ CH3(CH2)17CH2CHOH* + *` | C 端加氢 |
| R1458 | `CH3(CH2)17CH2CHOH* + H* ⇌ CH3(CH2)17CH2CH2OH* + *` | C 端加氢 |
| R1459 | `CH3(CH2)17CH2CO* + H* ⇌ CH3(CH2)17CH2COH* + *` | O 端加氢 |
| R1460 | `CH3(CH2)17CH2CHO* + H* ⇌ CH3(CH2)17CH2CHOH* + *` | O 端加氢 |
| R1461 | `CH3(CH2)17CH2CH2O* + H* ⇌ CH3(CH2)17CH2CH2OH* + *` | O 端加氢 |

### 5d CO 插入中间体 C–O 断裂

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1462 | `CH3(CH2)17CCO* + * ⇌ CH3(CH2)17CC* + O*` |  |
| R1463 | `CH3(CH2)17CCOH* + * ⇌ CH3(CH2)17CC* + OH*` |  |
| R1464 | `CH3(CH2)17CCHO* + * ⇌ CH3(CH2)17CCH* + O*` |  |
| R1465 | `CH3(CH2)17CCHOH* + * ⇌ CH3(CH2)17CCH* + OH*` |  |
| R1466 | `CH3(CH2)17CCH2O* + * ⇌ CH3(CH2)17CCH2* + O*` |  |
| R1467 | `CH3(CH2)17CCH2OH* + * ⇌ CH3(CH2)17CCH2* + OH*` |  |
| R1468 | `CH3(CH2)17CHCO* + * ⇌ CH3(CH2)17CHC* + O*` |  |
| R1469 | `CH3(CH2)17CHCOH* + * ⇌ CH3(CH2)17CHC* + OH*` |  |
| R1470 | `CH3(CH2)17CHCHO* + * ⇌ CH3(CH2)17CHCH* + O*` |  |
| R1471 | `CH3(CH2)17CHCHOH* + * ⇌ CH3(CH2)17CHCH* + OH*` |  |
| R1472 | `CH3(CH2)17CHCH2O* + * ⇌ CH3(CH2)17CHCH2* + O*` |  |
| R1473 | `CH3(CH2)17CHCH2OH* + * ⇌ CH3(CH2)17CHCH2* + OH*` |  |
| R1474 | `CH3(CH2)17CH2CO* + * ⇌ CH3(CH2)18C* + O*` |  |
| R1475 | `CH3(CH2)17CH2COH* + * ⇌ CH3(CH2)18C* + OH*` |  |
| R1476 | `CH3(CH2)17CH2CHO* + * ⇌ CH3(CH2)18CH* + O*` |  |
| R1477 | `CH3(CH2)17CH2CHOH* + * ⇌ CH3(CH2)18CH* + OH*` |  |
| R1478 | `CH3(CH2)17CH2CH2O* + * ⇌ CH3(CH2)18CH2* + O*` |  |
| R1479 | `CH3(CH2)17CH2CH2OH* + * ⇌ CH3(CH2)18CH2* + OH*` |  |

### 5e 链增长: 烯醇缩合(Storch–Anderson)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1480 | `CH3(CH2)16CH2COH* + HCOH* ⇌ CH3(CH2)17CCOH* + H2O*` | RC(OH)* + HCOH* 脱水缩合 |

### 5f 链增长: 烯基机理(Maitlis)

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1481 | `CH3(CH2)16CHCH* + CH2* ⇌ CH3(CH2)16CHCHCH2* + *` | 烯基 + CH2 → 烯丙基 |
| R1482 | `CH3(CH2)16CHCHCH2* ⇌ CH3(CH2)17CHCH*` | 烯丙基 ⇌ 烯基 (1,3-H 迁移) |
| R1483 | `CH3(CH2)16CHCHCH2* + H* ⇌ CH3(CH2)17CHCH2* + *` | 烯丙基加氢 → 吸附 α-烯烃 |

### 5g C_n 中间体加氢/脱氢

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1484 | `CH3(CH2)17CC* + H* ⇌ CH3(CH2)17CCH* + *` | α-C 加氢 |
| R1485 | `CH3(CH2)17CC* + H* ⇌ CH3(CH2)17CHC* + *` | β-C 加氢 |
| R1486 | `CH3(CH2)17CCH* + H* ⇌ CH3(CH2)17CCH2* + *` | α-C 加氢 |
| R1487 | `CH3(CH2)17CCH* + H* ⇌ CH3(CH2)17CHCH* + *` | β-C 加氢 |
| R1488 | `CH3(CH2)17CCH2* + H* ⇌ CH3(CH2)17CCH3* + *` | α-C 加氢 |
| R1489 | `CH3(CH2)17CCH2* + H* ⇌ CH3(CH2)17CHCH2* + *` | β-C 加氢 |
| R1490 | `CH3(CH2)17CCH3* + H* ⇌ CH3(CH2)17CHCH3* + *` | β-C 加氢 |
| R1491 | `CH3(CH2)17CHC* + H* ⇌ CH3(CH2)17CHCH* + *` | α-C 加氢 |
| R1492 | `CH3(CH2)17CHC* + H* ⇌ CH3(CH2)18C* + *` | β-C 加氢 |
| R1493 | `CH3(CH2)17CHCH* + H* ⇌ CH3(CH2)17CHCH2* + *` | α-C 加氢 |
| R1494 | `CH3(CH2)17CHCH* + H* ⇌ CH3(CH2)18CH* + *` | β-C 加氢 |
| R1495 | `CH3(CH2)17CHCH2* + H* ⇌ CH3(CH2)17CHCH3* + *` | α-C 加氢 |
| R1496 | `CH3(CH2)17CHCH2* + H* ⇌ CH3(CH2)18CH2* + *` | 逆反应 = 烷基 β-H 消除 → α-烯烃(链终止) |
| R1497 | `CH3(CH2)17CHCH3* + H* ⇌ CH3(CH2)18CH3* + *` | β-C 加氢 |
| R1498 | `CH3(CH2)18C* + H* ⇌ CH3(CH2)18CH* + *` | α-C 加氢 |
| R1499 | `CH3(CH2)18CH* + H* ⇌ CH3(CH2)18CH2* + *` | α-C 加氢 |
| R1500 | `CH3(CH2)18CH2* + H* ⇌ CH3(CH2)18CH3* + *` | 链终止: 烷基加氢(还原消除) → 正构烷烃 |

### 5h 产物脱附/再吸附

| 编号 | 基元步骤 | 说明 |
|---|---|---|
| R1501 | `CH3(CH2)18CH3* ⇌ C20H42(g) + *` | 正构烷烃脱附 |
| R1502 | `CH3(CH2)17CHCH2* ⇌ C20H40(g) + *` | α-烯烃脱附; 逆向 = 再吸附(二次加氢/再引发) |
