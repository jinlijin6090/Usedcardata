# 反应网络 PPT

`CO2加氢制航煤_反应网络.pptx`：把当前 OX-ZEO 集总反应网络画成网络图，共 8 页。网络图都用原生 PowerPoint 形状绘制，可以直接编辑；第 5 页裂解/齐聚弧线图是图片。通量来自 `../tandem_kinetics/calibration.json` 中的 G1 校准参数。

重新生成（需要 node 和 pptxgenjs，以及 pptx skill 的 apply_theme.js）：

```bash
python fluxes.py && python prep.py && python arcs.py
NODE_PATH=<node_modules> node build.js deck.pptx
```

## 完整反应网络（含自由基）

`CO2加氢C1催化完整反应网络_含自由基.pptx`：参考"零碳 / 低碳催化反应网络"总图的画法，去掉了氨 / 尿素（C–N 偶联）部分，共 3 页，幻灯片尺寸 20 × 12.5 英寸。

1. 总图，65 个节点、75 条反应边，分为以下几个区：
   * 电化学还原；
   * C₁ 表面中间体核心：甲酸盐路线、羧基 / RWGS、CO 加氢与解离、CHₓ 链；
   * 与表面物种一一对应的气相自由基区：HOCO•、HCOO•、HCO•、CH₃O•、•CH₂OH，以及 POM / OCM 中的 •CH₃、CH₃OO•、H•、HO₂•、OH•；
   * 亚硝酸甲酯羰基化制乙二醇；
   * 分子筛（OX-ZEO）：表面甲氧基、乙烯酮、烃池、烯烃循环；
   * 金属位链增长（费托 / FTO、CO 插入制醇）与酸位升级到航煤 / 柴油。

   节点形状表示物种类型，连线颜色表示反应类别，节点右下的圆点表示经过该中间体的反应类别。
2. 全部自由基清单：每个自由基的生成反应、消耗反应和对应的表面物种。
3. 共享中间体表，以及节点和反应边的统计。

```bash
NODE_PATH=<node_modules> node full_network.js out.pptx
```
