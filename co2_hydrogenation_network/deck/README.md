# 反应网络 PPT

`CO2加氢制航煤_反应网络.pptx`：把当前 OX-ZEO 集总反应网络画成网络图，共 8 页。网络图都用原生 PowerPoint 形状绘制，可以直接编辑；第 5 页裂解/齐聚弧线图是图片。通量来自 `../tandem_kinetics/calibration.json` 中的 G1 校准参数。

重新生成（需要 node 和 pptxgenjs，以及 pptx skill 的 apply_theme.js）：

```bash
python fluxes.py && python prep.py && python arcs.py
NODE_PATH=<node_modules> node build.js deck.pptx
```
