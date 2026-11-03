# 交互网页：C1–C12 全同分异构体反应网络

在线版本：https://claude.ai/artifact/HZveP1trkXVhuSSLBMFgQs （私有链接，需要在页面的 Share 菜单里共享后，其他人才能打开）

- 共 88 139 个物种、286 827 条基元反应，由 `../scripts/isomer_scale.py` 用 RDKit 显式枚举。
- 布局与单页 pptx 相同：列 = 碳数，行 = 物种类别，背景色 = 催化剂形式。
- 用 WebGL2 一次画出全部点和边。支持拖动、滚轮或双指缩放；放大后显示每个点的 SMILES；悬停查看信息，单击查看结构式和按反应类型分组的全部相邻物种。
- 可按名称或 SMILES 搜索，可按物种类别、反应类型、碳数范围筛选。

| 文件 | 说明 |
|---|---|
| `index.html` | 页面；结构式用 smiles-drawer 2.4.1（jsDelivr） |
| `network.json` | 布局和数据（base64 类型数组，7.1 MB） |
| `web_export.py` | 由图数据生成 `network.json` |

重建步骤：
```
DUMP=../data/isomer_graph_C12.json.gz python ../scripts/isomer_scale.py 12
python web_export.py ../data/isomer_graph_C12.json.gz network.json
```
本地预览：在本目录运行 `python -m http.server`，然后打开 index.html。页面没有 `<html>/<head>` 骨架，发布时由 Artifact 平台补上，浏览器也能直接渲染。
