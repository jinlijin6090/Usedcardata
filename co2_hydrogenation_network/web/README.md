# 交互网页：CO₂ 一步热催化加氢 C1–C16 反应网络

这里有两个页面。

| 页面 | 网络 | 在线链接（私有，需在页面 Share 菜单共享后别人才能打开） |
|---|---|---|
| `index.html`（主页面） | 推荐版：C1–C12 全部构造异构体，C13–C16 按 ZSM-5 择形。128 384 个物种、367 369 条基元反应 | https://claude.ai/artifact/HZveP1trkXVhuSSLBMFgQs |
| `lod.html`（分级加载） | 全量版：C1–C16 全部构造异构体（规模见 `../data/net_C16_full/stats.json`） | 见下文“全量版” |

## 主页面功能

- **布局**：列为碳数，行为物种类别，背景色表示所需的催化剂形式。C13–C16 区域加了浅色底，表示经过择形筛选。
- **浏览**：WebGL2 一次画出全部点和连线，支持拖动、滚轮/双指缩放。放大后显示每个点的 SMILES。
- **全部反应**：单击或搜索一个物种后，右侧按反应类型列出它参与的**每一条基元反应方程式**，方程式里的物种可以点击跳转。
- **从 CO₂ + H₂ 的路径**：
  - 把每条可逆反应拆成正、逆两条有向超边，从 CO₂ 出发做超图分层 BFS，得到每个物种的最少步数。H₂、H*、H⁺、H₂O、[H]、H• 视为始终可用。
  - 只把“多甲基苯 → 烃池”当作单向反应。
  - “最短路径”的计法：每一步只跟随步数恰好少 1 的含碳主反应物。页面给出最短路径的总数，列出前 40 条，并在画布上用橙色编号画出所选路径。
  - 含可逆反应的网络里“全部路径”有无穷多条，所以只枚举最短路径。
- **筛选**：可按物种类别、反应类型、碳数范围筛选，也可按名称或 SMILES 搜索。

## 离线版

`offline/CO2加氢反应网络_离线版.zip`：解压后双击 HTML 就能打开，数据和结构式绘制库都内置在文件里，不需要联网。

重建命令：`python build_offline.py`。它会把 `nodes.json`、`reactions.json` 和 `vendor/smiles-drawer.min.js` 内联进 HTML。

## 分级加载页面（全量 C16）

全量网络有几百万个物种，不能一次下载。`lod.html` 按缩放程度分四级显示，右上角会标出当前是哪一级：

| 级别 | 何时出现 | 需要下载的数据 |
|---|---|---|
| ① 格子 + 格间反应流量 | 打开页面即显示 | 只有 `meta.json`：各格的物种数，以及服务器端汇总的格间反应条数 |
| ② 单个物种的点 | 放大后 | 不需要下载：位置由物种编号和格子网格直接推算 |
| ③ 标签 | 放大到能分辨单个物种 | 下载这一区域的 `l/<k>.bin`：类别、到 CO₂ 的最少步数、最短路径数（log10）、SMILES |
| ④ 每条反应 | 再放大，或选中某个物种 | 下载 `x/<k>.bin`：该区域每个物种参与的全部反应，记录反应类型、自身在哪一侧，以及其他参与者编号（编码为相对自身编号的差值） |

- 数据块按物种编号切分，每块约 4 万个物种，用 gzip 压缩，浏览器里用 `DecompressionStream` 解压。
- 路径搜索只沿路径逐步下载需要的数据块。最少步数和最短路径数已经由 `lod_build.py` 在服务器端算好。

重建命令：
```
python lod_build.py ../data/net_C16_full lod_full 40000
```
然后把 `lod.html` 放进 `lod_full/` 目录（命名为 index.html）发布，或者在该目录运行 `python -m http.server` 后本地打开。

## 文件

| 文件 | 说明 |
|---|---|
| `index.html` / `nodes.json` / `reactions.json` | 主页面及其数据（物种编号与 `../data/net_C16_shape/species.tsv.gz` 的 id 一致） |
| `web_export.py` | 从 `network_gen.py` 的输出生成 `nodes.json` 和 `reactions.json` |
| `build_offline.py`、`offline/` | 单文件离线版 |
| `lod.html`、`lod_build.py` | 分级加载页面和数据块生成脚本 |
| `vendor/` | smiles-drawer 2.4.1（MIT 许可证） |

重建主页面数据：
```
python ../scripts/network_gen.py --cmax 16 --full-upto 12 --out ../data/net_C16_shape
python web_export.py ../data/net_C16_shape .
```
