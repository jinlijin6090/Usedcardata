# 汽车之家二手车数据爬虫

基于 Python 的汽车之家（che168.com）二手车数据爬虫，支持按城市、关键词筛选，结果导出为 CSV 和 JSON 格式。

## 功能特性

- 支持按**城市**和**关键词**（车型/品牌）筛选
- 支持多页爬取，可自定义最大页数
- 自动随机延时，降低被封 IP 风险
- 支持 HTTP 代理
- 数据同时导出为 **CSV** 和 **JSON**
- 字段包含：车型名称、价格、年份、里程、城市、变速箱、经销商、详情页链接等

## 环境要求

- Python 3.10+

## 安装

```bash
pip install -r requirements.txt
```

## 使用方法

### 基本用法

```bash
# 爬取全国二手车（默认 5 页）
python scraper.py

# 按城市 + 关键词爬取
python scraper.py --city 北京 --keyword 大众

# 指定更多页数
python scraper.py --city 上海 --keyword 宝马3系 --pages 10

# 自定义输出路径
python scraper.py --city 成都 --output data/chengdu_cars
```

### 完整参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--city` | `全国` | 目标城市 |
| `--keyword` | _(空)_ | 搜索关键词，如"大众"、"SUV" |
| `--pages` | `5` | 最多爬取页数 |
| `--output` | `output/used_cars` | 输出路径（不含扩展名） |
| `--proxy` | _(空)_ | HTTP 代理，如 `http://127.0.0.1:7890` |
| `--timeout` | `15` | 单次请求超时（秒） |

### 支持城市

| 中文 | 编码 |
|------|------|
| 全国 | china |
| 北京 | beijing |
| 上海 | shanghai |
| 广州 | guangzhou |
| 深圳 | shenzhen |
| 成都 | chengdu |
| 杭州 | hangzhou |
| 武汉 | wuhan |
| 南京 | nanjing |
| 重庆 | chongqing |
| 西安 | xian |
| 苏州 | suzhou |
| 天津 | tianjin |
| 郑州 | zhengzhou |
| 长沙 | changsha |

## 输出格式

运行后在 `output/` 目录生成两个文件：

- `used_cars.csv` — Excel 可直接打开（UTF-8 BOM 编码）
- `used_cars.json` — 结构化 JSON 数据

### CSV 字段说明

| 字段 | 说明 |
|------|------|
| `title` | 车型名称 |
| `price` | 价格（万元） |
| `year` | 年份 |
| `mileage` | 里程（万公里） |
| `city` | 所在城市 |
| `transmission` | 变速箱类型 |
| `displacement` | 排量 |
| `color` | 颜色 |
| `dealer` | 经销商/卖家 |
| `url` | 详情页链接 |
| `img_url` | 封面图片链接 |

## 注意事项

1. **仅供学习研究**，请遵守目标网站的 robots.txt 及相关法律法规。
2. 网站结构可能变化，如爬取数据为空，需根据实际 HTML 更新 CSS 选择器（见 `_extract_car` 和 `_parse_listings` 方法）。
3. 建议使用代理 IP 轮换（`--proxy`）进行大量数据采集，避免 IP 被封。
4. 默认延时 1.5～3.5 秒/页，可在 `_random_sleep` 调用处调整。

## 项目结构

```
Usedcardata/
├── scraper.py        # 主爬虫程序
├── requirements.txt  # Python 依赖
├── README.md         # 本文档
└── output/           # 爬取结果（自动创建）
    ├── used_cars.csv
    └── used_cars.json
```
