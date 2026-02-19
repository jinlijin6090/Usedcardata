#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
汽车之家二手车数据爬虫
目标网站: https://www.che168.com (汽车之家旗下二手车平台)
"""

import csv
import json
import time
import random
import logging
import argparse
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import requests
from bs4 import BeautifulSoup

# ─────────────────────────── 日志配置 ───────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ─────────────────────────── 数据结构 ───────────────────────────
@dataclass
class CarListing:
    """二手车信息"""
    title: str = ""           # 车型名称
    price: str = ""           # 价格（万元）
    year: str = ""            # 年份
    mileage: str = ""         # 里程（万公里）
    city: str = ""            # 所在城市
    transmission: str = ""    # 变速箱类型
    displacement: str = ""    # 排量
    color: str = ""           # 颜色
    dealer: str = ""          # 经销商/卖家
    url: str = ""             # 详情页链接
    img_url: str = ""         # 封面图片链接


# ─────────────────────────── 爬虫配置 ───────────────────────────
# 常用请求头，模拟浏览器访问
DEFAULT_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/122.0.0.0 Safari/537.36"
    ),
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "zh-CN,zh;q=0.9,en;q=0.8",
    "Accept-Encoding": "gzip, deflate, br",
    "Connection": "keep-alive",
    "Referer": "https://www.che168.com/",
}

# 省会城市 -> che168 城市编码（可按需扩充）
CITY_CODES = {
    "全国": "china",
    "北京": "beijing",
    "上海": "shanghai",
    "广州": "guangzhou",
    "深圳": "shenzhen",
    "成都": "chengdu",
    "杭州": "hangzhou",
    "武汉": "wuhan",
    "南京": "nanjing",
    "重庆": "chongqing",
    "西安": "xian",
    "苏州": "suzhou",
    "天津": "tianjin",
    "郑州": "zhengzhou",
    "长沙": "changsha",
}


# ─────────────────────────── 工具函数 ───────────────────────────
def _build_url(city_code: str, page: int, keyword: str = "") -> str:
    """
    构造 che168 列表页 URL。

    che168 列表页示例:
      https://www.che168.com/beijing/a0_0msdgscncgpi1ltocsp{page}exi0/
    关键词搜索:
      https://www.che168.com/beijing/a0_0msdgscncgpi1ltocsp1exi0/?keyWord=大众
    """
    base = f"https://www.che168.com/{city_code}/a0_0msdgscncgpi1ltocsp{page}exi0/"
    if keyword:
        base += f"?keyWord={keyword}"
    return base


def _safe_text(tag, default: str = "") -> str:
    """安全提取标签文本"""
    return tag.get_text(strip=True) if tag else default


def _random_sleep(min_s: float = 1.0, max_s: float = 3.0) -> None:
    """随机等待，降低被封风险"""
    time.sleep(random.uniform(min_s, max_s))


# ─────────────────────────── 核心爬虫 ───────────────────────────
class Che168Scraper:
    """汽车之家(che168)二手车列表爬虫"""

    def __init__(
        self,
        city: str = "全国",
        keyword: str = "",
        max_pages: int = 5,
        request_timeout: int = 15,
        proxy: Optional[str] = None,
    ):
        """
        参数:
            city        : 目标城市（见 CITY_CODES）或自定义城市编码
            keyword     : 搜索关键词，如"大众"、"宝马3系"（可为空，爬全部）
            max_pages   : 最多爬取页数
            request_timeout: 请求超时（秒）
            proxy       : HTTP 代理，如 "http://127.0.0.1:7890"
        """
        self.city_code = CITY_CODES.get(city, city)
        self.keyword = keyword
        self.max_pages = max_pages
        self.timeout = request_timeout

        self.session = requests.Session()
        self.session.headers.update(DEFAULT_HEADERS)
        if proxy:
            self.session.proxies = {"http": proxy, "https": proxy}

    # ── 请求单页 ──────────────────────────────────────────────
    def _fetch_page(self, url: str) -> Optional[BeautifulSoup]:
        """获取页面并解析为 BeautifulSoup 对象"""
        try:
            resp = self.session.get(url, timeout=self.timeout)
            resp.raise_for_status()
            resp.encoding = "gbk"          # che168 常用 GBK 编码
            return BeautifulSoup(resp.text, "html.parser")
        except requests.exceptions.RequestException as exc:
            logger.error("请求失败 [%s]: %s", url, exc)
            return None

    # ── 解析列表页 ────────────────────────────────────────────
    def _parse_listings(self, soup: BeautifulSoup) -> list[CarListing]:
        """从列表页 soup 中提取所有车辆信息"""
        listings: list[CarListing] = []

        # che168 列表项容器：<ul class="ucar-list"> <li> ...
        items = soup.select("ul.ucar-list > li") or soup.select("ul.list-ul > li")

        if not items:
            logger.warning("未找到车辆列表，页面结构可能已更新")
            return listings

        for item in items:
            try:
                car = self._extract_car(item)
                if car.title:
                    listings.append(car)
            except Exception as exc:  # noqa: BLE001
                logger.debug("解析单条数据出错: %s", exc)

        return listings

    def _extract_car(self, item) -> CarListing:
        """从单个列表 <li> 标签提取车辆字段"""
        car = CarListing()

        # --- 标题 ---
        title_tag = item.select_one(".card-name") or item.select_one(".tit-name")
        car.title = _safe_text(title_tag)

        # --- 价格 ---
        price_tag = (
            item.select_one(".pirce em")
            or item.select_one(".price-num")
            or item.select_one(".sell-price")
        )
        car.price = _safe_text(price_tag)

        # --- 年份 / 里程 / 城市（通常在同一行） ---
        info_tags = item.select(".param-list li") or item.select(".card-attr li")
        info_texts = [_safe_text(t) for t in info_tags]

        if len(info_texts) >= 1:
            car.year = info_texts[0]
        if len(info_texts) >= 2:
            car.mileage = info_texts[1]
        if len(info_texts) >= 3:
            car.city = info_texts[2]
        if len(info_texts) >= 4:
            car.transmission = info_texts[3]

        # --- 经销商 / 卖家 ---
        dealer_tag = item.select_one(".person-name") or item.select_one(".shop-name")
        car.dealer = _safe_text(dealer_tag)

        # --- 详情页链接 ---
        link_tag = item.select_one("a[href]")
        if link_tag:
            href = link_tag["href"]
            car.url = href if href.startswith("http") else f"https://www.che168.com{href}"

        # --- 封面图 ---
        img_tag = item.select_one("img[src]") or item.select_one("img[data-src]")
        if img_tag:
            car.img_url = img_tag.get("data-src") or img_tag.get("src", "")

        return car

    # ── 判断是否还有下一页 ────────────────────────────────────
    def _has_next_page(self, soup: BeautifulSoup, current_page: int) -> bool:
        """检查是否存在下一页"""
        next_btn = soup.select_one(".page-next:not(.disabled)") or soup.select_one(
            "a.nextpage"
        )
        if next_btn:
            return True
        # 备用：检查页码按钮
        page_links = soup.select(".page-item a, .pages a")
        page_numbers = []
        for a in page_links:
            try:
                page_numbers.append(int(a.get_text(strip=True)))
            except ValueError:
                pass
        return bool(page_numbers) and max(page_numbers, default=0) > current_page

    # ── 主入口 ────────────────────────────────────────────────
    def run(self) -> list[CarListing]:
        """执行爬取，返回所有车辆信息列表"""
        all_cars: list[CarListing] = []

        for page in range(1, self.max_pages + 1):
            url = _build_url(self.city_code, page, self.keyword)
            logger.info("正在爬取第 %d 页: %s", page, url)

            soup = self._fetch_page(url)
            if soup is None:
                logger.warning("第 %d 页请求失败，停止爬取", page)
                break

            cars = self._parse_listings(soup)
            logger.info("第 %d 页解析到 %d 条数据", page, len(cars))
            all_cars.extend(cars)

            if not cars or not self._has_next_page(soup, page):
                logger.info("已到最后一页，共爬取 %d 条数据", len(all_cars))
                break

            _random_sleep(1.5, 3.5)   # 礼貌爬取，避免封 IP

        return all_cars


# ─────────────────────────── 数据导出 ───────────────────────────
def save_to_csv(cars: list[CarListing], path: str) -> None:
    """保存为 CSV 文件"""
    if not cars:
        logger.warning("无数据可保存")
        return

    fieldnames = list(asdict(cars[0]).keys())
    Path(path).parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows([asdict(c) for c in cars])

    logger.info("已保存 %d 条数据到 CSV: %s", len(cars), path)


def save_to_json(cars: list[CarListing], path: str) -> None:
    """保存为 JSON 文件"""
    if not cars:
        logger.warning("无数据可保存")
        return

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump([asdict(c) for c in cars], f, ensure_ascii=False, indent=2)

    logger.info("已保存 %d 条数据到 JSON: %s", len(cars), path)


# ─────────────────────────── CLI 入口 ───────────────────────────
def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="汽车之家(che168)二手车数据爬虫",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--city", default="全国",
        help=f"目标城市，可选: {', '.join(CITY_CODES.keys())}",
    )
    parser.add_argument(
        "--keyword", default="",
        help="搜索关键词，如：大众、宝马3系、SUV（留空则爬取全部）",
    )
    parser.add_argument(
        "--pages", type=int, default=5,
        help="最多爬取页数",
    )
    parser.add_argument(
        "--output", default="output/used_cars",
        help="输出文件路径（不含扩展名），同时生成 .csv 和 .json",
    )
    parser.add_argument(
        "--proxy", default="",
        help="HTTP 代理地址，如 http://127.0.0.1:7890",
    )
    parser.add_argument(
        "--timeout", type=int, default=15,
        help="单次请求超时（秒）",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    logger.info("═" * 50)
    logger.info("汽车之家二手车爬虫启动")
    logger.info("城市: %s | 关键词: %s | 最大页数: %d", args.city, args.keyword or "全部", args.pages)
    logger.info("═" * 50)

    scraper = Che168Scraper(
        city=args.city,
        keyword=args.keyword,
        max_pages=args.pages,
        request_timeout=args.timeout,
        proxy=args.proxy or None,
    )

    cars = scraper.run()

    if cars:
        save_to_csv(cars, f"{args.output}.csv")
        save_to_json(cars, f"{args.output}.json")
        logger.info("爬取完成，共获取 %d 条二手车数据", len(cars))
    else:
        logger.warning("未获取到任何数据，请检查网络或目标网站结构是否变化")


if __name__ == "__main__":
    main()
