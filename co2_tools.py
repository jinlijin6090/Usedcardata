"""CO2 加氢网络 · 一键工具（菜单）。由 一键工具.bat / 一键工具.command / 一键工具.sh 启动，也可以直接 python co2_tools.py。

  1 打开物种级网页（全量 C16）      不需要任何依赖
  2 打开物种级网页（推荐版）
  3 提取子网络（按提示输入起点、终点、多走几步、碳数上限），完成后自动打开结果网页
  4 查询一个物种（节点信息 + 它参与的全部反应，导出 CSV）
  5 检查环境
"""
from __future__ import annotations

import sys
import time
import webbrowser
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NET = ROOT / "co2_hydrogenation_network"
SCRIPTS = NET / "scripts"
WEB = {"1": NET / "web" / "species_full" / "index.html", "2": NET / "web" / "species_shape" / "index.html"}
DB = {"1": NET / "data" / "net_C16_species", "2": NET / "data" / "net_C16_shape_species"}
OUT = ROOT / "结果"


def ask(prompt, default=""):
    try:
        v = input(f"{prompt}{f'（回车 = {default}）' if default else ''}：").strip()
    except EOFError:
        v = ""
    return v or default


def open_page(p: Path):
    if not p.exists():
        print(f"找不到 {p}。请确认仓库是完整克隆的（网页数据在 co2_hydrogenation_network/web/ 下）。")
        return
    print("正在用浏览器打开：", p)
    webbrowser.open(p.as_uri())


def need(mods):
    miss = []
    for m in mods:
        try:
            __import__(m)
        except ImportError:
            miss.append(m)
    if miss:
        print(f"缺少依赖：{', '.join(miss)}。请先运行“一键安装环境”，或 pip install {' '.join('rdkit' if m == 'rdkit' else m for m in miss)}")
        return False
    return True


def choose_db():
    v = ask("用哪个数据库？1 = 全量 C16（最完整，首次读取 2–3 分钟、需约 8 GB 内存），2 = 推荐版（小、快）", "1")
    return DB.get(v, DB["1"])


def extract():
    if not need(["numpy"]):
        return
    sys.path.insert(0, str(SCRIPTS))
    import subnet_extract
    db = choose_db()
    print("物种可以写常用名（CO2、甲醇、乙烯、丙烯、苯、甲苯、对二甲苯、正十六烷、异辛烷……）或 SMILES；多个用空格分隔。")
    src = ask("起点", "CO2").split()
    dst = ask("终点（必填）", "").split()
    if not dst:
        print("没有输入终点。")
        return
    k = ask("比最少步数多走几步", "0")
    mc = ask("只用碳数不超过多少的物种（回车 = 不限；建议填目标的碳数）", "")
    OUT.mkdir(exist_ok=True)
    out = OUT / f"子网络_{'_'.join(dst)[:40]}_k{k}{'_C' + mc if mc else ''}_{time.strftime('%m%d_%H%M%S')}"
    args = [str(db), "--from", *src, "--to", *dst, "--extra", k, "--out", str(out)]
    if mc:
        args += ["--max-c", mc]
    try:
        subnet_extract.main(args)
    except SystemExit as e:
        if e.code not in (0, None):
            print(e.code)
        return
    page = out / "网页.html"
    print(f"\n结果在：{out}")
    if page.exists():
        open_page(page)
    else:
        print("子网络太大，没有生成网页；表格已导出（物种.csv、反应.csv）。")


def query():
    if not need(["duckdb"]):
        return
    sys.path.insert(0, str(SCRIPTS))
    import species_query as q
    db = choose_db()
    s = ask("物种（SMILES 或名称，例如 CC(C)=C、对二甲苯）", "")
    if not s:
        return
    names = {"对二甲苯": "Cc1ccc(C)cc1", "甲苯": "Cc1ccccc1", "苯": "c1ccccc1", "乙烯": "C=C", "丙烯": "C=CC", "甲醇": "CO",
             "CO2": "O=C=O", "正十六烷": "CCCCCCCCCCCCCCCC", "异辛烷": "CC(C)CC(C)(C)C"}
    s = names.get(s, s)
    con = q.connect(db)
    ids = q.node_ids(con, s)
    if not ids:
        print("没有找到这个物种（SMILES 需要是合法写法）。")
        return
    il = ",".join(map(str, ids))
    q.show(con.execute(f"SELECT node_id, smiles, kind, formula, states FROM nodes WHERE node_id IN ({il})"))
    OUT.mkdir(exist_ok=True)
    dst = OUT / f"反应_{s.replace('/', '_')[:40]}.csv"
    con.execute(f"COPY (SELECT id, array_to_string(families, ';') AS families, equation FROM reactions "
                f"WHERE list_has_any(reactants, [{il}]) OR list_has_any(products, [{il}]) ORDER BY id) TO '{dst}' (HEADER, DELIMITER ',')")
    n = sum(1 for _ in open(dst, encoding="utf-8")) - 1
    print(f"它参与的反应共 {n:,} 条，已导出到 {dst}")


def check():
    print("Python", sys.version.split()[0], sys.executable)
    for m in ("numpy", "rdkit", "duckdb"):
        try:
            mod = __import__(m)
            print(f"  {m:7s} OK {getattr(mod, '__version__', '')}")
        except ImportError:
            print(f"  {m:7s} 缺少（打开网页不需要；提取需要 numpy，查询需要 duckdb，按 SMILES 识别物种建议装 rdkit）")
    for k, p in WEB.items():
        print(f"  网页 {'全量' if k == '1' else '推荐版'}：{'OK' if p.exists() else '缺少'}  {p}")
    for k, p in DB.items():
        print(f"  数据库 {'全量' if k == '1' else '推荐版'}：{'OK' if (p / 'nodes.tsv.gz').exists() else '缺少'}  {p}")


def main():
    acts = {"1": lambda: open_page(WEB["1"]), "2": lambda: open_page(WEB["2"]), "3": extract, "4": query, "5": check}
    if len(sys.argv) > 1 and sys.argv[1] in acts:          # e.g. python co2_tools.py 1
        acts[sys.argv[1]]()
        return
    while True:
        print("\n==== CO2 加氢网络 · 一键工具 ====\n"
              "1 打开物种级网页（全量 C16）\n2 打开物种级网页（推荐版）\n3 提取子网络（起点、终点、多走几步）\n"
              "4 查询一个物种（全部反应导出 CSV）\n5 检查环境\n0 退出")
        c = ask("请选择", "0")
        if c == "0":
            return
        if c in acts:
            try:
                acts[c]()
            except KeyboardInterrupt:
                print("\n已中断。")
            except Exception as e:      # keep the menu alive
                print("出错：", e)
        else:
            print("请输入 0–5。")


if __name__ == "__main__":
    main()
