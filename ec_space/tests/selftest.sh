#!/usr/bin/env bash
# 自检：在临时目录里跑一遍完整流程（重原子 ≤ 6），核对闭壳层分子数与独立枚举的已知值，并检查数据库、查询和网页。
set -euo pipefail
PY="${1:-python3}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ec_selftest.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}"
echo "[selftest] 临时目录 $TMP"
"$PY" -u -m ecspace.main --config "$HERE/config.ini" --set run.out="$TMP/out" --set space.max_heavy=6 \
      --set space.allow_peroxide=no --set write.derive_max_heavy=6 --set viewer.max_heavy=6 --set run.workers=2 \
      --stages enumerate,derive,stats,viewer > "$TMP/log.txt" 2>&1 || { tail -40 "$TMP/log.txt"; echo "[selftest] 失败"; exit 1; }
"$PY" - "$TMP/out" <<'EOF'
import json, sys
from pathlib import Path
out = Path(sys.argv[1])
expected = [2, 5, 17, 77, 379, 2270]          # independent enumeration, no O-O bonds
got = [json.loads((out / "closed" / f"L{n}" / "DONE.json").read_text())["count"] for n in range(1, 7)]
assert got == expected, f"closed-shell counts {got} != {expected}"
st = json.loads((out / "stats.json").read_text())
assert st["species_total"] > 50_000 and st["reactions_total"] > 50_000, st
assert (out / "viewer" / "index.html").exists() and (out / "viewer" / "meta.js").exists()
from ecspace.query import connect
con = connect(out)
n_sp = con.execute("SELECT count(*) FROM species").fetchone()[0]
ring = con.execute("SELECT count(*) FROM reactions WHERE family='radical_anion_fragmentation' AND reactants='O=C1OCCO1|-1' AND products='[CH2]COC(=O)[O-]'").fetchone()[0]
assert ring == 1, "EC radical-anion ring opening missing"
print(f"[selftest] 通过：闭壳层分子 {got}，物种 {st['species_total']:,}（DuckDB 读到 {n_sp:,} 行），反应 {st['reactions_total']:,}；EC•⁻ 开环反应存在；网页已生成")
EOF
