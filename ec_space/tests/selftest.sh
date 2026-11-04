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
      --set solvent.solvents=EC,DMC --stages enumerate,derive,solvent,stats,viewer > "$TMP/log.txt" 2>&1 || { tail -40 "$TMP/log.txt"; echo "[selftest] 失败"; exit 1; }
"$PY" - "$TMP/out" <<'EOF'
import json, sys
from pathlib import Path
out = Path(sys.argv[1])
expected = [2, 5, 17, 77, 379, 2270]          # independent enumeration, no O-O bonds
got = [json.loads((out / "closed" / f"L{n}" / "DONE.json").read_text())["count"] for n in range(1, 7)]
assert got == expected, f"closed-shell counts {got} != {expected}"
st = json.loads((out / "stats.json").read_text())
assert st["species_total"] > 50_000 and st["reactions_total"] > 50_000, st
assert (out / "viewer_EC+DMC" / "index.html").exists() and (out / "viewer_EC+DMC" / "meta.js").exists()
assert st["solvent"]["reactions_total"] > 50_000, st["solvent"]
from ecspace.query import connect
con = connect(out)
n_sp = con.execute("SELECT count(*) FROM species").fetchone()[0]
ring = con.execute("SELECT count(*) FROM reactions WHERE family='radical_anion_fragmentation' AND reactants='O=C1OCCO1|-1' AND products='[CH2]COC(=O)[O-]'").fetchone()[0]
assert ring == 1, "EC radical-anion ring opening missing"
sv = con.execute("SELECT count(*) FROM reactions WHERE family='solvent_h_abstraction' AND reactants='O=C1O[CH]CO1 + COC(=O)OC'").fetchone()[0]
sv += con.execute("SELECT count(*) FROM reactions WHERE family='solvent_h_abstraction' AND reactants='[CH2]OC(=O)OC + O=C1OCCO1'").fetchone()[0]
assert sv == 1, f"EC/DMC cross H-abstraction should appear exactly once, got {sv}"
me = con.execute("SELECT count(*) FROM reactions WHERE family='solvent_deprotonation' AND reactants='C[O-] + COC(=O)OC' AND products='CO + [CH2-]OC(=O)OC'").fetchone()[0]
assert me == 1, "deprotonation of DMC by CH3O- missing"
print(f"[selftest] 通过：闭壳层分子 {got}，物种 {st['species_total']:,}（DuckDB 读到 {n_sp:,} 行），反应 {st['reactions_total']:,}；EC•⁻ 开环反应存在；EC/DMC 溶剂反应 {st['solvent']['reactions_total']:,} 条（交叉抽 H 只记一次、CH₃O⁻ 夺 DMC 质子存在）；网页已生成")
EOF
