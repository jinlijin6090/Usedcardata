#!/usr/bin/env bash
# 一键运行 EC 化学空间硬枚举。
#
#   bash run.sh                         # 用同目录下的 config.ini
#   bash run.sh /path/to/config.ini     # 指定配置文件
#   bash run.sh config.ini --set space.max_heavy=10 --set write.derive_max_heavy=8
#   bash run.sh config.ini --stages estimate           # 只做规模估算
#   bash run.sh --selftest                              # 自检（约 1 分钟）
#
# 第一次运行会把 env/ecenv.tar.gz（或分卷 env/ecenv.tar.gz.part-*）解压成独立的 Python 环境，
# 不依赖服务器自带的 Python / conda，也不需要联网。中断后重新运行会从已完成的层继续。
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${EC_ENV_DIR:-$HERE/env/ecenv}"

unpack_env() {
  if [ -x "$ENV_DIR/bin/python" ] && [ -f "$ENV_DIR/.unpacked" ]; then return; fi
  local tarball="$HERE/env/ecenv.tar.gz"
  if [ ! -f "$tarball" ]; then
    if compgen -G "$HERE/env/ecenv.tar.gz.part-*" > /dev/null; then
      echo "[run.sh] 合并环境分卷 ..."
      cat "$HERE"/env/ecenv.tar.gz.part-* > "$tarball.tmp" && mv "$tarball.tmp" "$tarball"
    else
      echo "[run.sh] 找不到 env/ecenv.tar.gz（或分卷 env/ecenv.tar.gz.part-*），请按 README 放好环境包" >&2
      exit 1
    fi
  fi
  echo "[run.sh] 首次运行：解压 Python 环境到 $ENV_DIR（只做一次）..."
  local tmp="$ENV_DIR.unpack.$$"
  rm -rf "$tmp"
  mkdir -p "$tmp"
  tar -xzf "$tarball" -C "$tmp"
  ( set +u; source "$tmp/bin/activate"; conda-unpack ) > /dev/null
  if [ -d "$ENV_DIR" ] && [ ! -f "$ENV_DIR/.unpacked" ]; then rm -rf "$ENV_DIR"; fi
  if [ ! -d "$ENV_DIR" ]; then mv "$tmp" "$ENV_DIR"; else rm -rf "$tmp"; fi   # another job may have finished first
  touch "$ENV_DIR/.unpacked"
  "$ENV_DIR/bin/python" -c "import rdkit, numpy, duckdb; print('[run.sh] 环境就绪：RDKit', rdkit.__version__, '| numpy', numpy.__version__, '| DuckDB', duckdb.__version__)"
}

unpack_env
export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}" PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RDK_BUILD_THREAD_SAFE_MODULES=1

if [ "${1:-}" = "--selftest" ]; then
  exec bash "$HERE/tests/selftest.sh" "$ENV_DIR/bin/python"
fi

CONFIG="$HERE/config.ini"
if [ $# -gt 0 ] && [ "${1#--}" = "$1" ]; then CONFIG="$1"; shift; fi
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"
LOG_DIR="$(dirname "$CONFIG")"
LOG="$LOG_DIR/ec_space_$(date +%Y%m%d_%H%M%S).log"
echo "[run.sh] 配置：$CONFIG    日志：$LOG"
"$ENV_DIR/bin/python" -u -m ecspace.main --config "$CONFIG" "$@" 2>&1 | tee -a "$LOG"
