#!/usr/bin/env bash
cd "$(dirname "$0")"
PY=""
for p in "$HOME/miniconda3/envs/co2ec/bin/python" "$HOME/anaconda3/envs/co2ec/bin/python" "$HOME/miniforge3/envs/co2ec/bin/python" \
         "/opt/homebrew/Caskroom/miniconda/base/envs/co2ec/bin/python" "/opt/miniconda3/envs/co2ec/bin/python"; do
  [ -z "$PY" ] && [ -x "$p" ] && PY="$p"
done
if [ -z "$PY" ] && command -v conda >/dev/null 2>&1 && conda env list | grep -q "^co2ec "; then
  conda run -n co2ec --no-capture-output python co2_tools.py "$@"; exit $?
fi
[ -z "$PY" ] && PY="$(command -v python3 || command -v python)"
if [ -z "$PY" ]; then echo "没有找到 Python。请先运行 一键安装环境（需要先装 Miniconda）。"; exit 1; fi
"$PY" co2_tools.py "$@"
