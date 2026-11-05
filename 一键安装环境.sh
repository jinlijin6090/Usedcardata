#!/usr/bin/env bash
cd "$(dirname "$0")"
if ! command -v conda >/dev/null 2>&1; then
  echo "需要先安装 Miniconda：https://docs.conda.io/en/latest/miniconda.html"; exit 1
fi
if conda env list | grep -q "^co2ec "; then
  echo "环境 co2ec 已存在，正在更新…"; conda env update -n co2ec -f environment.yml
else
  echo "正在创建环境 co2ec（约 5–10 分钟，需要联网）…"; conda env create -f environment.yml
fi
echo "完成。现在可以运行 一键工具。"
