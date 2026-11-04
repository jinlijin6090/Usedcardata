#!/bin/bash
# SLURM 提交示例。用你自己的提交文件时，只需要把下面最后一行（bash .../run.sh ...）放进去。
#SBATCH -J ec_space
#SBATCH -N 1
#SBATCH -c 64                 # 核数：程序会自动用满分配到的核
#SBATCH --mem=200G            # 内存：每核约 1–3 GB（见 README 的资源表）
#SBATCH -t 72:00:00
#SBATCH -o ec_space_%j.log

bash /path/to/ec_space/run.sh /path/to/ec_space/config.ini
