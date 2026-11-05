@echo off
chcp 65001 >nul
cd /d "%~dp0"
where conda >nul 2>nul
if errorlevel 1 (
  echo 需要先安装 Miniconda（免费）：https://docs.conda.io/en/latest/miniconda.html
  echo 安装后在“开始”菜单打开 Anaconda Prompt，进入本文件夹，再运行：conda env create -f environment.yml
  pause
  exit /b 1
)
conda env list | findstr /b /c:"co2ec " >nul
if errorlevel 1 (
  echo 正在创建环境 co2ec（约 5–10 分钟，需要联网）…
  call conda env create -f environment.yml
) else (
  echo 环境 co2ec 已存在，正在按 environment.yml 更新…
  call conda env update -n co2ec -f environment.yml
)
echo 完成。现在可以双击“一键工具.bat”。
pause
