@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist "co2_hydrogenation_network\web\species_shape\index.html" (
  echo 找不到 co2_hydrogenation_network\web\species_shape\index.html，请确认仓库完整克隆。
  pause
  exit /b 1
)
start "" "co2_hydrogenation_network\web\species_shape\index.html"
