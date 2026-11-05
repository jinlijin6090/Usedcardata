@echo off
chcp 65001 >nul
cd /d "%~dp0"
set "PY="
for %%P in ("%USERPROFILE%\miniconda3\envs\co2ec\python.exe" "%USERPROFILE%\anaconda3\envs\co2ec\python.exe" "%USERPROFILE%\miniforge3\envs\co2ec\python.exe" "%LOCALAPPDATA%\miniconda3\envs\co2ec\python.exe" "%ProgramData%\miniconda3\envs\co2ec\python.exe" "%ProgramData%\anaconda3\envs\co2ec\python.exe") do (
  if not defined PY if exist "%%~P" set "PY=%%~P"
)
if defined PY goto run
where conda >nul 2>nul
if not errorlevel 1 (
  conda run -n co2ec --no-capture-output python co2_tools.py %*
  goto end
)
where python >nul 2>nul
if not errorlevel 1 (
  set "PY=python"
  goto run
)
echo 没有找到 Python。请先双击“一键安装环境.bat”（需要先装 Miniconda：https://docs.conda.io/en/latest/miniconda.html）。
goto end
:run
"%PY%" co2_tools.py %*
:end
pause
