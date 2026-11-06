# source this file:  source env/activate.sh
# Puts the bundled Python + LAMMPS first on PATH and isolates it from any user/system Python setup.
_NICOKE_ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export NICOKE_ROOT="$(dirname "$_NICOKE_ENV_DIR")"
export NICOKE_PY="$_NICOKE_ENV_DIR/nicoke_env/bin/python"
if [ ! -x "$NICOKE_PY" ]; then
    echo "ERROR: bundled environment not installed. Run: bash $_NICOKE_ENV_DIR/install_env.sh" >&2
    return 1 2>/dev/null || exit 1
fi
export PATH="$_NICOKE_ENV_DIR/nicoke_env/bin:$PATH"
export LD_LIBRARY_PATH="$_NICOKE_ENV_DIR/nicoke_env/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONNOUSERSITE=1          # ignore ~/.local site-packages
unset PYTHONPATH PYTHONHOME
export MPLBACKEND=Agg
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-false}"
