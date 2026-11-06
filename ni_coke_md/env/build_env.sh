#!/usr/bin/env bash
# (Re)build env/nicoke_env.tar.gz on any x86_64 Linux machine WITH internet access.
# You normally don't need this: the archive is shipped with the package.
#   bash env/build_env.sh
# Sources: conda-forge (Python 3.12 + conda-pack, via a standalone micromamba binary) and PyPI wheels
# (lammps 2025.7.22.4.0 incl. REAXFF+OPENMP, mpich, numpy, scipy, networkx, matplotlib, pillow).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
curl -fsSL -o mm.tar.bz2 https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.10.0-0.tar.bz2
mkdir mm && tar -xjf mm.tar.bz2 -C mm
export MAMBA_ROOT_PREFIX="$WORK/root"
./mm/bin/micromamba create -y -q -p "$WORK/nicoke_env" -c conda-forge --override-channels \
    "python=3.12" pip conda-pack
"$WORK/nicoke_env/bin/pip" install --no-cache-dir \
    lammps==2025.7.22.4.0 mpich numpy scipy networkx matplotlib pillow
"$WORK/nicoke_env/bin/conda-pack" -p "$WORK/nicoke_env" -o "$HERE/nicoke_env.tar.gz" \
    --ignore-missing-files --ignore-editable-packages -f
echo "built $HERE/nicoke_env.tar.gz"
