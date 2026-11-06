#!/usr/bin/env bash
# Check the bundled environment and measure ReaxFF speed on this machine.
#   bash env/check_env.sh [threads]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/activate.sh"
THREADS="${1:-$OMP_NUM_THREADS}"
"$NICOKE_PY" - "$THREADS" <<'EOF'
import sys, os, time, tempfile
import numpy, scipy, networkx, matplotlib, PIL
from lammps import lammps
th = int(sys.argv[1])
root = os.environ["NICOKE_ROOT"]
L = lammps(cmdargs=["-log", "none", "-screen", "none", "-nocite"] + (["-sf", "omp", "-pk", "omp", str(th)] if th > 1 else []))
pk = L.installed_packages
print(f"   python {sys.version.split()[0]}, numpy {numpy.__version__}, scipy {scipy.__version__}, "
      f"networkx {networkx.__version__}, matplotlib {matplotlib.__version__}")
print(f"   LAMMPS {L.version()}  REAXFF={'REAXFF' in pk}  OPENMP={'OPENMP' in pk}")
assert "REAXFF" in pk and "OPENMP" in pk
for c in ["units real", "atom_style charge", "boundary p p p",
          f"read_data {root}/structures/ni116.data",
          "pair_style reaxff NULL checkqeq yes safezone 3.0 mincap 200",
          f"pair_coeff * * {root}/ff/ffield.reax.FC C Ni O",
          "fix qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff", "timestep 0.25",
          "velocity all create 1200 1 mom yes rot yes", "fix nvt all nvt temp 1200 1200 50",
          "run 200"]:
    L.command(c)
t = time.time(); L.command("run 2000"); dt = time.time() - t
print(f"   ReaxFF Ni116 benchmark: {2000/dt:.0f} steps/s with {th} thread(s) "
      f"-> {4000/(2000/dt):.1f} s per ps (bare cluster; expect ~2x slower once ~120 C are deposited)")
print("OK")
EOF
