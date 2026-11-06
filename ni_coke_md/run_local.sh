#!/usr/bin/env bash
# Run cases on a workstation without a scheduler: N cases in parallel, T threads each.
#   bash run_local.sh [-j N] [-t T] [task ...]       e.g.  bash run_local.sh -j 4 -t 4 task1_flux
# Finished cases are skipped, interrupted ones resume from their checkpoint.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$ROOT/env/activate.sh" || exit 1
J=1; T="${OMP_NUM_THREADS:-4}"
while getopts "j:t:" o; do case $o in j) J=$OPTARG;; t) T=$OPTARG;; *) exit 1;; esac; done
shift $((OPTIND-1))
"$NICOKE_PY" "$ROOT/campaign/campaign.py" make "$@"
"$NICOKE_PY" "$ROOT/campaign/campaign.py" list "$@" | xargs -P "$J" -I{} bash "$ROOT/run_case.sh" {} "$T"
