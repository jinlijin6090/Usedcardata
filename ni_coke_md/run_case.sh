#!/usr/bin/env bash
# Run (or resume) one case.
#   bash run_case.sh runs/task1_flux/i20_s1 [threads] [max_hours]
# max_hours: stop cleanly at a checkpoint before the job's walltime (exit code 3 -> resubmit).
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$ROOT/env/activate.sh" || exit 1
CASE="$1"
THREADS="${2:-$OMP_NUM_THREADS}"
MAXH="${3:-}"
export OMP_NUM_THREADS="$THREADS"
ARGS=(--threads "$THREADS")
[ -n "$MAXH" ] && ARGS+=(--max-hours "$MAXH")
"$NICOKE_PY" "$ROOT/md/driver.py" "$CASE" "${ARGS[@]}" >> "$CASE/driver.out" 2>&1
rc=$?
echo "$rc" > "$CASE/last_exit"
case $rc in
  0) echo "$CASE finished";;
  3) echo "$CASE stopped at checkpoint (max_hours) - resubmit to continue";;
  *) echo "$CASE FAILED (exit $rc), see $CASE/driver.out";;
esac
exit $rc
