#!/usr/bin/env bash
# Submit every unfinished case that is not already queued/running. Safe to run repeatedly.
#
#   bash pbs/submit.sh                     # all tasks (task2 only once its task1 inputs exist)
#   bash pbs/submit.sh task1_flux task3_codep
#   bash pbs/submit.sh --dry-run           # show the qsub commands only
#
# Cases are packed CORES_PER_NODE/THREADS_PER_CASE per job (one node), longest first, so cases of
# similar length share a node. Each job stops its cases at a checkpoint before the walltime and
# resubmits the unfinished ones itself (pbs/job.pbs).
# Internal use by job.pbs:  --caselist FILE --resub N   |   --if-ready (silent when nothing to do)
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/env/activate.sh" || exit 1
source "$ROOT/pbs/pbs.conf"

DRY=0; CASELIST=""; RESUB=0; QUIET=0; TASKS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) DRY=1;;
        --caselist) CASELIST="$2"; shift;;
        --resub) RESUB="$2"; shift;;
        --if-ready) QUIET=1;;
        -h|--help) sed -n '2,13p' "$0"; exit 0;;
        *) TASKS+=("$1");;
    esac
    shift
done
say() { [ "$QUIET" = 1 ] || echo "$@"; }

JOBDIR="$ROOT/runs/.jobs"
mkdir -p "$JOBDIR"

# ---- one submitter at a time (several jobs may finish at once)
LOCK="$ROOT/runs/.submit.lock"
for i in $(seq 120); do
    mkdir "$LOCK" 2>/dev/null && break
    if [ -n "$(find "$LOCK" -maxdepth 0 -mmin +10 2>/dev/null)" ]; then rmdir "$LOCK" 2>/dev/null; fi
    sleep 1
done
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

# ---- walltime
wt="$WALLTIME"
if [ -z "$wt" ] && command -v qstat >/dev/null; then
    wt=$(qstat -Qf "$QUEUE" 2>/dev/null | awk -F= '/resources_max.walltime/{gsub(/ /,"",$2); print $2; exit}')
fi
[ -z "$wt" ] && wt="24:00:00"
WALL_H=$(echo "$wt" | awk -F: '{ if (NF==3) print $1+$2/60+$3/3600; else if (NF==2) print $1+$2/60; else print $1/3600 }')

queued() {   # is the job id stored in CASE/.pbs_job still known to PBS?
    local f="$1/.pbs_job" jid
    [ -f "$f" ] || return 1
    jid=$(cat "$f")
    [ -n "$jid" ] && [ "$jid" != "local" ] && qstat "$jid" >/dev/null 2>&1
}

# ---- collect cases: "<remaining_ps> <dir>" per line, grouped by task
declare -A GROUP
if [ -n "$CASELIST" ]; then
    while read -r d; do
        [ -n "$d" ] && GROUP[resub]+="0 $d"$'\n'
    done < "$CASELIST"
else
    [ ${#TASKS[@]} -eq 0 ] && TASKS=(task1_flux task2_anneal task3_codep)
    "$NICOKE_PY" "$ROOT/campaign/campaign.py" make "${TASKS[@]}" >/dev/null
    for t in "${TASKS[@]}"; do
        nw=0
        while read -r rem ready d; do
            [ -z "$d" ] && continue
            if [ "$ready" != 1 ]; then nw=$((nw+1)); continue; fi
            queued "$d" && continue
            GROUP[$t]+="$rem $d"$'\n'
        done < <("$NICOKE_PY" "$ROOT/campaign/campaign.py" pending "$t")
        [ $nw -gt 0 ] && say "$t: $nw case(s) waiting for their input (task1 not finished yet)$( [ "$AUTO_TASK2" = 1 ] && echo ' - will be submitted automatically')"
    done
fi

SLOTS=$(( CORES_PER_NODE / THREADS_PER_CASE ))
STAMP=$(date +%Y%m%d_%H%M%S)
n_sub=0
for g in "${!GROUP[@]}"; do
    mapfile -t sorted < <(printf '%s' "${GROUP[$g]}" | sort -k1,1 -g -r | awk '{print $2}')
    [ ${#sorted[@]} -eq 0 ] && continue
    for ((k=0; k<${#sorted[@]}; k+=SLOTS)); do
        chunk=("${sorted[@]:k:SLOTS}")
        tag=$( [ "$g" = resub ] && basename "$(dirname "${chunk[0]}")" || echo "$g" )
        list="$JOBDIR/${tag}_${STAMP}_$((k/SLOTS)).list"
        printf '%s\n' "${chunk[@]}" > "$list"
        nc=$(( ${#chunk[@]} * THREADS_PER_CASE ))
        [ "$FULL_NODE" = 1 ] && nc=$CORES_PER_NODE
        args=(-N "nicoke_${tag#task}" -q "$QUEUE" -l "select=1:ncpus=${nc}:mpiprocs=1:ompthreads=${THREADS_PER_CASE}"
              -l "walltime=$wt" -j oe -o "$JOBDIR/"
              -v "NICOKE_ROOT=$ROOT,CASELIST=$list,WALL_H=$WALL_H,RESUB=$RESUB")
        [ -n "$PBS_SOFTWARE" ] && args+=(-l "software=$PBS_SOFTWARE")
        # shellcheck disable=SC2206
        [ -n "$EXTRA_QSUB" ] && args+=($EXTRA_QSUB)
        if [ "$DRY" = 1 ]; then
            echo "qsub ${args[*]} $ROOT/pbs/job.pbs"
            sed 's#^#    #' "$list"
            continue
        fi
        jid=$(qsub "${args[@]}" "$ROOT/pbs/job.pbs") || { echo "qsub failed for $list" >&2; continue; }
        for d in "${chunk[@]}"; do echo "$jid" > "$d/.pbs_job"; done
        echo "submitted $jid: ${#chunk[@]} case(s), ncpus=$nc, walltime=$wt  [$list]"
        n_sub=$((n_sub+1))
    done
done
[ "$DRY" = 1 ] || [ $n_sub -gt 0 ] || say "nothing to submit (all cases finished, queued, or waiting for inputs)"
