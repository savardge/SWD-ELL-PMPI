#!/bin/bash
#SBATCH --job-name=basel_grp
#SBATCH --partition=shared-cpu,public-cpu
#SBATCH --time=01:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
# One rjMcMC inversion per array task: run dir = line (TASK_ID+1) of runs.list.
# The engine samples until killed and flushes the sample file every NKEEP
# rows, so a clean `timeout` at RUNSEC seconds loses at most one buffer.
# Submit from the campaign root on the cluster (the directory holding runs/):
#   sbatch --array=0-14 --output=logs/%x_%A_%a.out run_matrix.sbatch
# Partitions are named explicitly: bamboo also has *-arm partitions where this
# x86 build cannot run.
set -eu
module load foss/2025b
REPO=${REPO:-$HOME/Codes/SWD-ELL-PMPI}
BIN="$REPO/src/bin/prjmh_temper_rf"
ROOT=${SLURM_SUBMIT_DIR}
LIST=${LIST:-runs/runs.list}
RUNSEC=${RUNSEC:-3600}
RUN=$(sed -n "$(( SLURM_ARRAY_TASK_ID + 1 ))p" "$ROOT/$LIST")
[ -n "$RUN" ] || { echo "no run for task $SLURM_ARRAY_TASK_ID"; exit 1; }
cd "$ROOT/runs/$RUN"
FB=$(tail -1 filebase.txt | tr -d ' \n')
rm -f "${FB}_voro_sample.txt"
echo "$(date -Is) $RUN on $(hostname): $SLURM_NTASKS ranks, ${RUNSEC}s"
set +e
timeout --signal=INT --kill-after=30 "$RUNSEC" srun --mpi=pmix "$BIN" > "${FB}_run.log" 2>&1
rc=$?
set -e
ROWS=$(wc -l < "${FB}_voro_sample.txt" 2>/dev/null || echo 0)
echo "$(date -Is) $RUN finished (rc=$rc, 124 = timeout as intended): $ROWS sample rows"
[ "$ROWS" -gt 0 ] || { echo "NO SAMPLES -- see ${FB}_run.log"; tail -40 "${FB}_run.log"; exit 1; }
