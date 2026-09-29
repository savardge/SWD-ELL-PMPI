#!/bin/bash
#SBATCH --job-name=swdell_imap
#SBATCH --partition=debug-cpu
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
# Linux-build regression: example1_partial_coupling/SWD in IMAP mode must give
# logL = 160.74617935123956 (README) within 1e-4 (cross-platform single-precision
# noise is ~1e-5). Runs in a scratch copy so the committed
# example is untouched.   sbatch --output=logs/%x_%j.out imap_check.sbatch
set -eu
module load foss/2025b
REPO=${REPO:-$HOME/Codes/SWD-ELL-PMPI}
WORK=${SCRATCH_ROOT:-/srv/beegfs/scratch/users/s/savardg/swd_group_basel}/imap_check
rm -rf "$WORK"; mkdir -p "$WORK"
cp "$REPO/example1_partial_coupling/SWD/"* "$WORK/"
cd "$WORK"
sed -i '1s/^0 /1 /' HON_parameter.dat
srun --mpi=pmix "$REPO/src/bin/prjmh_temper_rf" > imap.log 2>&1 || true
L=$(grep -E "^ logL " imap.log | head -1 | awk -F= '{print $2}')
echo "logL = $L (Mac arm64 reference 160.74617935123956; x86 GCC 14 gives 160.74615570403239: single-precision forward-model noise, 1.5e-7 relative)"
if awk -v l="$L" 'BEGIN{d=l-160.74617935123956; if (d<0) d=-d; exit !(d<1e-4)}'; then echo "IMAP REGRESSION PASS"; else echo "IMAP REGRESSION FAIL"; tail -30 imap.log; exit 1; fi
