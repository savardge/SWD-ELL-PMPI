#!/bin/bash
#SBATCH --job-name=swdell_build
#SBATCH --partition=debug-cpu
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
# Build SWD-ELL-PMPI on bamboo (UNIGE) INSIDE a job -- never on a login node.
#   sbatch --output=logs/%x_%j.out build.sbatch
# Then the Linux regression: sbatch imap_check.sbatch (must print logL = 160.74617935123956).
set -eu
module load foss/2025b
REPO=${REPO:-$HOME/Codes/SWD-ELL-PMPI}
cd "$REPO/src"
echo "$(date -Is) building in $REPO/src on $(hostname) with $(mpif90 --version | head -1)"
make clean >/dev/null
make 2>&1 | grep -Ev "^(gfortran|mpif90|ar |a - |r - |mkdir|/bin/ln)" || true
ls -la bin/prjmh_temper_rf
# the standalone forward driver used by the harness (not needed on the cluster, but cheap)
cd "$REPO/validation"
gfortran -O2 -fallow-argument-mismatch -std=legacy -fdec-format-defaults -w -fno-realloc-lhs \
    -o disp_driver disp_driver.f90 ../src/swd/libswd.a
echo "$(date -Is) done"
