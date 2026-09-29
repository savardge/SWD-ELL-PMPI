# Synthetic Basel-1 / Otterbach-2 tests of the group-velocity inversion

Tests the per-slot velocity type (`GRP_OF`: phase and/or group velocity per
Rayleigh mode, README keyword block) on synthetic Rayleigh dispersion curves
forward-modelled from the Basel-1 / Otterbach-2 well model of Michel et al.
(2017) (`/Volumes/T7blue/riehen-data/well-data/Michel2016_gpdc.model`: 17
layers to 3 km, Vs 0.78-3.37 km/s, several low-velocity zones) at 30
log-spaced periods between 0.5 and 6 s.

## Pipeline

```bash
PY=/opt/anaconda3/envs/bayesbay_dev/bin/python        # numpy, matplotlib, disba
$PY make_synthetics.py     # data/: curves + sd, truth nodes, checks.txt; figures/fig_forward_check.png
$PY stage_runs.py          # runs/<name>/ (15 engine run dirs) + runs/runs.list
# ---- bamboo (UNIGE): never run anything on a login node ----
git push; ssh bamboo 'cd ~/Codes/SWD-ELL-PMPI && git pull'
ssh bamboo 'cd ~/Codes/SWD-ELL-PMPI/validation/basel_group && sbatch --output=logs/%x_%j.out build.sbatch'
ssh bamboo 'cd ~/Codes/SWD-ELL-PMPI/validation/basel_group && sbatch --output=logs/%x_%j.out imap_check.sbatch'   # logL 160.74617935123956
rsync -a runs/ bamboo:/srv/beegfs/scratch/users/s/savardg/swd_group_basel/runs/
scp run_matrix.sbatch bamboo:/srv/beegfs/scratch/users/s/savardg/swd_group_basel/
ssh bamboo 'cd /srv/beegfs/scratch/users/s/savardg/swd_group_basel && sbatch --array=0-14 --output=logs/%x_%A_%a.out run_matrix.sbatch'
# ---- back on the Mac ----
rsync -a bamboo:/srv/beegfs/scratch/users/s/savardg/swd_group_basel/runs/ runs/
$PY analyze.py             # runs/<name>/diagnostics.png, results/summary.csv, figures/overlay_*.png
```

## Design

- **Truth in model space.** The Michel model is written as 17 engine nodes
  (depth, `dVs = Vs - 2.1`, `dVpVs = Vp/Vs - 1.85`) around a constant
  reference (`vel_ref`: Vs 2.1 km/s, Vp/Vs 1.85, tail rows at 5, 8, 12 km at
  the same values, so the half-space continues the deepest layer). Prior:
  Vs 0.7-3.5 km/s, Vp/Vs 1.70-2.00, `hmx` 5 km, `hmin` 20 m, `NLMX` 30,
  Poisson(10) on k. The engine's own density relation applies, so the
  synthetic "engine" curves are exactly reproducible by the sampler.
- **Noise**: Gaussian, 2 % of the datum on phase, 4 % on group (phase and
  group are measured independently), one realisation (seed 20260929) shared
  by every run; `sd` files carry it, so with `ICOV_SWD 3` the sampled sigma
  of each slot is a unit-scale multiplier (prior 0.3-3).
- **Scan**: `SWD_SCAN 0.6 3.6 0.005 0.001`, cold (the truth has LVZs, so no
  `DVSCON`/`DVSMONO` and no warm start).

## Matrix (`stage_runs.py`)

| run | slots | periods |
|---|---|---|
| `R0p_*`, `R0g_*`, `R0pg_*`, `R01pg_*` | R0 phase; R0 group; both; + R1 phase and group | `full` 0.5-6 s (30), `short` 0.5-2 s (17), `long` 2-6 s (13) |
| `R0pg_mismatch` | R0 phase 0.5-6 s + R0 group 1-4 s | different grids: no root-solve reuse |
| `R0pg_modelerr` | R0 phase + group | curves from the literal Michel stack (true Vp, Gardner density), outside the model space |
| `R0pg_full_icov1` | R0 phase + group | `ICOV_SWD 1`, `IMAGSCALE 1`, per-slot `SDMN_SWD`/`SDMX_SWD` instead of sd files |

Each run: 12 MPI ranks (6 chains at T = 1, `dTlog` 1.15), one hour, one seed.

## Forward-model checks (`data/checks.txt`, `figures/fig_forward_check.png`)

Every one of these must PASS before the inversions are meaningful:

- the engine in IMAP mode reproduces the `disp_driver` curves of the truth
  slot by slot (phase and group, R0 and R1) to single precision, and reports
  that the group slots reuse the phase slots' root solve;
- DISPER80 phase velocity vs disba 0.7.0 (roots matched by velocity, since
  disba duplicates roots on this LVZ model and its mode indices shift):
  <= 0.02 m/s for R0 and R1 on both stacks;
- the analytic group velocity vs a finite difference of our own phase curve:
  <= 1 m/s (0.2 m/s typical; the R1 cut-off period is excluded, the branch
  bends too sharply there for a 0.2 % finite difference). disba's
  `GroupDispersion` is reported for information only: it disagrees by
  hundreds of m/s on this model, as it did in `validation/` (fig 5);
- our R0 **group** velocity vs the gpdc curve shipped with the model
  (`Michel2016_gpdc.disp`, 0.2-4.3 s, an independent code): 0.02 % median,
  0.08 % max -- which also settles that the file holds group slowness.

## Results

`analyze.py` writes one `diagnostics.png` per run (posteriors with the true
profile, convergence per chain, data fits per slot) and
`results/summary.csv`; see the README section "Group-velocity synthetic test"
for the outcome of the campaign.
