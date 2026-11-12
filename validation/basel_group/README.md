# Synthetic Basel-1 / Otterbach-2 tests of the group-velocity inversion

Tests the per-slot velocity type (`GRP_OF`: phase and/or group velocity per
Rayleigh mode, README keyword block) on synthetic Rayleigh dispersion curves
forward-modelled from the Basel-1 / Otterbach-2 well model of Michel et al.
(2017) (`Michel2016_gpdc.model`, not distributed; set `SWD_WELL_DIR` to its folder: 17
layers to 3 km, Vs 0.78-3.37 km/s, several low-velocity zones) at 30
log-spaced periods between 0.5 and 6 s.

## Pipeline

```bash
PY=python                  # an environment with numpy, matplotlib, disba
$PY make_synthetics.py     # data/: curves + sd, truth nodes, checks.txt; figures/fig_forward_check.png
$PY stage_runs.py          # runs/<name>/ (15 engine run dirs) + runs/runs.list
# ---- cluster (SCRATCH_ROOT = your scratch run dir): never run anything on a login node ----
git push; ssh bamboo 'cd SWD-ELL-PMPI && git pull'
ssh bamboo 'cd SWD-ELL-PMPI/validation/basel_group && sbatch --output=logs/%x_%j.out build.sbatch'
ssh bamboo 'cd SWD-ELL-PMPI/validation/basel_group && sbatch --output=logs/%x_%j.out imap_check.sbatch'   # logL 160.74617935123956
rsync -a runs/ bamboo:$SCRATCH_ROOT/runs/
scp run_matrix.sbatch bamboo:$SCRATCH_ROOT/
ssh bamboo 'cd $SCRATCH_ROOT && sbatch --array=0-14 --output=logs/%x_%A_%a.out run_matrix.sbatch'
# ---- back on the workstation ----
rsync -a bamboo:$SCRATCH_ROOT/runs/ runs/
$PY analyze.py             # runs/<name>/diagnostics.png, results/summary.csv, figures/overlay_*.png
                           # (built on tools/swdell; any run dir: python tools/plot_rjhist.py RUN_DIR)
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
profile, convergence per chain, data fits per slot) and `results/summary.csv`.
15 runs of 1 h x 12 ranks on bamboo (2026-09-29).

- **The implementation is correct end to end.** Every run fits its data with
  reduced chi-square 0.6-0.9 per slot, the hierarchical sigma multiplier of
  every slot comes back at 0.93-1.17 (expected 1), `R0pg_full_icov1`
  (`ICOV_SWD 1`, `IMAGSCALE 1`, per-slot bounds) recovers 1.95 % / 3.77 %
  (truth 2 / 4 %), and `R0pg_mismatch` (group on a different grid, no
  root-solve reuse) behaves like `R0pg_full`. The engine's
  IMAP predictions equal `disp_driver`'s slot by slot, and the analytic group
  velocity matches the gpdc curve shipped with the model to 0.02 %.
- **What group velocity adds (measured, full band).** Mean 68 % half-width
  of Vs in 0-1 km: R0 phase 0.172, R0 group 0.165, R0 phase + group
  0.141 km/s (1-3 km: 0.249 / 0.214 / 0.194); mean |median - truth| in
  0-1 km 0.34 / 0.32 / 0.32 km/s. So group velocity narrows the bands by
  ~20 % and the joint run is the tightest, but the posterior median hardly
  moves: the extra information sharpens the average structure, it does not
  recover the 26-100 m layers. All R0 runs agree from 1 km down (the
  2.0-2.4 km/s section is smoothed, the 3 km half-space is exact). Adding R1
  phase + group starts to resolve the 1.1-1.6 km layering (5/16 true
  interfaces found by a node-density heuristic vs 2-3), but those runs are
  6-8x slower (overtone scan) and were not converged in 1 h (stationarity
  L1 0.21; re-run for 4-6 h). `R0pg_mismatch` (group on 1-4 s only) is
  consistent with `R0pg_full` with a wider band (0.21 km/s in 0-1 km).
- **Resolution, not the sampler, limits the recovery -- and the bands are
  over-confident.** k median 7-9 vs the true 17: the fine layers of the log
  are below what 0.5-6 s dispersion can see, and the trans-D posterior
  replaces them by their average. Measured against the blocky truth, the
  68 % bands contain it at only 13-31 % of depths in 0-1 km and 36-54 % in
  1-3 km (all runs): the independent-noise likelihood with a hierarchical
  sigma treats the unresolved layering as absent, not as uncertainty, the
  usual rjMcMC over-confidence. Interpret the band as the uncertainty of the
  *smoothed* profile.
- **Model error.** With curves computed from the literal Michel stack (true
  Vp and Gardner density, up to 0.28 km/s different from the engine's own
  density relation), the sampler absorbs the difference into Vs (the profile
  shifts by less than the band width) rather than into sigma (0.94 / 0.96):
  the density relation matters at the percent level on this model.
- Convergence in 1 h: stationarity / between-chain / burn-in L1 metrics
  0.02-0.11 for every one- and two-slot run (<= 0.05 comfortable, >= 0.15
  re-run), 0.21 for the four-slot runs.
