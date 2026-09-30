# CLAUDE.md

Guidance for Claude Code when working in this repository. `README.md` is the
authoritative user documentation (inputs, keywords, validation numbers); keep
the two consistent.

## What this is

Fortran 90/77 + MPI research code for trans-dimensional (reversible-jump
MCMC, parallel-tempered) Bayesian inversion of surface-wave dispersion (SWD:
multi-mode Rayleigh phase and/or group velocity) and Rayleigh ellipticity
(ELL) for 1-D Vs structure. Fork of Dettmer's RF-SWD-ELL-MT code with the RF
and MT machinery removed (see git history / upstream for those).

There is no unit-test suite, no linter, no package manager. "Running" the
code means launching the MPI sampler in a directory of input files.

## Build

gfortran + Open MPI, macOS arm64 or Linux (`src/Makefile.compiler`; the
original NVHPC flags are kept in `Makefile.compiler.nvhpc`):

```bash
cd src && make clean && make      # -> src/bin/prjmh_temper_rf
```

- `src/Makefile` picks the LAPACK library by platform: Accelerate on macOS,
  `-lflexiblas` otherwise (override with `make LIB="-llapack -lblas"`). On the
  UNIGE clusters use `module load foss/2025b` and build inside a SLURM job,
  never on a login node.
- `-fcheck=bounds -g` are on by default in `FFLAGS_COMMON` (cheap insurance);
  drop them for a ~20 % faster production build.
- `src/swd/` builds `obj/libswd.a` (DISPER80: `raydsp.f`, `raymrx.f`,
  `dispersion.f90`). `src/*.mod` are symlinks into `src/swd/`. Build products
  (`*.o`, `*.mod`, `*.a`, `src/bin/`, `src/obj/`) are gitignored.
- There are no `postpred`/`postlog` sources in this repo any more.

## Run

Everything comes from the current working directory (no arguments):

```bash
cd example1_partial_coupling/SWD
mpirun -np 12 ../../src/bin/prjmh_temper_rf     # Ctrl+C to stop; samples are flushed every NKEEP rows
```

`-np` must exceed `NPTCHAINS1 + 1`. Set `IMAP 1` (line 1 of the parameter
file) to predict the data of `<base>_map_voro.dat` and exit; this is the
regression check (`example1_partial_coupling/SWD` must give
`logL = 160.74617935123956`) and the way synthetic data are generated.

## Input / output conventions

`filebase.txt` holds the prefix length and the prefix `<base>`; every other
file is `<base>_<something>` (see the table in README). SWD curve *slots* are
identified by `(MODE_OF, GRP_OF)` = (Rayleigh mode number, 0 phase / 1
group); the file name of a slot is built by `SWD_SLOT_FILE` in
`read_input.f90` (`_SWD`, `_SWD_M<m>`, `_SWDG`, `_SWDG_M<m>`, and `_sd…`
for `ICOV_SWD = 3`).

### `<base>_parameter.dat`: 48 positional lines + keyword tail

`READPARFILE` reads the first 48 lines with unlabelled sequential `READ(20,*)`
(the trailing `!!` text is a comment); inserting or reordering a line
silently shifts every later value. Everything after line 48 is parsed as
`KEYWORD value(s)` lines (`DVSCON`, `DVSMONO`, `VP_BROCHER`, `RHO_BROCHER`,
`MODE_OF`, `GRP_OF`, `IGRP`, `SDMN_SWD`, `SDMX_SWD`, `SWD_SCAN`, `SWD_WARM`);
unknown lines are ignored. A new setting goes in the keyword tail, never as a
new positional line, and must be echoed in `PRINTPAR2` and documented in the
README keyword block. `<base>_covparameter.dat` (35 lines) is positional too.

`plotting_scripts/` is the group's MATLAB toolbox for the *original* 81-line
format; it does not read this repo's parameter files. Python post-processing
for the sample-file layout lives in the masw-das repository and in
`validation/basel_group/analyze.py`.

## Architecture

- **`rjmcmc_com.f90`** — global state: `objstruc` (one model state: Voronoi
  nodes `voro`, layer parameters `par`, per-slot `sdparSWD`/`arparSWD`,
  `Dobs/Dpred/DresSWD(NMODE,NDAT_SWD)`) plus every tunable read from the
  parameter files. Module globals (`MODE_OF`, `GRP_OF`, …) need no MPI
  broadcast: every rank reads the input files itself.
- **`prjmh_temper_rf.f90`** — main program / MPI driver: tempering ladder,
  RJMCMC moves, proposals (`PROPOSAL*`, `PROPOSAL_SDSWD`, `PROPOSAL_ARSWD`),
  chain exchange, `SAVESAMPLE`, `SAVEREPLICA` (IMAP output, one row per slot).
- **`alloc_obj.f90`** — allocates `objstruc` and builds the MPI derived type
  (`MAKE_MPI_STRUC_SP`); must stay in sync with the `objstruc` definition.
- **`loglhood.f90`** — `LOGLHOOD` = `LOGLHOOD_SWD` + `LOGLHOOD_ELL`.
  `LOGLHOOD_SWD` builds the layer stack (+ deep tail from `_vel_ref.txt`),
  applies the `DVSCON`/`DVSMONO` indicator priors, forwards every slot with
  `dispersion_cu` (one root search gives phase and group velocity; a phase
  and a group slot of the same mode on the same period grid share it),
  rejects any model that cannot produce an observed mode at an observed
  period, then sums the per-slot likelihoods (`ICOV_SWD` 0/1/2/3).
- **`swd/`** — DISPER80 (Saito). `RAYDSPN` finds the n-th mode by counting
  sign changes along a c-scan, optionally warm-started (see README).
- **`COVmatEst.f90` / `UpdateCOV.f90`** — iterative data-covariance
  estimation; single-curve only (refused when `NMODE > 1`).
- **`validation/`** — forward-model validation against disba
  (`compare_disba.py`, rerun after any change to `swd/`) and the synthetic
  Basel-1 group-velocity test harness (`basel_group/`).

## Gotchas

- **ELL is not buildable**: the gpell/Geopsy forward call is commented out
  in `LOGLHOOD_ELL`; `I_ELL = 1` will not work.
- The `SWD_SCAN` step is an accuracy setting, not a speed knob: a coarse
  `dc` can step over a root pair and reject valid models (README, "Root-scan
  step and high Vp/Vs").
- The one-sided priors `dVs`/`dVpVs` are widths around `<base>_vel_ref.txt`;
  a constant reference gives a uniform absolute prior.
- `.asv` files are MATLAB autosaves; ignore them.
