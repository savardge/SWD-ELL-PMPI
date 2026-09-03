# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Fortran 90/77 + MPI research code for probabilistic (Bayesian) 1-D multiphysics inversion of receiver function (RF), surface wave dispersion (SWD), Rayleigh wave ellipticity (ELL), and magnetotelluric (MT) data. Sampling is reversible-jump MCMC (trans-dimensional: the number of layers is unknown) with parallel tempering. Companion code for Shahsavari et al. (2025), arXiv:2510.15779.

There is no test suite, no linter, and no package manager. "Running" the code means launching an MPI sampler in a working directory full of input files and letting it run for hours-to-days.

## Build

Requires Linux + NVIDIA HPC SDK (`mpifort`) with LAPACK/BLAS. Will not build on macOS.

```bash
cd src
make clean
make          # -> src/bin/prjmh_temper_rf
```

- Compiler and flags live in `src/Makefile.compiler` (`FC`, `FCSERIAL`, `FCSWD`, `FCRS`). `PASH` points into `$NVCOMPILERS/Linux_x86_64/22.5/compilers/lib` — bump this for a different SDK version. `Makefile.compiler_terrawulf` is an Intel/MKL variant kept for a specific cluster; swap it in by copying over `Makefile.compiler`.
- The top-level `src/Makefile` recursively builds `raysum/` and `swd/` into `obj/libraysum.a` and `obj/libswd.a`, and `ray3d/` into loose `.o` files. `make clean` also cleans those subdirectories.
- Post-processing binaries `postpred` and `postlog` are **not** in the Makefile — build them with `./make_postpred.sh` (edit the hard-coded library paths first; the script also links `libgpell.a`, which does not exist here — see ELL note below).
- Compiled artifacts (`obj/*.o`, `bin/`, `*.mod`, `ray3d.o`, prebuilt `postpred`/`postlog`) are committed to the repo. `src/*.mod` is a symlink into `src/swd/`. Expect `git status` noise after any build.

## Run

The binary takes no arguments — everything comes from files in the **current working directory**:

```bash
cd example1_partial_coupling/MT_SWD_RF
mpirun -np 12 ../../src/bin/prjmh_temper_rf
# Ctrl+C to stop sampling (results written incrementally)
```

`RF.sh` in each example directory is the SLURM equivalent (12 tasks, module loads for gcc + nvhpc).

Number of MPI ranks = number of parallel-tempering chains (`NPTCHAINS = NTHREAD`), so `-np` must exceed `NPTCHAINS1` (the number of chains held at T=1). The paper's results used 12.

Example directories: `example1_partial_coupling/` and `example2_full_coupling/`, each with one subdirectory per data combination (`RF`, `SWD`, `MT`, `MT_SWD`, `SWD_RF`, `MT_SWD_RF`). Each subdirectory is a self-contained run: inputs, outputs, and previously produced figures all live together, so **re-running overwrites the committed results in place**.

### Plotting

MATLAB. Put `plotting_scripts/` on the MATLAB path, then run `rf_plot_rjhist_varpar2.m` (it hard-codes `filename='HON_sample.mat'` at the top; edit to point elsewhere). `plotting_scripts/` is a large shared personal toolbox from the group — most of it is unrelated to this project; the relevant entry points are `rf_plot_*`, `rf_read_parfile.m`, `batch_samples.m`, `batch_map.m`.

## Input / output file conventions

`filebase.txt` in the run directory holds two lines: the length of the file prefix, then the prefix itself (e.g. `3` / `HON`). **Every** other filename is `<prefix>_<something>`, built in `read_input.f90` (lines ~150–190 and ~790–860). The station geometry file is the one fixed name: `sample.geom`.

Key inputs: `<base>_parameter.dat`, `<base>_covparameter.dat`, `<base>_RF.txt` / `_SWD.dat` / `_ELL.dat` / `_MT.dat`, covariance matrices `<base>_Cd*.dat` / `_Cdi*.dat`, reference model `<base>_vel_ref.txt`, and (for synthetic tests) `<base>_map_voro_true.dat`.

Key outputs: `<base>_voro_sample.txt` (the posterior sample chain), `<base>_map_voro.dat` (MAP model), `<base>_mappred*.dat` / `_obs*.dat` (data fit), `<base>_stepsize.txt`, `<base>_seeds.log`, `<base>_RJMH.log`. `postpred` consumes `<base>_sample_postpred.dat` and writes `<base>_postpred.dat`.

### `<base>_parameter.dat` is strictly positional

`READPARFILE` in `read_input.f90` reads it with unlabelled sequential `READ(20,*)` statements — the trailing `!!` text on each line is just a comment. Inserting, removing, or reordering a line silently shifts every subsequent value. Any change must be mirrored in three places:

1. `READPARFILE` / `READCOVPARFile` in `src/read_input.f90`
2. `PRINTPAR2` (echoes the settings at startup)
3. `plotting_scripts/rf_read_parfile.m` (the MATLAB reader used by all plotting scripts)

The same applies to `<base>_covparameter.dat`.

Which datasets are inverted is controlled by the `I_RV` (`-1` = RF, `1` = radial+vertical waveforms), `I_SWD`, `I_ELL`, `I_MT` flags — this is the only difference between the `MT/`, `SWD_RF/`, `MT_SWD_RF/` … example directories.

## Architecture

**`rjmcmc_com.f90`** — the global state module. Defines `objstruc`, the single struct carrying an entire model state (Voronoi nodes `voro`, derived layer parameters `par`, per-dataset hierarchical std devs `sdpar*`, AR error parameters `arpar*`, and `Dobs*`/`Dpred*`/`Dres*` arrays for each data type), plus every tunable global read from the parameter files. Essentially all other files `USE RJMCMC_COM`, so adding a field to `objstruc` or a global here ripples widely.

**`prjmh_temper_rf.f90`** — the main program and MPI driver. Master rank (`src`) owns the chains and the sample file; slave ranks run Metropolis-Hastings steps. Contains the tempering ladder (`beta = 1/dTlog^(it-1)`, `NPTCHAINS1` chains at T=1), the RJMCMC moves (`BIRTH_FULL`/`DEATH_FULL`, `BIRTH_VARPAR`/`DEATH_VARPAR` for variable per-layer complexity, `BIRTH_SINGLE`/`DEATH_SINGLE`), the perturbation proposals (`PROPOSAL*`), chain exchange (`TEMPSWP_MH`), sample output (`SAVESAMPLE`), and time-based checkpointing (`TCHCKPT`).

**`alloc_obj.f90`** — allocates `objstruc` members *and* builds the MPI derived datatypes (`MAKE_MPI_STRUC_SP`) used to ship objects between ranks. These must stay in sync with the `objstruc` definition; a mismatch produces silent corruption rather than a compile error.

**`loglhood.f90`** — `LOGLHOOD` dispatches on the `I_*` flags to `LOGLHOOD_RF` / `LOGLHOOD_RV` / `LOGLHOOD_SWD` / `LOGLHOOD_ELL` / `LOGLHOOD_MT`; the joint likelihood is their sum, which is what makes this "multiphysics". Also holds `MAKE_CURMOD` (Voronoi nodes → layered velocity model) and `INTERPLAYER` (the variable-complexity layer logic). `loglhood2.f90` is an unused older variant, as are `backup.f`, `respknt_all_in_one.f`, and the `raysum_27Nov14/` / `raysum_28Oct14/` / `respkennett/` directories.

**Forward solvers** (each a separate legacy code, mostly fixed-form Fortran 77):
- `raysum/` — RF synthetics (default; `iraysum=1`)
- `ray3d/` — alternative RF/waveform solver (`iraysum=0`), also handles dipping layers
- `swd/` — Rayleigh/Love dispersion (`raydsp.f`, `lovdsp.f`, `dispersion.f90`)
- `MT1DforwardB3.f90` — 1-D MT impedance / apparent resistivity + phase

**Iterative data-covariance estimation** — `COVmatEst.f90`, `ZCOVmatEst.f90`, `UpdateCOV.f90`. The main program's outermost `DO` loop is over `icovIter`: sample, re-estimate the data error covariance from the residual ensemble, check convergence per dataset, repeat until `cov_converged` or `MAXcovIter`. Governed entirely by `<base>_covparameter.dat`. Setting all `Icov_iterUpdate_*` to 0 reduces this to the original single-pass fixed-covariance behaviour.

**Numerical Recipes support** — `nr.f90`, `nrtype.f90`, `nrutil.f90`, `svdcmp.f90`, `four1.f90`, `realft.f90`, `convlv.f90`, etc. Treat as vendored; don't refactor.

## Gotchas

- **ELL is not currently buildable.** Despite the repo name and the `I_ELL` parameter, `Makefile` contains no ELL/gpell targets (it is byte-identical to `Makefile_withoutELL`), the `gpell/` directory referenced by `config_sharedlib` and `make_postpred.sh` is absent, and both the solver call in `LOGLHOOD_ELL` (`loglhood.f90:1570`) and its initialisation (`prjmh_temper_rf.f90:204`) are commented out. Running with `I_ELL=1` will not work without restoring the missing gpell/Geopsy dependency.
- Helper Python scripts (`make_map_voro.py`, `makeFilesseis2.py`, `make_std.py`) and the MATLAB simulation scripts (`makesim.m`, `makesim_correlated_noise.m`) are one-off utilities with hard-coded filebases and parameters at the top — copies are duplicated into each example directory and drift from each other. Edit the copy in the directory you are working in.
- `.asv` files are MATLAB autosaves; ignore them.
