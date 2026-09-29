# SWD-ELL-PMPI

Trans-dimensional (reversible-jump MCMC, parallel-tempered, MPI) Bayesian
inversion of **surface-wave dispersion (SWD) and Rayleigh-wave ellipticity
(ELL)** for 1-D shear-wave velocity structure, with **multi-mode Rayleigh
dispersion** (any set of mode branches, each on its own period grid) and an
optional **adjacent-layer contrast constraint**.

Repository: [savardge/SWD-ELL-PMPI](https://github.com/savardge/SWD-ELL-PMPI)
(renamed from RF-SWD-ELL-MT-PMPI), a fork of
[pejman-sh86/RF-SWD-ELL-MT-PMPI](https://github.com/pejman-sh86/RF-SWD-ELL-MT-PMPI)
(Jan Dettmer's rjMcMC code line). Relative to that upstream (all on `main`):

- the receiver-function (RF) and magnetotelluric (MT) machinery has been
  **removed** (raysum/ray3d forward codes, MT1D forward, their data paths,
  proposals, covariance iterations and parameter-file lines). The original
  joint RF-SWD-ELL-MT code is in the git history and in the upstream repo;
- the ellipticity path is kept as in the upstream code. NOTE: its forward
  call (`ellipticity_gpell`, geopsy's gpell library) is commented out in the
  upstream source, so ELL inversion currently needs that library to be wired
  back in;
- ported from the NVIDIA HPC SDK to **gfortran + Open MPI** (macOS arm64 and
  Linux); see `src/Makefile.compiler`.

Everything below the "Multi-mode SWD" heading was ported from
`receiver_rjmcmc_varpar_sourceinv_joint` (branch `multimode-raydsp`) and is
validated to reproduce it bit-for-bit (see Validation).

## Build

```
git clone git@github.com:savardge/SWD-ELL-PMPI.git
cd SWD-ELL-PMPI/src
make            # gfortran/mpif90; LAPACK via -framework Accelerate (macOS)
```
On Linux replace `LIB = -framework Accelerate` in `src/Makefile` by
`-llapack -lblas`. Flags are in `src/Makefile.compiler` (the NVHPC original is
kept as `Makefile.compiler.nvhpc`).

## Run

```
cd example1_partial_coupling/SWD
mpirun -np 12 ../../src/bin/prjmh_temper_rf
```
`-np` must exceed `NPTCHAINS1 + 1`. A run directory needs
`filebase.txt` (two lines: name length, name `<base>`) and:

| file | content |
|---|---|
| `<base>_parameter.dat` | 48 positional lines + optional keyword lines (below) |
| `<base>_covparameter.dat` | 35 lines: iterative covariance-estimation settings (SWD, ELL); set `ICOVest 0` and `Icov_iterUpdate_* 0` for a plain rjMcMC run |
| `<base>_SWD.dat` | fundamental-mode phase-velocity curve: `period(s) velocity(km/s)`, ascending period |
| `<base>_SWD_M<m>.dat` | higher-mode phase curves, one per mode number `m` listed in `MODE_OF` |
| `<base>_SWDG.dat`, `<base>_SWDG_M<m>.dat` | **group**-velocity curves (slots with `GRP_OF = 1`), same format |
| `<base>_sdSWD.dat`, `<base>_sdSWD_M<m>.dat`, `<base>_sdSWDG*.dat` | `ICOV_SWD = 3` only: per-point standard deviations (km/s), one file per curve slot named like the data files, same row order and count as the data file. The likelihood divides residuals by these AND by the curve's sampled `sdparSWD`, so `sdparSWD` becomes a dimensionless scale (set `sdmn`/`sdmx` accordingly, or `ISD_SWD = 0` to fix it at the `map_voro` value, typically 1). A missing or short file, or a non-positive sd, is fatal at read time (2026-09-08; the earlier build read mode slot 1 only and left other slots at zero) |
| `<base>_ELL.dat` | ellipticity data (if `I_ELL = 1`) |
| `<base>_vel_ref.txt` | reference Vs model when `I_VREF = 1` (node velocities are perturbations around it) |
| `<base>_map_voro.dat` | starting model: `k`, `NLMX*NPL` node triplets (depth km, dVs, dVpVs; unused slots 0), `sdparSWD(NMODE)`, `sdparELL(NMODE_ELL)`, `arparSWD(NMODE)`, `arparELL(NMODE_ELL)` |

Output `<base>_voro_sample.txt`: one row per kept sample =
`logL, logPr, tcmp, k, voro(NLMX*NPL), sdparSWD, sdparELL, arparSWD, arparELL,
acc_ratio, iaccept_bd, ireject_bd, iaccept_bds, chain, source_rank`.
Split burn-in per `source_rank` (last column).

### Parameter file (positional lines)

```
 1 IMAP        1 = predict data for the map_voro model and exit
 2 IMAGSCALE   1 = magnitude-scaled error model
 3 ENOS        1 = even-numbered order-statistics prior on node depths
 4 IPOIPR      1 = Poisson prior on k (rate = lambda)
 5 IAR         1 = autoregressive error model
 6 I_VARPAR    1 = variable layer complexity (trans-D)
 7 IBD_SINGLE  1 = birth/death for single parameters onto nodes
 8 I_SWD       1 = invert SWD
 9 I_ELL       1 = invert ELL
10 I_VREF      1 = perturbations around <base>_vel_ref.txt
11 I_VPVS      1 = sample Vp/Vs, -1 = Vp = 1.75 Vs
12 ISMPPRIOR   1 = sample the prior
13 ISETSEED    1 = fixed random-seed table
14 IEXCHANGE   1 = parallel-tempering exchange moves
15 NDAT_SWD    max number of data per SWD curve (array width)
16 NMODE       number of SWD curves
17 NDAT_ELL    number of ELL data
18 NMODE_ELL   number of ELL modes
19 NLMN        min number of nodes
20 NLMX        max number of nodes
21 ICHAINTHIN  chain thinning interval
22 NKEEP       samples buffered before each write
23 NPTCHAINS1  number of T = 1 chains
24 dTlog       tempering increment (T_i = dTlog^i)
25 lambda      Poisson-prior rate for k
26 hmx         max partition depth [km]
27 hmin        min layer thickness [km]
28 armxSWD     max AR prediction size, SWD
29 armxELL     max AR prediction size, ELL
30 TCHCKPT     checkpoint interval [s] (inert)
31 dVs         one-sided Vs prior half-width around the reference [km/s]
32 dVpVs       one-sided Vp/Vs prior half-width
33 sdmn        hierarchical-sigma prior lower bounds: SWD ELL [km/s]
34 sdmx        hierarchical-sigma prior upper bounds: SWD ELL
35 ISD_SWD     1 = sample hierarchical sigma of the SWD curves
36 ISD_ELL     1 = sample hierarchical sigma of the ELL curves
37 ICOV_SWD    SWD likelihood: 0 implicit sigma, 1 hierarchical sigma per curve, 2 Cdi file, 3 sd file
38 ICOV_ELL    ELL likelihood (same coding)
39-48          ELL_verbose ELL_prec I_ABS_ELL I_LOG10_ELL I_SAMPLING_TYPE_ELL I_SET_STEP_ELL STEP_SIZE_ELL I_SET_COUNT_ELL COUNT_ELL I_SET_RANGE_ELL
```
Trailing lines that are not keywords are ignored.

### Multi-mode SWD: keyword lines (anywhere after line 48)

```
DVSCON   0.100                 max |adjacent-layer dVs| in km/s: indicator prior evaluated
                               on the final layer stack BEFORE the forward call (Kennett
                               2023/2026 Seismica; BayHunter lvz/hvz parity). Absent/<0 = off
DVSMONO  0.050                 one-sided DVSCON: max ALLOWED adjacent-layer Vs DECREASE with
                               depth in km/s. 0 = strictly non-decreasing; a small tolerance
                               admits the metre-scale softening real boring logs show while
                               excluding a fast lid over a much slower layer. Same indicator
                               prior, checked before the forward. Absent/<0 = off. Also
                               enables the warm start (no LVZ can exist, so it is exact)
VP_BROCHER 1                   Vp from Vs by Brocher (2005, BSSA 95, eq. 9), valid 0 < Vs < 4.5 km/s.
                               The sampled Vp/Vs is then ignored (leave dVpVs tiny). Also applied
                               to the deep half-space tail. Default 0
RHO_BROCHER 1                  density from Vp by Brocher (2005, eq. 1, Nafe-Drake fit), Vp clamped
                               to the polynomial's 1.5-8.5 km/s range. Default 0 = the legacy
                               2.35 + 0.036 (Vp-3)^2 g/cc, a CRUSTAL relation that returns ~2.6 g/cc
                               for a 175 m/s soil (the flat density in vel_ref only feeds the tail)
MODE_OF  0 2                   Rayleigh mode number of each curve slot (NMODE integers
                               >= 0). Files are named by mode (_SWD.dat, _SWD_M2.dat);
                               "0 2" fits the fundamental + second higher mode with no R1.
                               Absent = 0 1 ... NMODE-1
IGRP     0                     0 = phase velocity (default), 1 = group velocity: the type
                               of every slot unless GRP_OF is given
GRP_OF   0 1                   velocity type per curve slot (NMODE values, 0 phase, 1 group).
                               "MODE_OF 0 0" + "GRP_OF 0 1" inverts R0 phase + R0 group
                               jointly, each with its own period grid, hierarchical sigma
                               and AR parameter (phase and group are measured independently
                               and carry different noise). Group slots read _SWDG*.dat /
                               _sdSWDG*.dat. The (mode, type) pairs must be unique; order is
                               free. A phase and a group slot of the same mode on the same
                               period grid share ONE root search: DISPER80 evaluates the
                               analytic group velocity (energy integrals) at every phase root
SDMN_SWD 0.3 0.3               per-slot lower / upper bounds of the hierarchical-sigma prior
SDMX_SWD 3.0 3.0               (NMODE values each; absent = the scalar sdmn/sdmx of lines
                               33/34 for every slot). The map_voro start value of every slot
                               must lie inside its bounds when ISD_SWD = 1 (fatal otherwise)
SWD_SCAN 0.08 1.6 0.005 0.001  DISPER80 root-scan window cmin cmax and step dc in km/s,
                               optional overtone step dc_over (default dc/5). Default
                               2.0 6.5 0.05 (crustal); near-surface work needs the
                               values shown (give dc_over explicitly for bit-reproducible
                               runs across builds)
SWD_WARM 1                     warm-started root scan (see below): 1 on, 0 off,
                               -1 (default) on iff DVSCON > 0 or DVSMONO >= 0
```

The n-th Rayleigh mode is found by counting sign changes of the DISPER80
secular function along the c-scan (`swd/raydsp.f`, `RAYDSPN`). A model that
cannot produce an observed mode at an observed period is rejected (dropping
the point instead would let the likelihood reward vanishing modes). Each
curve carries its own hierarchical sigma.

Limitation: `ICOV_SWD = 2` (inverse covariance from a file) and the iterative
covariance estimation (`Icov_iterUpdate_SWD = 1`) hold a single
`NDAT_SWD x NDAT_SWD` matrix for all curves and are refused when `NMODE > 1`.

### Warm-started root scan (speed)

DISPER80 brackets a root by scanning phase velocity from `cmin` in steps of
`dc` and counting sign changes, and the original code repeats that scan from
`cmin` for every period and every mode: at the overtone step (dc = 1 m/s)
about 850 propagator calls per period. Because a mode's phase velocity rises
with period, the previous period's root is a lower bound for the next one, so
the scan can start just below it (`SWD_WARM`, on by default when `DVSCON` is
active):

| | propagator calls | speed-up |
|---|---|---|
| fundamental | scan from cmin -> from the previous root | 2.3-3.2x |
| overtones | " | 7.4-7.6x |

The scan GRID is untouched (the start is the grid point `cmin + K*dc` below
the bound), so an accepted warm scan returns the same bracket, and therefore a
**bit-identical** root, as the cold scan. Three checks keep that true:

1. **sign-parity certificate** (in `RAYDSPN`): the secular function changes
   sign at every root, so the parity of `sign F(cstart)/sign F(cmin)` counts
   the roots skipped. It must be the mode number; if it is not, the routine
   reverts to the cold scan by itself, at the cost of one propagator call.
2. **no crossing / implausible jump** (in `dispersion`): the period is redone
   cold if the warm scan finds nothing above its start, or lands more than
   `CWFWD` above the previous root.
3. **monotonicity audit**: if the finished curve ever decreases with period --
   the signature of a model whose roots move down faster than the warm window
   -- the whole curve is recomputed cold.

Limit: a strongly inverse-dispersive model (a fast lid over a slow channel)
can move a root down past TWO roots at once, which parity cannot see. Such
models only exist when the adjacent-layer contrast is unconstrained, so with
`SWD_WARM -1` (the default) the warm start is enabled only when `DVSCON > 0`
or `DVSMONO >= 0` (a monotonic profile cannot have a channel at all);
`SWD_WARM 1` forces it on, `SWD_WARM 0` off.

**Root-scan step and high Vp/Vs.** The fundamental root is found by a sign-change
scan in steps of `dc` (`SWD_SCAN`). For soft, high-Vp/Vs stacks (saturated soil,
Vp/Vs 4-8) the default 5 m/s step can step over the root pair and the model is
rejected as if no mode existed: one 175/300/900 m/s test model was rejected at
`dc` 0.005 and fine at 0.002, with the same logL to 4 digits at 0.002 and 0.001.
On 60 posterior-like structures under Brocher Vp/density, 59 were valid at 0.005
and all 60 at 0.003, so use `SWD_SCAN 0.08 2.5 0.003 0.001` (1.7x the cost) for
any run that leaves the Vp/Vs = 2 regime. Validation of the Brocher path: the
layer stack the engine prints matches an independent implementation to four
digits (Vp to 0.0001 m/s, density exactly) and the predicted dispersion matches
disba on that stack to 0.00 m/s; with both keywords off the predicted data are
byte-identical to the previous binary (`bin/prjmh_temper_rf.pre_brocher`).

Validation (`swd/test_warm_driver.f90`, `tools/export_test_models.py`): warm
and cold curves compared on 15 000 models drawn from the real HVC posteriors
(v2 R0+R1+R2, v3 R0..R3, and the NLMX = 30 / hmin = 2 m run), all three modes
-- **zero differing periods, zero differing validity patterns**. On 4 000
models of the unconstrained fast-lid posterior the fundamental and first
overtone are still exact and the second overtone differs at 0.2 % of periods,
which is the regime the default switch excludes.

### Converting older inputs

`tools/convert_legacy_inputs.py {rfswdellmt|receiver} SRC_DIR DST_DIR`
converts run directories from the original 81-line RF-SWD-ELL-MT format or
from `receiver_rjmcmc_varpar_sourceinv_joint` (44/46-line) format.

## Validation (2026-09-02, macOS arm64, gfortran 16 + Open MPI 5)

- `example1_partial_coupling/SWD` IMAP: logL = 160.74617935123956, and the
  predicted curve is identical to the pre-port build and to the authors'
  shipped predictions.
- HVC dam inputs (masw-das campaign), IMAP logL identical to the last digit to
  `receiver_rjmcmc_varpar_sourceinv_joint/multimode-raydsp`: R0 only
  (-19.155037748234786), R0+R1+R2 with DVSCON (-260.50604563286140), R0+R2 via
  `MODE_OF 0 2` (-267.50804914725632), NLMX = 30 / hmin = 2 m
  (-260.50604563286140).
- Full-length seeded run (160k kept samples, HVC v2 R0+R1+R2 with DVSCON)
  against the receiver code's posterior on identical inputs: depth-median
  Vs offset 0.032 in 68%-half-width units (replicate gate 0.2), medians equal
  to 1 m/s at every depth, identical 68% bands, per-curve sigma medians
  24/39/25 m/s in both, identical logL and k medians. The random streams
  differ because the codes consume the RNG differently.

## Post-processing

MATLAB scripts of the original code are in `plotting_scripts/`
(`rf_plot_rjhist_varpar3.m` draws the interface-probability / Vs / Vp-Vs
panels). Python equivalents for the sample-file layout above live in the
masw-das repository (`scripts/rjmcmc_rjhist_panels.py`,
`scripts/rjmcmc_dam_posterior.py`).

## References

- Dettmer, J., Dosso, S. E., Holland, C. W. (2010-2015): trans-dimensional
  Bayesian inversion papers underlying this code.
- Kennett, B. L. N. (2023, 2026), Seismica: interacting waveguides and the
  representation of gradient structures with higher modes (basis of DVSCON).
