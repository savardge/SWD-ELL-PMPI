# SWD-ELL-PMPI

Trans-dimensional (reversible-jump MCMC, parallel-tempered, MPI) Bayesian
inversion of **surface-wave dispersion (SWD) and Rayleigh-wave ellipticity
(ELL)** for 1-D shear-wave velocity structure, with **multi-mode Rayleigh
dispersion** (any set of mode branches, each on its own period grid) and an
optional **adjacent-layer contrast constraint**.

This is a fork of [pejman-sh86/RF-SWD-ELL-MT-PMPI](https://github.com/pejman-sh86/RF-SWD-ELL-MT-PMPI)
(Jan Dettmer's rjMcMC code line). Branch `multimode-raydsp`:

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
cd src
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
| `<base>_SWD.dat` | fundamental-mode curve: `period(s) phase_velocity(km/s)`, ascending period |
| `<base>_SWD_M<m>.dat` | higher-mode curves, one per mode number `m` listed in `MODE_OF` |
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
MODE_OF  0 2                   Rayleigh mode number of each curve slot (NMODE ascending
                               integers). Files are named by mode (_SWD.dat, _SWD_M2.dat);
                               "0 2" fits the fundamental + second higher mode with no R1.
                               Absent = 0 1 ... NMODE-1
IGRP     0                     0 = phase velocity (default), 1 = group velocity
SWD_SCAN 0.08 1.6 0.005 0.001  DISPER80 root-scan window cmin cmax and step dc in km/s,
                               optional overtone step dc_over (default dc/5). Default
                               2.0 6.5 0.05 (crustal); near-surface work needs the
                               values shown (give dc_over explicitly for bit-reproducible
                               runs across builds)
```

The n-th Rayleigh mode is found by counting sign changes of the DISPER80
secular function along the c-scan (`swd/raydsp.f`, `RAYDSPN`). A model that
cannot produce an observed mode at an observed period is rejected (dropping
the point instead would let the likelihood reward vanishing modes). Each
curve carries its own hierarchical sigma.

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
