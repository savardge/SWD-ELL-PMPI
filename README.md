# SWD-ELL-PMPI

Trans-dimensional (rjMcMC, parallel-tempered, MPI) Bayesian inversion of
surface-wave dispersion (SWD) and Rayleigh ellipticity (ELL) for 1-D Vs, with
multi-mode phase and group curves (each on its own period grid) and optional
adjacent-layer contrast constraints.

Fork of [pejman-sh86/RF-SWD-ELL-MT-PMPI](https://github.com/pejman-sh86/RF-SWD-ELL-MT-PMPI)
(Jan Dettmer's code line). Changes: receiver-function and MT paths removed
(still in the git history); ported from NVHPC to gfortran + Open MPI; the
keywords below added. `MODE_OF`, `DVSCON` and `SWD_SCAN` were ported from
`receiver_rjmcmc_varpar_sourceinv_joint` and reproduce it bit-for-bit; the rest
is new here. The ELL forward call (geopsy `gpell`) is commented out upstream,
so ELL inversion needs that library wired back in.

## References

Original code (cite when using it):

- Shahsavari, P., Dettmer, J., Unsworth, M. J., & Schaeffer, A. (2025). Multiphysics inversion with variable complexity of receiver-function, surface-wave dispersion and magnetotelluric data reduces uncertainty for lithosphere structure. arXiv:2510.15779. https://doi.org/10.48550/arXiv.2510.15779
- Shahsavari, P., Dettmer, J., Unsworth, M. J., & Schaeffer, A. (2022). Data and model appraisal for joint inversion of magnetotelluric, receiver function, Rayleigh wave dispersion, and horizontal to vertical amplitude ratio (HV) data. AGU Fall Meeting 2022, S15B-03.

Methods:

- Dettmer, J., Molnar, S., Steininger, G., Dosso, S. E., & Cassidy, J. F. (2012). Trans-dimensional inversion of microtremor array dispersion data with hierarchical autoregressive error models. *Geophys. J. Int.*, 188(2), 719–734. https://doi.org/10.1111/j.1365-246X.2011.05302.x
- Dettmer, J., & Dosso, S. E. (2012). Trans-dimensional matched-field geoacoustic inversion with hierarchical error models and interacting Markov chains. *J. Acoust. Soc. Am.*, 132(4), 2239–2250. https://doi.org/10.1121/1.4746016
- Dettmer, J., Dosso, S. E., Bodin, T., Stipčević, J., & Cummins, P. R. (2015). Direct-seismogram inversion for receiver-side structure with uncertain source–time functions. *Geophys. J. Int.*, 203(2), 1373–1387. https://doi.org/10.1093/gji/ggv375
- Kennett, B. L. N. (2023). Interacting seismic waveguides: multimode surface waves and leaking modes. *Seismica*, 2(1). https://doi.org/10.26443/seismica.v2i1.282 (`DVSCON`)
- Kennett, B. L. N. (2026). Higher modes in the near surface: the influence of gradients. *Seismica*, 5(2). https://doi.org/10.26443/seismica.v5i2.2404 (`DVSCON`)

## Build

```
cd src && make     # gfortran/mpif90
```
LAPACK: Accelerate on macOS, `-lflexiblas` elsewhere; override with
`make LIB="-llapack -lblas"`. Flags in `src/Makefile.compiler`.

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
| `<base>_sdSWD.dat`, `<base>_sdSWD_M<m>.dat`, `<base>_sdSWDG*.dat` | `ICOV_SWD = 3` only: per-point sd (km/s), one file per slot, same rows as its data file. `sdparSWD` then becomes a dimensionless scale (bound it ~1, or `ISD_SWD = 0`). Missing/short file or sd <= 0 is fatal |
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

### Keyword lines (anywhere after line 48)

| keyword (example) | meaning; default |
|---|---|
| `MODE_OF 0 2` | Rayleigh mode number of each curve slot; files named by mode (`_SWD_M2.dat`). Default `0 1 ... NMODE-1` |
| `IGRP 0` | 0 phase, 1 group, for every slot unless `GRP_OF` is given |
| `GRP_OF 0 1` | velocity type per slot. `MODE_OF 0 0` + `GRP_OF 0 1` = joint R0 phase + group, each with its own grid, sigma and AR term; (mode, type) pairs must be unique. Phase and group of one mode on one grid share a root search |
| `SDMN_SWD 0.3 0.3`, `SDMX_SWD 3.0 3.0` | per-slot bounds of the hierarchical-sigma prior. Default: lines 33/34. The `map_voro` start must lie inside them |
| `DVSCON 0.100` | max \|adjacent-layer dVs\| (km/s), checked before the forward call (BayHunter lvz/hvz parity). Default off |
| `DVSMONO 0.050` | max allowed Vs *decrease* with depth (km/s); 0 = non-decreasing. Default off |
| `SWD_SCAN 0.08 1.6 0.005 0.001` | root scan `cmin cmax dc [dc_over]` (km/s). Default `2.0 6.5 0.05`, dc_over = dc/5; near-surface needs smaller values. Use `dc` 0.003 when Vp/Vs > 2 (see `docs/root_scan.md`) |
| `SWD_WARM 1` | warm-started root scan, 2-8x faster, bit-identical (see `docs/root_scan.md`). Default on iff `DVSCON > 0` or `DVSMONO >= 0` |
| `VP_BROCHER 1` | Vp from Vs, Brocher (2005) eq. 9; sampled Vp/Vs ignored. Default 0 |
| `RHO_BROCHER 1` | density from Vp, Brocher (2005) eq. 1. Default: legacy crustal 2.35 + 0.036 (Vp-3)^2 |
| `RHO_LINEAR 0.32 0.77` | density = a·Vp + b (g/cc); 0.32 0.77 = BayHunter. Overrides `RHO_BROCHER`. Default off |

Mode *n* is the *n*-th sign change of the DISPER80 secular function along the
c-scan; a model that cannot produce an observed mode at an observed period is
rejected. `ICOV_SWD = 2` and `Icov_iterUpdate_SWD = 1` hold one matrix for all
curves and are refused when `NMODE > 1`.

`tools/convert_legacy_inputs.py {rfswdellmt|receiver} SRC DST` converts the old
81-line RF-SWD-ELL-MT or 44/46-line receiver formats.

## Validation

- `example1_partial_coupling/SWD` IMAP logL 160.74617935123956 (macOS arm64,
  gfortran 16); x86 GCC 14: 160.74615570403239 (single-precision noise).
- HVC dam inputs: IMAP logL identical to the last digit to the receiver code
  (R0, R0+R1+R2 with `DVSCON`, `MODE_OF 0 2`), and a 160k-sample seeded run
  reproduces its posterior (medians to 1 m/s, same 68 % bands).
- Forward model vs disba: `validation/README.md`. Phase + group synthetic
  (Basel well model, `GRP_OF`): `validation/basel_group/README.md`.

## Post-processing (`tools/`, Python)

```bash
python tools/plot_rjhist.py RUN_DIR      # posteriors, convergence, data fits -> diagnostics.png
python tools/plot_datafit.py RUN_DIR     # observed vs IMAP prediction
python tools/print_map.py RUN_DIR        # max-logL sample -> <base>_map_voro.dat
python tools/plot_overlay.py A=RUN_A B=RUN_B --out fig.png
```
Library: `tools/swdell` (see its docstring). Legacy MATLAB scripts:
`plotting_scripts/legacy/`.
