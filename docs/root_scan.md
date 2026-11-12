# DISPER80 root scan: warm start and step size

Detail behind the `SWD_WARM` and `SWD_SCAN` keywords (see the main README).

## Warm-started root scan (speed)

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

## Root-scan step and high Vp/Vs

The fundamental root is found by a sign-change
scan in steps of `dc` (`SWD_SCAN`). For soft, high-Vp/Vs stacks (saturated soil,
Vp/Vs 4-8) a 5 m/s step (`dc` 0.005) can step over the root pair and the model is
rejected as if no mode existed: one 175/300/900 m/s test model was rejected at
`dc` 0.005 and fine at 0.002, with the same logL to 4 digits at 0.002 and 0.001.
On 60 posterior-like structures under Brocher Vp/density, 59 were valid at 0.005
and all 60 at 0.003, so use `SWD_SCAN 0.08 2.5 0.003 0.001` (1.7x the cost) for
any run that leaves the Vp/Vs = 2 regime. Validation of the Brocher path: the
layer stack the engine prints matches an independent implementation to four
digits (Vp to 0.0001 m/s, density exactly) and the predicted dispersion matches
disba on that stack to 0.00 m/s; with both keywords off the predicted data are
byte-identical to the pre-Brocher build (`ac335c6`).

## Validation

Validation (`swd/test_warm_driver.f90`, `tools/export_test_models.py`): warm
and cold curves compared on 15 000 models drawn from the real HVC posteriors
(v2 R0+R1+R2, v3 R0..R3, and the NLMX = 30 / hmin = 2 m run), all three modes
-- **zero differing periods, zero differing validity patterns**. On 4 000
models of the unconstrained fast-lid posterior the fundamental and first
overtone are still exact and the second overtone differs at 0.2 % of periods,
which is the regime the default switch excludes.
