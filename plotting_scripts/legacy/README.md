# Legacy MATLAB post-processing (RF-SWD-ELL-MT format)

These scripts were written for the **original 81-line `<base>_parameter.dat`**
of RF-SWD-ELL-MT-PMPI (RF / MT lines, `sample.geom`, `NRF`, `NTIME`, ...) and
for sample rows that carried RF sigma / AR blocks. They do **not** read this
repository's files: `rf_read_parfile.m` takes lines 1-51 positionally (on a
48-line file every value after `IMAP` comes from the wrong line and the
keyword tail is never seen), the sample-column arithmetic assumes
`3*NRF` RF sigma columns, `data_fit.m` / `plot_postpred*.m` read a single
`_SWD.dat` and a one-row `_mappredSWD.dat`, and nothing knows the
`(MODE_OF, GRP_OF)` curve slots or the `_SWDG*` / `_sdSWD*` files.

Post-processing for the current formats is the Python package
`tools/swdell/` (`tools/plot_rjhist.py`, `plot_datafit.py`, `plot_overlay.py`,
`print_map.py`; see the README "Post-processing" section). Kept here for the
figure style and for runs of the old code line.

| script | role (old format) |
|---|---|
| `rf_plot_rjhist_varpar2.m` | main figure script (the most complete version; `varpar3` is a near-surface fork with a 53-line reader); reads `<base>_sample.mat` |
| `rf_read_parfile.m`, `rf_read_parfile3.m` | 81-line parameter file + `sample.geom` |
| `misc/convert_sample.m`, `convert_samples.m` | `_voro_sample.txt` -> `.mat` with burn-in (the plural version overwrites its input file) |
| `rf_convert_sample_laynode_to_lay.m`, `rf_laynode_to_lay.m`, `rf_voro_to_lay.m`, `rf_getref.m`, `plprof.m` | nodes -> layers, reference interpolation (`rf_getref` matches the engine's `GETREF`) |
| `rf_print_map.m`, `rf_print_map2.m`, `rf_abs_to_refrel.m` | MAP / true model -> `_map_voro.dat` |
| `data_fit.m`, `plot_postpred*.m`, `make_sample_postpred.m` | data fit and posterior-predictive plots (the `postpred` program no longer exists) |
| `make_input_files_cdi_*.m`, `rf_est_covmat.m`, `plot_cov_mats.m` | empirical data-covariance estimation (`ICOV = 2`) |
| `plot_postlog.m`, `plot_interface_prob*.m`, `batch_*.m`, `rf_plot_datafit.m`, `rf_plot_raysum.m`, `rf_plot_trace_azimuth.m` | dataset comparison / RF figures |

Known problems if you run them on old-format data: `rf_plot_rjhist_varpar2.m`
uses `XTick_Vs` / `XTick_VpVs` / `XTick_res` that are commented out
(l. 1481-1483), its chain panels test an undefined `NTH`, and lines starting
with `!` are shell escapes. They need `misc/get_loc.m`, `misc/hpd.m` and
`cloudPlot/` on the path (`addpath(genpath('plotting_scripts'))`).
