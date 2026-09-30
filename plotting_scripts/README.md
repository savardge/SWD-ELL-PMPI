# plotting_scripts

The group's shared MATLAB toolbox, most of which is unrelated to this
repository (`ffi_*`, `cmt_*`, `auv_*`, `refl_*`, oases, tsunami, beach-ball,
... scripts of other projects). It is kept as is.

The scripts that post-processed **this** rjMcMC code line were written for
the original RF-SWD-ELL-MT file formats and have been moved to `legacy/`
(see `legacy/README.md`). They cannot read the current 48-line parameter
file, the `(MODE_OF, GRP_OF)` curve slots or the current sample layout.

Post-processing for the current formats is in Python:

```bash
python tools/plot_rjhist.py RUN_DIR          # posteriors / convergence / data fits
python tools/plot_datafit.py RUN_DIR         # observed vs the MAP prediction of an IMAP run
python tools/print_map.py RUN_DIR            # max-logL sample -> <base>_map_voro.dat
python tools/plot_overlay.py A=RUN_A B=RUN_B --out fig.png
```

(`tools/swdell/`: readers mirroring `src/read_input.f90`, model
reconstruction mirroring `GETREF` / `INTERPLAYER` / `MAKE_CURMOD`, forward
through `validation/disp_driver`, figures.)
