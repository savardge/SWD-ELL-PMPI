"""Post-processing of SWD-ELL-PMPI run directories (current file formats).

Readers mirror src/read_input.f90 (READPARFILE, READDATA, SWD_SLOT_FILE) and
the sample/map layouts of src/prjmh_temper_rf.f90; the model reconstruction
mirrors GETREF / INTERPLAYER / MAKE_CURMOD and the deep tail of LOGLHOOD_SWD
(src/loglhood.f90). Entry points:

    from swdell import io, model, forward, posterior, plots
    run = io.read_run("path/to/run_dir")            # inputs + sample matrix
    plots.rjhist_figure(run, out="diagnostics.png")   # posterior / convergence / fits

Command-line front ends: tools/plot_rjhist.py, tools/plot_datafit.py,
tools/plot_overlay.py, tools/print_map.py.
"""
__version__ = "1.0"
