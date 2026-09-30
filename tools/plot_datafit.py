#!/usr/bin/env python
"""Observed curves vs the MAP prediction of an IMAP run, per slot
(the replacement of plotting_scripts/legacy/data_fit.m).

  python tools/plot_datafit.py RUN_DIR [--out FILE]

RUN_DIR must hold <base>_mappredSWD.dat, i.e. the engine was run with IMAP 1
on its <base>_map_voro.dat (tools/print_map.py writes one from a sample file).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from swdell import io, plots  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("run_dir")
    p.add_argument("--out", default=None)
    a = p.parse_args()
    run = io.read_run(a.run_dir, with_sample=False)
    out = a.out or os.path.join(a.run_dir, "datafit.png")
    fig = plots.data_fit_figure(run)
    plots.save_close(fig, out)
    print("wrote", out)


if __name__ == "__main__":
    main()
