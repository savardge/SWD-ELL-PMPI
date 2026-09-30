#!/usr/bin/env python
"""Write the maximum-logL sample of a run as <base>_map_voro.dat (the
replacement of plotting_scripts/legacy/rf_print_map.m), so the engine can
predict its data with IMAP 1 (-> <base>_mappredSWD.dat, tools/plot_datafit.py).

  python tools/print_map.py RUN_DIR [--k K] [--burn 0.3] [--out FILE]

--k restricts the search to samples with K nodes. The previous map file is
kept as <base>_map_voro.dat.bak.
"""
import argparse
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from swdell import io, posterior  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("run_dir")
    p.add_argument("--k", type=int, default=None)
    p.add_argument("--burn", type=float, default=0.3)
    p.add_argument("--out", default=None)
    a = p.parse_args()
    run = io.read_run(a.run_dir)
    post, _ = io.burn_in_split(run.sample, a.burn)
    ks = post[:, 3].astype(int)
    for k in np.unique(ks):
        sub = post[ks == k]
        print(f"k = {k:2d}: {len(sub):7d} samples, max logL {sub[:, 0].max():.4f}")
    row = posterior.map_row(post, a.k)
    out = a.out or run.path("_map_voro.dat")
    if os.path.exists(out):
        shutil.copy(out, out + ".bak")
    np.savetxt(out, io.sample_row_to_map(row, run.cfg)[None, :], fmt="%.10e")
    print(f"wrote {out}: k = {int(row[3])}, logL = {row[0]:.6f}")


if __name__ == "__main__":
    main()
