#!/usr/bin/env python
"""Vs(z) posterior median + 16-84 % band of several runs side by side.

  python tools/plot_overlay.py LABEL=RUN_DIR [LABEL=RUN_DIR ...] --out FILE
      [--burn 0.3] [--nprof 4000] [--zmax KM] [--truth-map-voro FILE]

The truth file (map_voro layout) is drawn with the first run's reference model.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from swdell import io, model, plots, posterior  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("entries", nargs="+", metavar="LABEL=RUN_DIR")
    p.add_argument("--out", required=True)
    p.add_argument("--burn", type=float, default=0.3)
    p.add_argument("--nprof", type=int, default=4000)
    p.add_argument("--zmax", type=float, default=None)
    p.add_argument("--truth-map-voro", default=None)
    a = p.parse_args()
    rng = np.random.default_rng(0)
    entries, zgrid, truth_vs, vlim = [], None, None, None
    for e in a.entries:
        label, rd = e.split("=", 1) if "=" in e else (os.path.basename(os.path.normpath(e)), e)
        run = io.read_run(rd)
        if zgrid is None:
            zgrid = posterior.depth_grid(run, zmax=a.zmax)
            vlim = plots.vs_limits(run, zgrid)
            if a.truth_map_voro:
                z, dvs, _ = plots.truth_from_map_voro(a.truth_map_voro, run.cfg)
                truth_vs = model.vs_on_grid(z, dvs, zgrid, run.cfg, run.vel_ref)
        post, _ = io.burn_in_split(run.sample, a.burn)
        vs = posterior.profiles(post, run, zgrid, a.nprof, rng)[0]
        entries.append((label, posterior.quantiles(vs)))
    fig = plots.overlay_figure(entries, zgrid, truth_vs=truth_vs, vs_lim=vlim)
    plots.save_close(fig, a.out)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
