#!/usr/bin/env python
"""Diagnostics figure of one or more SWD-ELL-PMPI run directories
(posteriors / convergence / data fits; see tools/swdell/plots.py).

  python tools/plot_rjhist.py RUN_DIR [RUN_DIR ...]
      [--burn 0.3] [--nfwd 300] [--nprof 4000] [--zmax KM] [--out FILE]
      [--truth-map-voro FILE] [--truth-curves CSV] [--truth-logl X] [--truth-label S]

Writes <run_dir>/diagnostics.png (or --out for a single run) and prints the
summary line (k median, convergence metrics, sigma and chi2 per slot).
Posterior-predictive bands need validation/disp_driver (see validation/README.md);
without it the MAP prediction of an IMAP run is drawn when present.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from swdell import io, plots  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("runs", nargs="+")
    p.add_argument("--burn", type=float, default=0.3)
    p.add_argument("--nfwd", type=int, default=300)
    p.add_argument("--nprof", type=int, default=4000)
    p.add_argument("--zmax", type=float, default=None, help="depth axis limit (km); default hmx")
    p.add_argument("--out", default=None, help="output file (single run only)")
    p.add_argument("--truth-map-voro", default=None, help="true model in map_voro layout (synthetic tests)")
    p.add_argument("--truth-curves", default=None, help="noise-free curves CSV (slot,mode,grp,period,true,...)")
    p.add_argument("--truth-logl", type=float, default=None)
    p.add_argument("--truth-label", default="truth")
    a = p.parse_args()
    if a.out and len(a.runs) > 1:
        sys.exit("--out applies to a single run")
    for rd in a.runs:
        run = io.read_run(rd)
        truth = None
        if a.truth_map_voro or a.truth_curves or a.truth_logl is not None:
            truth = {"label": a.truth_label}
            if a.truth_map_voro:
                truth["nodes"] = plots.truth_from_map_voro(a.truth_map_voro, run.cfg)
            if a.truth_curves:
                truth["curves"] = plots.truth_curves_from_csv(a.truth_curves, len(run.slots))
            if a.truth_logl is not None:
                truth["logL"] = a.truth_logl
        out = a.out or os.path.join(rd, "diagnostics.png")
        fig, summary, _ = plots.rjhist_figure(run, truth=truth, burn=a.burn, nfwd=a.nfwd, nprof=a.nprof,
                                              zmax=a.zmax, rng=np.random.default_rng(0))
        plots.save_close(fig, out)
        print(f"{run.name}: rows {summary['rows']}, post {summary['post']}, k median {summary['k_median']:.0f}, "
              f"L1 stat/chain/burn {summary['stationarity_L1']:.3f}/{summary['between_chain_L1']:.3f}/{summary['burn30v60_L1']:.3f}; "
              + ", ".join(f"{k[6:]} sigma {v:.3g}" for k, v in summary.items() if k.startswith("sigma_") and np.isfinite(v))
              + f" -> {out}")


if __name__ == "__main__":
    main()
