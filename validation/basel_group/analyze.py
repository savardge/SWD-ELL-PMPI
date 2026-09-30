#!/usr/bin/env python
"""Diagnostics of the synthetic Basel-1 inversion runs.

For every run directory (runs/<name>/ holding <base>_voro_sample.txt):
  runs/<name>/diagnostics.png   the swdell rjhist figure (tools/swdell/plots.py)
                  with the truth overlays: posteriors vs the well profile,
                  convergence per chain, data fits with the noise-free curves
  results/summary.csv   one row per run: the general summary of the figure
                  (k, convergence metrics, sigma and chi2 per slot) plus the
                  synthetic-only recovery numbers (median error in 68 %
                  half-width units and coverage per depth band, true
                  interfaces found by a node-density heuristic)
  figures/overlay_<band>.png   Vs(z) median + band of the input combinations
                  of one period band over the truth; overlay_variants.png for
                  mismatch / model-error / icov1 against R0pg_full

  python analyze.py                       # every run in runs/ that has samples
  python analyze.py runs/R0pg_full ...    # selected runs
Options: --burn 0.3 (burn-in fraction per chain), --nfwd 300, --nprof 4000.
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
from swdell import io, model, plots  # noqa: E402

RUNS = os.path.join(C.HERE, "runs")
DATA = os.path.join(C.HERE, "data")
FIGS = os.path.join(C.HERE, "figures")
RESULTS = os.path.join(C.HERE, "results")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("runs", nargs="*")
    p.add_argument("--burn", type=float, default=0.3)
    p.add_argument("--nfwd", type=int, default=300)
    p.add_argument("--nprof", type=int, default=4000)
    p.add_argument("--no-overlays", action="store_true")
    return p.parse_args()


def truth_for(run):
    t = np.loadtxt(os.path.join(DATA, "truth_nodes.csv"), delimiter=",", skiprows=1)
    truth = {"label": "Michel 2017", "nodes": (t[:, 0], t[:, 1], t[:, 2])}
    p = os.path.join(run.run_dir, "truth_curves.csv")
    if os.path.exists(p):
        truth["curves"] = plots.truth_curves_from_csv(p, len(run.slots))
    p = os.path.join(run.run_dir, "truth_logL.txt")
    if os.path.exists(p):
        truth["logL"] = float(open(p).read())
    cfg = run.cfg
    truth["sigma"] = [1.0 if cfg["ICOV_SWD"] == 3 else C.NOISE[g] for m, g in cfg["slots"]]
    return truth


def recovery(ex, truth):
    """Synthetic-only numbers: median error in half-width units and coverage
    per depth band; true interfaces found (node-depth density within +-3 %
    of the depth exceeding twice its mean)."""
    zgrid, q, vs_true, zh, zbins = ex["zgrid"], ex["q"], ex["vs_true"], ex["zh"], ex["zbins"]
    out = {}
    for z0, z1 in ((0.0, 0.5), (0.5, 2.0), (2.0, 5.0)):
        m = (zgrid >= z0) & (zgrid < z1)
        half = np.maximum((q[2, m] - q[0, m]) / 2.0, 1e-4)
        tag = f"{z0:g}_{z1:g}km"
        out[f"err_{tag}"] = float(np.mean(np.abs(q[1, m] - vs_true[m]) / half))
        out[f"cover_{tag}"] = float(np.mean((vs_true[m] >= q[0, m]) & (vs_true[m] <= q[2, m])))
    zc = 0.5 * (zbins[1:] + zbins[:-1])
    z_true = truth["nodes"][0]
    found = 0
    for zt in z_true[1:]:
        w = np.abs(zc - zt) <= max(0.03 * zt, 0.03)
        if w.any() and zh[w].max() > 2.0 * zh.mean():
            found += 1
    out["interfaces_found"] = f"{found}/{len(z_true) - 1}"
    return out


def analyze_run(run_dir, args, rng):
    run = io.read_run(run_dir)
    truth = truth_for(run)
    print(f"{run.name}: {len(run.sample)} rows, NMODE {run.cfg['NMODE']}")
    fig, summary, ex = plots.rjhist_figure(run, truth=truth, burn=args.burn, nfwd=args.nfwd, nprof=args.nprof,
                                           rng=rng, title=f"{run.name}: synthetic Basel-1 test, slots "
                                           + ", ".join(io.slot_label(m, g, short=True) for m, g in run.cfg["slots"]))
    out = os.path.join(run_dir, "diagnostics.png")
    plots.save_close(fig, out)
    print("  wrote", out)
    rec = recovery(ex, truth)
    row = {k: summary[k] for k in ("run", "rows", "post", "chains", "k_median", "move_frac")}
    row.update(rec)
    row.update({k: summary[k] for k in ("stationarity_L1", "between_chain_L1", "burn30v60_L1", "logL_median", "truth_logL")})
    row.update({k: v for k, v in summary.items() if k.startswith("sigma_") or k.startswith("chi2_")})
    return row, ex


def overlays(results, zgrid, vs_true, vs_lim):
    groups = {}
    for band in ("full", "short", "long"):
        groups[f"overlay_{band}"] = [n for n in results if n.endswith(f"_{band}") and n.split("_")[0] in ("R0p", "R0g", "R0pg", "R01pg")]
    groups["overlay_variants"] = [n for n in ("R0pg_full", "R0pg_mismatch", "R0pg_modelerr", "R0pg_full_icov1") if n in results]
    for fname, names in groups.items():
        if len(names) < 2:
            continue
        fig = plots.overlay_figure([(n, results[n]) for n in names], zgrid, truth_vs=vs_true, vs_lim=vs_lim,
                                   truth_label="Michel 2017")
        plots.save_close(fig, os.path.join(FIGS, fname + ".png"))
        print("wrote", os.path.join(FIGS, fname + ".png"))


def main():
    args = parse_args()
    runs = args.runs or sorted(os.path.join(RUNS, d) for d in os.listdir(RUNS)
                               if os.path.exists(os.path.join(RUNS, d, f"{C.BASE}_voro_sample.txt")))
    if not runs:
        sys.exit("no run directory with a sample file")
    rng = np.random.default_rng(0)
    rows, results, first = [], {}, None
    for r in runs:
        try:
            row, ex = analyze_run(r, args, rng)
        except Exception as e:  # keep going through the matrix
            print(f"  {r}: FAILED: {e}")
            continue
        rows.append(row)
        results[row["run"]] = ex["q"]
        first = first or ex
    os.makedirs(RESULTS, exist_ok=True)
    if rows:
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        with open(os.path.join(RESULTS, "summary.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
        print("wrote", os.path.join(RESULTS, "summary.csv"))
    if not args.no_overlays and first is not None:
        overlays(results, first["zgrid"], first["vs_true"], (0.6, 3.6))


if __name__ == "__main__":
    main()
