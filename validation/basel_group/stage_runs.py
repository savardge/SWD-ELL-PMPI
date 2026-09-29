#!/usr/bin/env python
"""Stage the synthetic Basel-1 inversion matrix: one engine run directory per
cell under runs/<name>/, plus runs.list (one name per line, the job-array
index order of run_matrix.sbatch).

Matrix (user decisions 2026-09-29): input combinations R0p, R0g, R0p+R0g,
R0p+R0g+R1p+R1g (p = phase, g = group) x period bands full (0.5-6 s), short
(0.5-2 s), long (2-6 s); plus R0pg_mismatch (phase 0.5-6 s, group 1-4 s: no
shared grid, so no root-solve reuse), R0pg_modelerr (data from the literal
Michel stack: true Vp and Gardner density, outside the engine's model space)
and R0pg_full_icov1 (ICOV_SWD 1 + IMAGSCALE 1 with per-slot sigma bounds
instead of sd files).

Every run dir also gets truth_map_voro.dat (the Michel model as engine
nodes), truth_curves.csv (the noise-free curves of its slots) and
truth_logL.txt (the engine's logL of the truth on the noisy data, from an
IMAP run in a scratch copy), for analyze.py.

  python stage_runs.py            # needs data/synthetic_curves.csv (make_synthetics.py)
"""
import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

DATA = os.path.join(C.HERE, "data")
RUNS = os.path.join(C.HERE, "runs")

INPUTS = {                      # name -> slots (mode, grp)
    "R0p": [(0, 0)],
    "R0g": [(0, 1)],
    "R0pg": [(0, 0), (0, 1)],
    "R01pg": [(0, 0), (0, 1), (1, 0), (1, 1)],
}


def matrix():
    cells = []
    for band in ("full", "short", "long"):
        for name, slots in INPUTS.items():
            cells.append({"name": f"{name}_{band}", "slots": slots, "set": "engine",
                          "bands": {s: band for s in slots}, "icov": 3})
    cells.append({"name": "R0pg_mismatch", "slots": [(0, 0), (0, 1)], "set": "engine",
                  "bands": {(0, 0): "full", (0, 1): "mismatch_g"}, "icov": 3})
    cells.append({"name": "R0pg_modelerr", "slots": [(0, 0), (0, 1)], "set": "michel",
                  "bands": {(0, 0): "full", (0, 1): "full"}, "icov": 3})
    cells.append({"name": "R0pg_full_icov1", "slots": [(0, 0), (0, 1)], "set": "engine",
                  "bands": {(0, 0): "full", (0, 1): "full"}, "icov": 1})
    return cells


def load_curves():
    rows = np.genfromtxt(os.path.join(DATA, "synthetic_curves.csv"), delimiter=",", names=True,
                         dtype=None, encoding="utf-8")
    return rows


def slot_rows(curves, setname, mode, grp, band):
    lo, hi = C.BANDS[band]
    m = (curves["set"] == setname) & (curves["mode"] == mode) & (curves["grp"] == grp) \
        & (curves["period"] >= lo * (1 - 1e-9)) & (curves["period"] <= hi * (1 + 1e-9))
    r = curves[m]
    return r[np.argsort(r["period"])]


def stage(cell, curves, z, dvs, dvpvs):
    run = os.path.join(RUNS, cell["name"])
    if os.path.isdir(run):
        shutil.rmtree(run)
    os.makedirs(run)
    slots = cell["slots"]
    nmode = len(slots)
    C.write_filebase(run)
    C.write_vel_ref(run)
    C.write_covparameter(run)
    ndat = 0
    with open(os.path.join(run, "truth_curves.csv"), "w") as ft:
        ft.write("slot,mode,grp,period,true,obs,sd\n")
        for i, (mode, grp) in enumerate(slots):
            r = slot_rows(curves, cell["set"], mode, grp, cell["bands"][(mode, grp)])
            if len(r) == 0:
                raise RuntimeError(f"{cell['name']}: no data for slot {(mode, grp)}")
            ndat = max(ndat, len(r))
            C.write_curve(os.path.join(run, C.slot_filename(mode, grp)), r["period"], r["obs"])
            if cell["icov"] == 3:
                np.savetxt(os.path.join(run, C.slot_filename(mode, grp, sd=True)), r["sd"], fmt="%.9e")
            for row in r:
                ft.write("%d,%d,%d,%.9f,%.9f,%.9f,%.9f\n" % (i + 1, mode, grp, row["period"], row["true"], row["obs"], row["sd"]))
    ndat_swd = ndat + 5
    mode_of = [m for m, g in slots]
    grp_of = [g for m, g in slots]
    if cell["icov"] == 3:
        # sd files carry the noise; sdparSWD is a unit-scale multiplier
        kw = dict(icov_swd=3, imagscale=0, sdmn=(0.3, 1e-3), sdmx=(3.0, 1e-1))
        sd_start = 1.0
    else:
        # ICOV_SWD 1 + IMAGSCALE 1: sigma is a fraction of the datum; per-slot bounds
        kw = dict(icov_swd=1, imagscale=1, sdmn=(0.002, 1e-3), sdmx=(0.2, 1e-1),
                  sdmn_swd=[0.002 if g == 0 else 0.005 for g in grp_of],
                  sdmx_swd=[0.10 if g == 0 else 0.20 for g in grp_of])
        sd_start = [0.02 if g == 0 else 0.04 for g in grp_of]
    with open(os.path.join(run, f"{C.BASE}_parameter.dat"), "w") as f:
        f.write(C.parameter_lines(nmode, mode_of, grp_of, ndat_swd, **kw))
    # starting model: 3 nodes at the reference (Vs 2.1 km/s), sigma multipliers 1
    C.write_map_voro(os.path.join(run, f"{C.BASE}_map_voro.dat"), [0.0, 0.5, 2.0], [0.0, 0.0, 0.0],
                     [0.0, 0.0, 0.0], nmode, sdpar_swd=sd_start)
    C.write_map_voro(os.path.join(run, "truth_map_voro.dat"), z, dvs, dvpvs, nmode, sdpar_swd=sd_start)
    # logL of the truth on this run's noisy data (IMAP in a scratch copy)
    with tempfile.TemporaryDirectory() as tmp:
        C.copy_run_inputs(run, tmp)
        shutil.copy(os.path.join(run, "truth_map_voro.dat"), os.path.join(tmp, f"{C.BASE}_map_voro.dat"))
        with open(os.path.join(tmp, f"{C.BASE}_parameter.dat"), "w") as f:
            f.write(C.parameter_lines(nmode, mode_of, grp_of, ndat_swd, imap=1, nptchains1=2, **kw))
        logl, out = C.run_engine_imap(tmp, np_ranks=4)
        if "WARNING" in out:
            print("   ", [l for l in out.splitlines() if "WARNING" in l][0])
        # chi2 per slot from the IMAP residuals: obs - pred
        pred = C.read_mappred(tmp)
    with open(os.path.join(run, "truth_logL.txt"), "w") as f:
        f.write(f"{logl:.10f}\n")
    chi2 = []
    for i, (mode, grp) in enumerate(slots):
        r = slot_rows(curves, cell["set"], mode, grp, cell["bands"][(mode, grp)])
        res = (r["obs"] - pred[i, :len(r)]) / r["sd"]
        chi2.append(float(np.sum(res ** 2) / len(r)))
    print(f"  {cell['name']:<18} slots {slots} n={ndat} truth logL = {logl:.3f}  "
          f"reduced chi2 of the truth per slot = " + ", ".join(f"{c:.2f}" for c in chi2))
    return cell["name"]


def main():
    curves = load_curves()
    z, dvs, dvpvs = C.truth_nodes()
    os.makedirs(RUNS, exist_ok=True)
    names = [stage(cell, curves, z, dvs, dvpvs) for cell in matrix()]
    with open(os.path.join(RUNS, "runs.list"), "w") as f:
        f.write("\n".join(names) + "\n")
    print(f"staged {len(names)} runs in {RUNS}; job-array indices 0-{len(names) - 1} of runs.list")


if __name__ == "__main__":
    main()
