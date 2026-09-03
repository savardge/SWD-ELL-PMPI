"""Verify the DISPER80 (raydsp/raydspn) Rayleigh forward model against disba.

For every test model in models.py, compares phase and group velocity of the
fundamental mode and the first two overtones, period by period, between:

  this code : dispersion() -> RAYDSPN -> RAYMRX   (Saito DISPER-80, single
              precision, mode selected by counting secular-function sign
              changes on the c-scan grid)
  disba     : PhaseDispersion / GroupDispersion   (Dunkin matrix, float64,
              independent root finder)

Both are given exactly the same layer stack, so any disagreement is numerical
or algorithmic, not a difference in the model.

Usage:  python compare_disba.py [--dc DC] [--dc-over DCOVER] [--outdir DIR]
"""

import argparse
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models as M  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DRIVER = os.path.join(HERE, "disp_driver")
MAXMODE = 2


def run_driver(model, periods, cmin, cmax, dc, dc_over, iwarm, tmpdir):
    """Run the Fortran driver on one model; return dict[mode] -> (ivalid, c, u)."""
    fmod = os.path.join(tmpdir, "model.txt")
    fper = os.path.join(tmpdir, "periods.txt")
    M.write_model_file(fmod, model)
    np.savetxt(fper, periods, fmt="%.9f")

    cmd = [
        DRIVER, fmod, fper, str(MAXMODE),
        f"{cmin:.9f}", f"{cmax:.9f}", f"{dc:.9f}", f"{dc_over:.9f}", str(iwarm),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"driver failed: {res.stderr[:2000]}")

    lines = res.stdout.strip().split("\n")[1:]  # drop header
    out = {m: {"ivalid": np.zeros(len(periods), int),
               "c": np.zeros(len(periods)),
               "u": np.zeros(len(periods))} for m in range(MAXMODE + 1)}
    for ln in lines:
        _, mode, iper, _, iv, c, u = ln.split(",")
        mode, iper = int(mode), int(iper) - 1
        out[mode]["ivalid"][iper] = int(iv)
        out[mode]["c"][iper] = float(c)
        out[mode]["u"][iper] = float(u)
    return out


def run_disba(model, periods, dc=1e-5, dt=0.005):
    """disba reference; returns dict[mode] -> (mask, c, u) aligned to periods."""
    from disba import GroupDispersion, PhaseDispersion

    # disba wants columns thickness, vp, vs, rho
    thickness = model[:, 0].copy()
    rho, vp, vs = model[:, 1], model[:, 2], model[:, 3]
    # give the half-space a nominal thickness; both codes ignore it
    thickness[-1] = max(thickness[-1], 1.0)

    pd = PhaseDispersion(thickness, vp, vs, rho, algorithm="dunkin", dc=dc)
    gd = GroupDispersion(thickness, vp, vs, rho, algorithm="dunkin", dc=dc, dt=dt)

    # disba needs periods ascending and returns only the periods where the
    # requested mode exists
    out = {}
    for mode in range(MAXMODE + 1):
        mask = np.zeros(len(periods), bool)
        c = np.zeros(len(periods))
        u = np.zeros(len(periods))
        for solver, arr in ((pd, c), (gd, u)):
            try:
                cur = solver(periods, mode=mode, wave="rayleigh")
            except Exception:
                continue
            if len(cur.period) == 0:
                continue
            idx = np.searchsorted(periods, cur.period)
            ok = (idx < len(periods)) & np.isclose(
                periods[np.clip(idx, 0, len(periods) - 1)], cur.period, rtol=1e-9
            )
            arr[idx[ok]] = cur.velocity[ok]
            if solver is pd:
                mask[idx[ok]] = True
        out[mode] = {"mask": mask, "c": c, "u": u}
    return out


def compare(name, model, periods, dc, dc_over, tmpdir, iwarm=1, disba_dt=0.005,
            disba_dc=1e-5):
    cmin, cmax = M.scan_window(model)
    fort = run_driver(model, periods, cmin, cmax, dc, dc_over, iwarm, tmpdir)
    # The reference must NOT depend on the setting under test: disba has its own
    # root-search step and is only converged below ~1e-4 km/s (at 5e-4 it skips
    # modes outright), so it is pinned at a converged value here.
    ref = run_disba(model, periods, dc=disba_dc, dt=disba_dt)

    rows = []
    for mode in range(MAXMODE + 1):
        f, r = fort[mode], ref[mode]
        both = (f["ivalid"] == 1) & r["mask"]
        only_f = (f["ivalid"] == 1) & ~r["mask"]
        only_r = (f["ivalid"] == 0) & r["mask"]
        # group velocity is only meaningful where disba actually produced one
        both_u = both & (r["u"] > 0)
        dc_abs = np.abs(f["c"][both] - r["c"][both]) * 1000.0  # m/s
        du_abs = np.abs(f["u"][both_u] - r["u"][both_u]) * 1000.0
        rows.append(
            dict(
                model=name, mode=mode, n_both=int(both.sum()),
                n_only_fortran=int(only_f.sum()), n_only_disba=int(only_r.sum()),
                c_max=float(dc_abs.max()) if dc_abs.size else np.nan,
                c_rms=float(np.sqrt((dc_abs**2).mean())) if dc_abs.size else np.nan,
                c_relmax=float((dc_abs / (r["c"][both] * 1000.0)).max())
                if dc_abs.size else np.nan,
                u_max=float(du_abs.max()) if du_abs.size else np.nan,
                u_rms=float(np.sqrt((du_abs**2).mean())) if du_abs.size else np.nan,
                u_relmax=float((du_abs / (r["u"][both_u] * 1000.0)).max())
                if du_abs.size else np.nan,
            )
        )
    return rows, fort, ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dc", type=float, default=0.001)
    ap.add_argument("--dc-over", type=float, default=0.0002)
    ap.add_argument("--disba-dt", type=float, default=0.005)
    ap.add_argument("--disba-dc", type=float, default=1e-5)
    ap.add_argument("--iwarm", type=int, default=1)
    ap.add_argument("--outdir", default=os.path.join(HERE, "out"))
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tmpdir = os.path.join(args.outdir, "tmp")
    os.makedirs(tmpdir, exist_ok=True)

    # 1.67 - 100 Hz: the band a shallow MASW/DAS survey actually resolves
    periods = np.logspace(np.log10(0.01), np.log10(0.6), 40)

    mods = M.build_models()
    all_rows = []
    detail = []
    for name, model in mods.items():
        rows, fort, ref = compare(
            name, model, periods, args.dc, args.dc_over, tmpdir,
            iwarm=args.iwarm, disba_dt=args.disba_dt,
            disba_dc=args.disba_dc,
        )
        all_rows.extend(rows)
        for mode in range(MAXMODE + 1):
            for i, T in enumerate(periods):
                detail.append(
                    (name, mode, T, fort[mode]["ivalid"][i], int(ref[mode]["mask"][i]),
                     fort[mode]["c"][i], ref[mode]["c"][i],
                     fort[mode]["u"][i], ref[mode]["u"][i])
                )

    # ---- summary table ----
    hdr = (f"{'model':<20}{'mode':>5}{'n':>5}{'onlyF':>7}{'onlyD':>7}"
           f"{'c_max(m/s)':>12}{'c_rms':>10}{'c_rel':>10}"
           f"{'u_max(m/s)':>12}{'u_rms':>10}{'u_rel':>10}")
    print(f"\nFortran scan dc={args.dc} dc_over={args.dc_over} km/s, warm={args.iwarm}"
          f"  |  disba reference dc={args.disba_dc} dt={args.disba_dt}")
    print(hdr)
    print("-" * len(hdr))
    for r in all_rows:
        print(f"{r['model']:<20}{r['mode']:>5}{r['n_both']:>5}"
              f"{r['n_only_fortran']:>7}{r['n_only_disba']:>7}"
              f"{r['c_max']:>12.4f}{r['c_rms']:>10.4f}{r['c_relmax']:>10.2e}"
              f"{r['u_max']:>12.4f}{r['u_rms']:>10.4f}{r['u_relmax']:>10.2e}")

    fin = [r for r in all_rows if np.isfinite(r["c_max"])]
    if fin:
        print(f"\nworst phase |diff| : {max(r['c_max'] for r in fin):.4f} m/s "
              f"({max(r['c_relmax'] for r in fin):.2e} relative)")
    finu = [r for r in all_rows if np.isfinite(r["u_max"])]
    if finu:
        print(f"worst group |diff| : {max(r['u_max'] for r in finu):.4f} m/s "
              f"({max(r['u_relmax'] for r in finu):.2e} relative)")

    np.savetxt(
        os.path.join(args.outdir, "detail.csv"),
        np.array([(d[0], d[1], f"{d[2]:.9f}", d[3], d[4],
                   f"{d[5]:.9f}", f"{d[6]:.9f}", f"{d[7]:.9f}", f"{d[8]:.9f}")
                  for d in detail], dtype=object),
        fmt="%s", delimiter=",",
        header="model,mode,period,ivalid_fortran,valid_disba,c_fortran,c_disba,u_fortran,u_disba",
        comments="",
    )
    print(f"\nper-period detail -> {os.path.join(args.outdir, 'detail.csv')}")


if __name__ == "__main__":
    main()
