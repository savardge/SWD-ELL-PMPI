"""Why do LVZ fundamentals disagree at high frequency?

Sweeps the c-scan step (dc) and shifts the scan origin (cmin) off the branch
point c = min(Vs), for the LVZ models, and checks whether the Fortran root
converges to disba's.
"""
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models as M  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DRIVER = os.path.join(HERE, "disp_driver")
TMP = os.path.join(HERE, "out", "diag")
os.makedirs(TMP, exist_ok=True)


def fort(model, periods, mode, cmin, cmax, dc):
    fm = os.path.join(TMP, "m.txt")
    fp = os.path.join(TMP, "p.txt")
    M.write_model_file(fm, model)
    np.savetxt(fp, periods, fmt="%.9f")
    r = subprocess.run(
        [DRIVER, fm, fp, str(mode), f"{cmin:.9f}", f"{cmax:.9f}",
         f"{dc:.9f}", f"{dc:.9f}", "0"],
        capture_output=True, text=True, check=True)
    c = np.zeros(len(periods))
    for ln in r.stdout.strip().split("\n")[1:]:
        p = ln.split(",")
        if int(p[1]) == mode:
            c[int(p[2]) - 1] = float(p[5]) if int(p[4]) else np.nan
    return c


def disba_c(model, periods, mode, dc):
    from disba import PhaseDispersion
    th = model[:, 0].copy(); th[-1] = max(th[-1], 1.0)
    pd = PhaseDispersion(th, model[:, 2], model[:, 3], model[:, 1],
                         algorithm="dunkin", dc=dc)
    out = np.full(len(periods), np.nan)
    cur = pd(periods, mode=mode, wave="rayleigh")
    idx = np.searchsorted(periods, cur.period)
    out[idx] = cur.velocity
    return out


mods = M.build_models()
periods = np.array([0.010, 0.0123, 0.0152, 0.020, 0.030])

for name in ["lvz_strong", "lvz_deep", "lid_over_channel"]:
    model = mods[name]
    vsmin = model[:, 3].min()
    cmin0, cmax = M.scan_window(model)
    print("=" * 100)
    print(f"{name}   min(Vs) = {vsmin*1000:.1f} m/s   "
          f"scan window {cmin0*1000:.1f} - {cmax*1000:.1f} m/s")
    # is the branch point c = min(Vs) exactly on the scan grid?
    for dc in [0.001, 0.0002]:
        k = (vsmin - cmin0) / dc
        print(f"   dc={dc*1000:6.3f} m/s : (min(Vs)-cmin)/dc = {k:.6f} "
              f"{'<-- branch point ON the grid' if abs(k-round(k)) < 1e-6 else ''}")

    print(f"\n{'':>12}" + "".join(f"{1/T:>12.1f}Hz" for T in periods))
    ref = disba_c(model, periods, 0, 0.0001)
    print(f"{'disba':>12}" + "".join(f"{v*1000:>14.4f}" for v in ref))
    for dc in [0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002]:
        c = fort(model, periods, 0, cmin0, cmax, dc)
        print(f"{'dc=' + f'{dc*1000:.3f}':>12}"
              + "".join(f"{v*1000:>14.4f}" for v in c))
    # shift the scan origin so the branch point is NOT on the grid
    for shift in [0.5, 0.37]:
        for dc in [0.001, 0.0002]:
            cmin_s = cmin0 + shift * dc
            c = fort(model, periods, 0, cmin_s, cmax, dc)
            print(f"{'+' + f'{shift}dc,{dc*1000:.2f}':>12}"
                  + "".join(f"{v*1000:>14.4f}" for v in c))
    print()
