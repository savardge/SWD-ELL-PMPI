"""Side-by-side list of every trapped Rayleigh root each code finds.

Also checks whether the apparent "missing modes" are just the top of the scan
window: modes pile up just below the half-space Vs as they approach cut-off, so
a cmax of 0.999*Vs_hs can cut off real roots.  The window is therefore pushed
to 0.999999*Vs_hs here.
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
NM = 10


def fort_roots(model, periods, cmin, cmax, dc):
    fm, fp = os.path.join(TMP, "mr.txt"), os.path.join(TMP, "pr.txt")
    M.write_model_file(fm, model)
    np.savetxt(fp, periods, fmt="%.9f")
    r = subprocess.run(
        [DRIVER, fm, fp, str(NM - 1), f"{cmin:.9f}", f"{cmax:.9f}",
         f"{dc:.9f}", f"{dc:.9f}", "0"],
        capture_output=True, text=True, check=True)
    c = np.full((NM, len(periods)), np.nan)
    for ln in r.stdout.strip().split("\n")[1:]:
        p = ln.split(",")
        if int(p[4]):
            c[int(p[1]), int(p[2]) - 1] = float(p[5])
    return c


def disba_roots(model, periods, dc=1e-5):
    from disba import PhaseDispersion
    th = model[:, 0].copy(); th[-1] = max(th[-1], 1.0)
    pd = PhaseDispersion(th, model[:, 2], model[:, 3], model[:, 1],
                         algorithm="dunkin", dc=dc)
    c = np.full((NM, len(periods)), np.nan)
    for m in range(NM):
        try:
            cur = pd(periods, mode=m, wave="rayleigh")
            if len(cur.period):
                c[m, np.searchsorted(periods, cur.period)] = cur.velocity
        except Exception:
            pass
    return c


mods = M.build_models()
periods = np.array([0.050, 0.124, 0.233, 0.438])

for name in ["lid_over_channel", "lvz_weak", "tailings_dry"]:
    model = mods[name]
    vs = model[:, 3]
    cmin = 0.80 * vs.min()
    print("=" * 104)
    print(f"{name}   Vs = {', '.join(f'{v*1000:.0f}' for v in vs)} m/s  "
          f"(half-space {vs[-1]*1000:.1f} m/s)")
    for pad, lbl in [(0.999, "cmax = 0.999*Vs_hs"),
                     (0.999999, "cmax = 0.999999*Vs_hs")]:
        cmax = pad * vs[-1]
        cf = fort_roots(model, periods, cmin, cmax, 0.0001)
        cd = disba_roots(model, periods)
        print(f"\n  {lbl}  ->  scan {cmin*1000:.1f} - {cmax*1000:.4f} m/s")
        for i, T in enumerate(periods):
            fl = [v * 1000 for v in cf[:, i] if not np.isnan(v)]
            dl = [v * 1000 for v in cd[:, i] if not np.isnan(v)]
            print(f"    T={T:.3f}s ({1/T:5.2f} Hz)  n_fort={len(fl)} n_disba={len(dl)}")
            print(f"        fortran: {' '.join(f'{v:8.2f}' for v in fl)}")
            print(f"        disba  : {' '.join(f'{v:8.2f}' for v in dl)}")
            # roots present in one list but not the other (1 m/s tolerance)
            miss_f = [v for v in dl if not any(abs(v - w) < 1.0 for w in fl)]
            extra_f = [v for v in fl if not any(abs(v - w) < 1.0 for w in dl)]
            if miss_f:
                print(f"        MISSED by fortran: {' '.join(f'{v:.2f}' for v in miss_f)}")
            if extra_f:
                print(f"        EXTRA in fortran : {' '.join(f'{v:.2f}' for v in extra_f)}")
    print()
