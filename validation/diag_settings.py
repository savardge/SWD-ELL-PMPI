"""(A) How close can cmax get to the half-space Vs before the propagator fails?
(B) Accuracy vs cost as a function of the c-scan step, for choosing SWD_SCAN.
"""
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models as M  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DRIVER = os.path.join(HERE, "disp_driver")
TMP = os.path.join(HERE, "out", "diag")
os.makedirs(TMP, exist_ok=True)


def run(model, periods, maxmode, cmin, cmax, dc, dcov, iwarm=1):
    fm, fp = os.path.join(TMP, "ms.txt"), os.path.join(TMP, "ps.txt")
    M.write_model_file(fm, model)
    np.savetxt(fp, periods, fmt="%.9f")
    t0 = time.time()
    r = subprocess.run(
        [DRIVER, fm, fp, str(maxmode), f"{cmin:.9f}", f"{cmax:.9f}",
         f"{dc:.9f}", f"{dcov:.9f}", str(iwarm)],
        capture_output=True, text=True, check=True)
    dt = time.time() - t0
    c = np.full((maxmode + 1, len(periods)), np.nan)
    nvalid = 0
    for ln in r.stdout.strip().split("\n")[1:]:
        p = ln.split(",")
        if int(p[4]):
            c[int(p[1]), int(p[2]) - 1] = float(p[5])
            nvalid += 1
    return c, nvalid, dt


def disba_ref(model, periods, nm, dc=1e-5):
    from disba import PhaseDispersion
    th = model[:, 0].copy(); th[-1] = max(th[-1], 1.0)
    pd = PhaseDispersion(th, model[:, 2], model[:, 3], model[:, 1],
                         algorithm="dunkin", dc=dc)
    c = np.full((nm, len(periods)), np.nan)
    for m in range(nm):
        try:
            cur = pd(periods, mode=m, wave="rayleigh")
            if len(cur.period):
                c[m, np.searchsorted(periods, cur.period)] = cur.velocity
        except Exception:
            pass
    return c


mods = M.build_models()
periods = np.logspace(np.log10(0.01), np.log10(0.6), 40)

print("#" * 92)
print("# (A) cmax vs the half-space Vs: where does the propagator give up?")
print("#" * 92)
print(f"{'model':<20}{'Vs_hs(m/s)':>12}  " + "".join(f"{p:>11}" for p in
      ["0.990", "0.999", "0.9999", "0.99999", "0.999999"]))
print("-" * 92)
for name in ["lvz_weak", "tailings_dry", "grad_smooth"]:
    model = mods[name]
    vshs = model[:, 3][-1]
    cmin = 0.8 * model[:, 3].min()
    cells = []
    for pad in [0.990, 0.999, 0.9999, 0.99999, 0.999999]:
        _, nv, _ = run(model, periods, 2, cmin, pad * vshs, 0.0001, 0.0001)
        cells.append(f"{nv:>11d}")
    print(f"{name:<20}{vshs*1000:>12.1f}  " + "".join(cells))
print("(numbers are valid (period,mode) predictions out of 120; a collapse means"
      "\n the propagator returned a hard error and the whole call was aborted)")

print()
print("#" * 92)
print("# (B) accuracy vs cost of the c-scan step (modes 0-2, 40 periods,")
print("#     scan window fitted to each model, warm start on)")
print("#" * 92)
hdr = (f"{'dc = dc_over [m/s]':<22}{'worst |dc| (m/s)':>18}{'median |dc|':>14}"
       f"{'cpu (ms/model)':>17}")
print(hdr); print("-" * len(hdr))
usable = [n for n in mods if n not in ("2layer_normal", "near_homogeneous")]
for dcv in [0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005]:
    worst, alld, tt = 0.0, [], 0.0
    for name in usable:
        model = mods[name]
        cmin, cmax = M.scan_window(model)
        c, _, dt = run(model, periods, 2, cmin, cmax, dcv, dcv)
        tt += dt
        ref = disba_ref(model, periods, 3)
        # compare only where the two agree on the mode branch (within 1 m/s at
        # the finest step), so mode-index divergence does not pollute the
        # precision measurement
        cfine, _, _ = run(model, periods, 2, cmin, cmax, 0.00002, 0.00002)
        ok = (~np.isnan(c) & ~np.isnan(ref) & ~np.isnan(cfine)
              & (np.abs(cfine - ref) * 1000 < 1.0))
        if ok.sum():
            d = np.abs(c[ok] - ref[ok]) * 1000
            worst = max(worst, d.max()); alld.extend(d.tolist())
    print(f"{dcv*1000:<22.3f}{worst:>18.4f}{np.median(alld):>14.4f}"
          f"{tt/len(usable)*1000:>17.1f}")
