"""Two targeted diagnostics for the cases that do not converge with the scan step.

(1) MODE IDENTITY.  For each period, find which disba mode index the Fortran's
    "mode m" curve actually corresponds to.  A clean m -> m match everywhere
    means the mode numbering agrees; a switch to m+1/m-1 at some period is a
    mode-counting divergence, not a precision problem.

(2) GROUP VELOCITY, independently.  disba's GroupDispersion finite-differences
    its own phase curve and is visibly unreliable in places (it returns U of
    10-30 m/s for modes whose phase velocity is 230-470 m/s, which is not
    physical).  So instead of trusting it, U is rebuilt here from the CONVERGED
    disba phase curve using the exact relation

        U = c / (1 + (T/c) dc/dT)

    with a central difference in T, and each estimate is only kept when the
    three phase values used are mutually consistent (no mode jump between
    them).  That is an independent check of DISPER80's analytic U.
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
NMODE_REF = 7


def fort_run(model, periods, maxmode, cmin, cmax, dc, dcov):
    fm, fp = os.path.join(TMP, "m.txt"), os.path.join(TMP, "p.txt")
    M.write_model_file(fm, model)
    np.savetxt(fp, periods, fmt="%.9f")
    r = subprocess.run(
        [DRIVER, fm, fp, str(maxmode), f"{cmin:.9f}", f"{cmax:.9f}",
         f"{dc:.9f}", f"{dcov:.9f}", "0"],
        capture_output=True, text=True, check=True)
    out = {m: {"c": np.full(len(periods), np.nan),
               "u": np.full(len(periods), np.nan)} for m in range(maxmode + 1)}
    for ln in r.stdout.strip().split("\n")[1:]:
        p = ln.split(",")
        m, i, ok = int(p[1]), int(p[2]) - 1, int(p[4])
        if ok:
            out[m]["c"][i] = float(p[5])
            out[m]["u"][i] = float(p[6])
    return out


def disba_curves(model, periods, dc=1e-5):
    from disba import PhaseDispersion
    th = model[:, 0].copy(); th[-1] = max(th[-1], 1.0)
    pd = PhaseDispersion(th, model[:, 2], model[:, 3], model[:, 1],
                         algorithm="dunkin", dc=dc)
    out = {}
    for m in range(NMODE_REF):
        v = np.full(len(periods), np.nan)
        try:
            cur = pd(periods, mode=m, wave="rayleigh")
            if len(cur.period):
                v[np.searchsorted(periods, cur.period)] = cur.velocity
        except Exception:
            pass
        out[m] = v
    return out, pd


def group_from_phase(pd, mode, periods, h=2e-3):
    """U = c / (1 + (T/c) dc/dT) from the converged disba phase curve."""
    U = np.full(len(periods), np.nan)
    for i, T in enumerate(periods):
        ts = np.array([T * (1 - h), T, T * (1 + h)])
        try:
            cur = pd(ts, mode=mode, wave="rayleigh")
        except Exception:
            continue
        if len(cur.period) != 3:
            continue
        cm, c0, cp = cur.velocity
        # reject if the three samples are not on one smooth branch
        if abs(cp - cm) > 0.05 * c0:
            continue
        dcdT = (cp - cm) / (ts[2] - ts[0])
        denom = 1.0 + (T / c0) * dcdT
        if abs(denom) > 1e-6:
            U[i] = c0 / denom
    return U


def main():
    mods = M.build_models()
    periods = np.logspace(np.log10(0.01), np.log10(0.6), 40)

    print("#" * 100)
    print("# (1) MODE IDENTITY: which disba mode does the Fortran's mode m match?")
    print("#" * 100)
    for name, mode in [("lvz_weak", 2), ("lid_over_channel", 1),
                       ("lid_over_channel", 2), ("lvz_weak", 1)]:
        model = mods[name]
        cmin, cmax = M.scan_window(model)
        f = fort_run(model, periods, 2, cmin, cmax, 0.0001, 0.0001)
        ref, _ = disba_curves(model, periods)
        print(f"\n{name}  Fortran mode {mode}   "
              f"(Vs = {', '.join(f'{v*1000:.0f}' for v in model[:,3])} m/s, "
              f"scan {cmin*1000:.0f}-{cmax*1000:.0f} m/s)")
        print(f"  {'T(s)':>8}{'f(Hz)':>8}{'c_fort':>11}   best disba match")
        prev = None
        for i, T in enumerate(periods):
            cf = f[mode]["c"][i]
            if np.isnan(cf):
                continue
            best, bd = None, 1e9
            for m in range(NMODE_REF):
                if not np.isnan(ref[m][i]) and abs(ref[m][i] - cf) < bd:
                    bd, best = abs(ref[m][i] - cf), m
            tag = ""
            if best != prev:
                tag = f"   <<< mode index {prev} -> {best}"
                prev = best
            nex = "none" if best is None else f"mode {best} (diff {bd*1000:8.3f} m/s)"
            if tag or i % 6 == 0:
                print(f"  {T:>8.4f}{1/T:>8.2f}{cf*1000:>11.3f}   {nex}{tag}")

    print()
    print("#" * 100)
    print("# (2) GROUP VELOCITY vs U rebuilt from the converged disba PHASE curve")
    print("#" * 100)
    hdr = f"{'model':<22}{'mode':>5}{'n':>5}{'u_max(m/s)':>13}{'u_rms':>10}{'u_rel':>11}"
    print(hdr); print("-" * len(hdr))
    worst = []
    for name, model in mods.items():
        cmin, cmax = M.scan_window(model)
        f = fort_run(model, periods, 2, cmin, cmax, 0.0001, 0.0001)
        _, pd = disba_curves(model, periods)
        for mode in range(3):
            uref = group_from_phase(pd, mode, periods)
            uf = f[mode]["u"]
            ok = ~np.isnan(uref) & ~np.isnan(uf)
            if ok.sum() == 0:
                print(f"{name:<22}{mode:>5}{0:>5}{'-':>13}{'-':>10}{'-':>11}")
                continue
            d = np.abs(uf[ok] - uref[ok]) * 1000
            rel = (d / (uref[ok] * 1000)).max()
            print(f"{name:<22}{mode:>5}{ok.sum():>5}{d.max():>13.4f}"
                  f"{np.sqrt((d**2).mean()):>10.4f}{rel:>11.2e}")
            worst.append((d.max(), rel, name, mode))
    worst.sort(reverse=True)
    print("\nworst group-velocity disagreements (Fortran vs U from disba phase):")
    for d, rel, n, m in worst[:5]:
        print(f"   {n:<22} mode {m}   {d:9.4f} m/s   ({rel:.2e} relative)")


if __name__ == "__main__":
    main()
