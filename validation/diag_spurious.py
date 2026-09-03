"""(A) Are the extra Fortran modes real?

At a given period, count how many trapped Rayleigh modes each code finds.  If
the Fortran reports a mode where disba finds no such root, the sign-change
count in RAYDSPN has picked up something that is not a zero of the dispersion
function.

(B) Is DISPER80's analytic group velocity self-consistent?

U is compared against a high-accuracy numerical derivative of DISPER80's OWN
phase-velocity curve, U = c / (1 + (T/c) dc/dT).  This uses no external code,
so it separates "our U is wrong" from "the disba reference is wrong".
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


def fort_run(model, periods, maxmode, cmin, cmax, dc, dcov, tag="a"):
    fm, fp = os.path.join(TMP, f"m{tag}.txt"), os.path.join(TMP, f"p{tag}.txt")
    M.write_model_file(fm, model)
    np.savetxt(fp, periods, fmt="%.9f")
    r = subprocess.run(
        [DRIVER, fm, fp, str(maxmode), f"{cmin:.9f}", f"{cmax:.9f}",
         f"{dc:.9f}", f"{dcov:.9f}", "0"],
        capture_output=True, text=True, check=True)
    c = np.full((maxmode + 1, len(periods)), np.nan)
    u = np.full((maxmode + 1, len(periods)), np.nan)
    for ln in r.stdout.strip().split("\n")[1:]:
        p = ln.split(",")
        m, i, ok = int(p[1]), int(p[2]) - 1, int(p[4])
        if ok:
            c[m, i] = float(p[5])
            u[m, i] = float(p[6])
    return c, u


def disba_count(model, periods, nmax=8, dc=1e-5):
    from disba import PhaseDispersion
    th = model[:, 0].copy(); th[-1] = max(th[-1], 1.0)
    pd = PhaseDispersion(th, model[:, 2], model[:, 3], model[:, 1],
                         algorithm="dunkin", dc=dc)
    c = np.full((nmax, len(periods)), np.nan)
    for m in range(nmax):
        try:
            cur = pd(periods, mode=m, wave="rayleigh")
            if len(cur.period):
                c[m, np.searchsorted(periods, cur.period)] = cur.velocity
        except Exception:
            pass
    return c


def part_a():
    print("#" * 100)
    print("# (A) number of trapped modes found by each code")
    print("#" * 100)
    mods = M.build_models()
    periods = np.array([0.010, 0.050, 0.124, 0.233, 0.438, 0.600])
    for name in ["lid_over_channel", "lvz_weak", "lvz_strong", "tailings_dry"]:
        model = mods[name]
        cmin, cmax = M.scan_window(model)
        NM = 7
        cf, _ = fort_run(model, periods, NM - 1, cmin, cmax, 0.0001, 0.0001)
        cd = disba_count(model, periods, nmax=NM)
        print(f"\n{name}   Vs = {', '.join(f'{v*1000:.0f}' for v in model[:,3])} m/s"
              f"   scan {cmin*1000:.0f}-{cmax*1000:.0f} m/s")
        print(f"  {'T(s)':>8}{'f(Hz)':>8}{'n_fortran':>11}{'n_disba':>9}   "
              f"fortran c [m/s]")
        for i, T in enumerate(periods):
            nf = int(np.sum(~np.isnan(cf[:, i])))
            nd = int(np.sum(~np.isnan(cd[:, i])))
            vals = " ".join(f"{v*1000:7.1f}" for v in cf[:, i] if not np.isnan(v))
            flag = "   <-- extra mode(s) in Fortran" if nf > nd else ""
            print(f"  {T:>8.4f}{1/T:>8.2f}{nf:>11}{nd:>9}   {vals}{flag}")


def part_b():
    print()
    print("#" * 100)
    print("# (B) DISPER80 analytic U vs derivative of DISPER80's OWN phase curve")
    print("#    (self-consistency; no external reference involved)")
    print("#" * 100)
    mods = M.build_models()
    periods = np.logspace(np.log10(0.01), np.log10(0.6), 40)
    h = 2e-3
    hdr = (f"{'model':<22}{'mode':>5}{'n':>5}{'u_max(m/s)':>13}{'u_rms':>10}"
           f"{'u_rel':>11}   worst at")
    print(hdr); print("-" * (len(hdr) + 10))
    worst = []
    for name, model in mods.items():
        cmin, cmax = M.scan_window(model)
        # phase curve at T, T(1-h), T(1+h) -- one driver call on the merged grid
        grid = np.sort(np.concatenate([periods * (1 - h), periods, periods * (1 + h)]))
        c_all, u_all = fort_run(model, grid, 2, cmin, cmax, 0.0001, 0.0001, tag="b")
        idx = {round(float(t), 12): j for j, t in enumerate(grid)}
        for mode in range(3):
            un, ua, tw = [], [], []
            for T in periods:
                jm = idx[round(float(T * (1 - h)), 12)]
                j0 = idx[round(float(T), 12)]
                jp = idx[round(float(T * (1 + h)), 12)]
                cm, c0, cp = c_all[mode, jm], c_all[mode, j0], c_all[mode, jp]
                if np.isnan(cm) or np.isnan(c0) or np.isnan(cp):
                    continue
                if abs(cp - cm) > 0.05 * c0:      # branch jump between samples
                    continue
                dcdT = (cp - cm) / (T * 2 * h)
                den = 1.0 + (T / c0) * dcdT
                if abs(den) < 1e-6 or np.isnan(u_all[mode, j0]):
                    continue
                un.append(c0 / den); ua.append(u_all[mode, j0]); tw.append(T)
            if not un:
                print(f"{name:<22}{mode:>5}{0:>5}{'-':>13}{'-':>10}{'-':>11}")
                continue
            un, ua, tw = np.array(un), np.array(ua), np.array(tw)
            d = np.abs(ua - un) * 1000
            k = int(np.argmax(d))
            rel = (d / (un * 1000)).max()
            print(f"{name:<22}{mode:>5}{len(un):>5}{d.max():>13.4f}"
                  f"{np.sqrt((d**2).mean()):>10.4f}{rel:>11.2e}"
                  f"   T={tw[k]:.4f}s ({1/tw[k]:.1f} Hz)")
            worst.append((d.max(), rel, name, mode, tw[k]))
    worst.sort(reverse=True)
    print("\nworst self-consistency failures:")
    for d, rel, n, m, T in worst[:6]:
        print(f"   {n:<22} mode {m}  {d:9.4f} m/s ({rel:.2e})  at T={T:.4f}s")


if __name__ == "__main__":
    part_a()
    part_b()
