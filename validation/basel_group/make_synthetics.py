#!/usr/bin/env python
"""Synthetic Rayleigh dispersion curves from the Basel-1 / Otterbach-2 well
model (Michel et al. 2017) for the group-velocity inversion tests.

Writes data/
  synthetic_curves.csv   set, mode, grp, period, true, obs, sd  (grp 0 phase, 1 group)
                         set "engine": truth forwarded with the ENGINE's own
                         stack (legacy density, Vp = Vs*Vp/Vs, deep tail) --
                         the truth lies inside the model space;
                         set "michel": the literal Michel stack (true Vp,
                         Gardner density) -- the model-error case
  truth_nodes.csv        the truth as engine nodes (depth km, dVs, dVpVs)
  checks.txt             every cross-check, with pass/fail
  imap_truth/            the IMAP run that verifies the engine reproduces the
                         disp_driver curves slot by slot (phase AND group)
figures/fig_forward_check.png   DISPER80 vs disba (phase, group, group rebuilt
                         from disba's phase curve) and vs gpdc on this model

Run with the disba environment:
  /opt/anaconda3/envs/bayesbay_dev/bin/python make_synthetics.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

DATA = os.path.join(C.HERE, "data")
FIGS = os.path.join(C.HERE, "figures")
MAXMODE = 1
SLOTS = [(0, 0), (0, 1), (1, 0), (1, 1)]     # (mode, grp)


def log(fh, msg):
    print(msg)
    fh.write(msg + "\n")


def forward_truth(stack, periods):
    """dict (mode, grp) -> (valid, values) from disp_driver."""
    res = C.run_disp_driver([stack], periods, maxmode=MAXMODE)
    out = {}
    for mode in range(MAXMODE + 1):
        valid, c, u = res[(0, mode)]
        out[(mode, 0)] = (valid, c)
        out[(mode, 1)] = (valid, u)
    return out


def imap_identity_check(z, dvs, dvpvs, periods, truth, fh):
    """Run the engine (IMAP 1) with the truth as map_voro and the four slots on
    dummy data; its predicted rows must equal the disp_driver curves."""
    run = os.path.join(DATA, "imap_truth")
    os.makedirs(run, exist_ok=True)
    C.write_filebase(run)
    C.write_vel_ref(run)
    C.write_covparameter(run)
    ndat = len(periods) + 10
    with open(os.path.join(run, f"{C.BASE}_parameter.dat"), "w") as f:
        f.write(C.parameter_lines(4, [m for m, g in SLOTS], [g for m, g in SLOTS], ndat,
                                  imap=1, icov_swd=0, isd_swd=0, nptchains1=2))
    C.write_map_voro(os.path.join(run, f"{C.BASE}_map_voro.dat"), z, dvs, dvpvs, 4)
    for mode, grp in SLOTS:
        valid, _ = truth[(mode, grp)]
        C.write_curve(os.path.join(run, C.slot_filename(mode, grp)), periods[valid], np.ones(valid.sum()))
    logl, out = C.run_engine_imap(run, np_ranks=4)
    rows = C.read_mappred(run)
    worst = 0.0
    for i, (mode, grp) in enumerate(SLOTS):
        valid, val = truth[(mode, grp)]
        n = valid.sum()
        d = np.max(np.abs(rows[i, :n] - val[valid]))
        worst = max(worst, d)
        log(fh, f"  IMAP slot {i + 1} (mode {mode}, {'group' if grp else 'phase'}): "
                f"max |engine - disp_driver| = {d:.2e} km/s over {n} periods")
    reuse = sum("reuses the root solve" in l for l in out.splitlines())
    log(fh, f"  IMAP reported {reuse} slot(s) reusing an earlier root solve (expect 2)")
    # the engine carries the stack in single precision (curmod2 is REAL(SP) in
    # m/s, then /1000), disp_driver reads 7-digit km/s: ~1e-7 relative on the
    # layer velocities, i.e. ~1e-6 km/s on c and ~1e-4 km/s on the (more
    # sensitive) group velocity
    ok = worst < 5e-4 and reuse == 2
    log(fh, f"[{'PASS' if ok else 'FAIL'}] engine IMAP == disp_driver on the engine stack to single precision "
            f"(worst {worst * 1e3:.3f} m/s)")
    return ok


def own_fd_group_check(name, stack, periods, fh, h=2e-3, tol_ms=1.0):
    """U = c / (1 + (T/c) dc/dT) from OUR OWN phase curve (three periods per
    datum, T(1-h), T, T(1+h)) against the analytic U of DISPER80. This is the
    reference for the group velocity: disba's GroupDispersion and its mode
    indexing are not reliable on LVZ models (validation/README.md, fig 5)."""
    pfd = np.sort(np.concatenate([periods * (1 - h), periods, periods * (1 + h)]))
    res = C.run_disp_driver([stack], pfd, maxmode=MAXMODE)
    i0 = np.searchsorted(pfd, periods)
    im = np.searchsorted(pfd, periods * (1 - h))
    ip = np.searchsorted(pfd, periods * (1 + h))
    out = {}
    vs_hs = stack[-1, 3]
    for mode in range(MAXMODE + 1):
        v, c, u = res[(0, mode)]
        ok = v[i0] & v[im] & v[ip]
        # at a mode's cut-off (c within 0.1 % of the half-space Vs) the branch
        # bends too sharply for a 0.2 % finite difference; not a test of U there
        cutoff = ok & (c[i0] > 0.999 * vs_hs)
        with np.errstate(all="ignore"):
            dcdt = (c[ip] - c[im]) / (pfd[ip] - pfd[im])
            ufd = c[i0] / (1.0 + (periods / c[i0]) * dcdt)
        d = np.where(ok, (ufd - u[i0]) * 1e3, np.nan)
        worst = np.nanmax(np.abs(np.where(cutoff, np.nan, d)))
        out[mode] = d
        log(fh, f"  {name} R{mode}: analytic U vs U from our own phase curve (h = {h:g}): "
                f"max |diff| = {worst:.3f} m/s, median {np.nanmedian(np.abs(d)):.3f} m/s ({ok.sum()} periods"
                + (f"; {cutoff.sum()} at the cut-off excluded, |diff| {np.nanmax(np.abs(d[cutoff])):.1f} m/s there)" if cutoff.any() else ")"))
        log(fh, f"[{'PASS' if worst <= tol_ms else 'FAIL'}] {name} R{mode} group velocity self-consistent within {tol_ms:g} m/s")
    return out


def disba_curves(stack, periods, maxmode=3, dc=1e-5, dt=0.005):
    """disba phase / group per mode index, plus U rebuilt from disba's own
    phase curve (diag_modes.py::group_from_phase)."""
    from disba import GroupDispersion, PhaseDispersion
    th = stack[:, 0].copy()
    th[-1] = max(th[-1], 1.0)
    rho, vp, vs = stack[:, 1], stack[:, 2], stack[:, 3]
    pd = PhaseDispersion(th, vp, vs, rho, algorithm="dunkin", dc=dc)
    gd = GroupDispersion(th, vp, vs, rho, algorithm="dunkin", dc=dc, dt=dt)
    out = {}
    for mode in range(maxmode + 1):
        c = np.full(len(periods), np.nan)
        u = np.full(len(periods), np.nan)
        ufp = np.full(len(periods), np.nan)
        for solver, arr in ((pd, c), (gd, u)):
            try:
                cur = solver(periods, mode=mode, wave="rayleigh")
            except Exception:
                continue
            idx = np.searchsorted(periods, cur.period)
            ok = (idx < len(periods)) & np.isclose(periods[np.clip(idx, 0, len(periods) - 1)], cur.period, rtol=1e-9)
            arr[idx[ok]] = cur.velocity[ok]
        h = 2e-3
        for i, T in enumerate(periods):
            ts = np.array([T * (1 - h), T, T * (1 + h)])
            try:
                cur = pd(ts, mode=mode, wave="rayleigh")
            except Exception:
                continue
            if len(cur.period) != 3:
                continue
            cm, c0, cp = cur.velocity
            if abs(cp - cm) > 0.05 * c0:
                continue
            dcdT = (cp - cm) / (ts[2] - ts[0])
            den = 1.0 + (T / c0) * dcdT
            if abs(den) > 1e-6:
                ufp[i] = c0 / den
        out[mode] = {"c": c, "u": u, "u_from_c": ufp}
    return out


def disba_check(name, stack, periods, ours, fh, tol_c=0.02e-3, tol_u=0.01):
    """Match every DISPER80 root to the nearest disba phase root (any mode
    index: disba duplicates roots in LVZ models, which shifts its numbering)
    and compare c, U (disba GroupDispersion) and U rebuilt from disba's phase."""
    try:
        db = disba_curves(stack, periods)
    except ImportError:
        log(fh, "[SKIP] disba not importable in this python; run under /opt/anaconda3/envs/bayesbay_dev")
        return None, None
    rep = {}
    allc = np.array([db[m]["c"] for m in db])           # (nmodes, nper)
    for mode in range(MAXMODE + 1):
        valid, c, u = ours[(0, mode)]
        dc = np.full(len(periods), np.nan)
        du = np.full(len(periods), np.nan)
        dufc = np.full(len(periods), np.nan)
        for i in np.where(valid)[0]:
            col = allc[:, i]
            if np.all(np.isnan(col)):
                continue
            j = int(np.nanargmin(np.abs(col - c[i])))
            dc[i] = db[j]["c"][i] - c[i]
            du[i] = db[j]["u"][i] - u[i]
            dufc[i] = db[j]["u_from_c"][i] - u[i]
        rep[mode] = (dc, du, dufc)
        wc = np.nanmax(np.abs(dc)) * 1e3 if np.any(np.isfinite(dc)) else np.nan
        wu = np.nanmax(np.abs(du)) * 1e3 if np.any(np.isfinite(du)) else np.nan
        wufc = np.nanmax(np.abs(dufc)) * 1e3 if np.any(np.isfinite(dufc)) else np.nan
        log(fh, f"  {name} R{mode}: worst |c - disba c| = {wc:.3f} m/s, "
                f"|U - disba GroupDispersion| = {wu:.3f} m/s, "
                f"|U - U(disba phase curve)| = {wufc:.3f} m/s  ({int(np.isfinite(dc).sum())} periods)")
        okc = wc <= tol_c * 1e3
        log(fh, f"[{'PASS' if okc else 'FAIL'}] {name} R{mode} phase within {tol_c * 1e3:.2f} m/s of disba (roots matched by velocity)")
        if wufc > tol_u * 1e3 or wu > tol_u * 1e3:
            log(fh, f"[INFO] {name} R{mode}: disba's group velocity disagrees with ours at some periods; disba "
                    f"duplicates roots on this LVZ model (its mode indices shift), so its GroupDispersion and a "
                    f"finite difference along its mode index are not a reference here -- see the own-phase-curve check")
    return db, rep


def gpdc_check(stack, fh):
    """Compare our curves with gpdc's Michel2016_gpdc.disp (frequency, slowness s/m)."""
    if not os.path.exists(C.MICHEL_DISP):
        log(fh, "[SKIP] gpdc .disp file not found")
        return None
    rows = np.loadtxt(C.MICHEL_DISP, comments="#")
    T = 1.0 / rows[:, 0]
    v = 1.0 / rows[:, 1] / 1e3                       # km/s
    order = np.argsort(T)
    T, v = T[order], v[order]
    res = C.run_disp_driver([stack], T, maxmode=0)
    valid, c, u = res[(0, 0)]
    rc = np.abs(c[valid] / v[valid] - 1)
    ru = np.abs(u[valid] / v[valid] - 1)
    which = "GROUP" if np.median(ru) < np.median(rc) else "PHASE"
    r = ru if which == "GROUP" else rc
    log(fh, f"  gpdc .disp ({len(T)} points, {T.min():.2f}-{T.max():.2f} s) is the {which} slowness of R0: "
            f"median |rel diff| vs ours = {np.median(r) * 100:.3f} %, max = {np.max(r) * 100:.3f} % "
            f"(the other interpretation: median {np.median(rc if which == 'GROUP' else ru) * 100:.1f} %)")
    ok = np.max(r) < 0.005
    log(fh, f"[{'PASS' if ok else 'FAIL'}] our R0 {which.lower()} velocity within 0.5 % of gpdc on the literal Michel stack")
    return T, v, c, u, valid, which


def figure(stack_m, periods, ours_m, db, rep, gp, truth_e):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIGS, exist_ok=True)
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.8))
    # Vs profile
    z = np.r_[0.0, np.cumsum(stack_m[:-1, 0])]
    zz = np.repeat(np.r_[z, 6.0], 2)[1:-1]
    vv = np.repeat(stack_m[:, 3], 2)
    ax[0].plot(vv, zz, "k-", lw=1.5, label="Vs (Michel 2017)")
    ax[0].plot(np.repeat(stack_m[:, 2], 2), zz, "-", color="0.6", lw=1, label="Vp")
    ax[0].set_ylim(5.2, 0); ax[0].set_xlabel("velocity (km/s)"); ax[0].set_ylabel("depth (km)")
    ax[0].set_title("Basel-1 / Otterbach-2"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    cols = {0: "C0", 1: "C3"}
    for mode in range(MAXMODE + 1):
        valid, c, u = ours_m[(0, mode)]
        ax[1].plot(periods[valid], c[valid], "-", color=cols[mode], lw=2, label=f"R{mode} phase, DISPER80")
        ax[2].plot(periods[valid], u[valid], "-", color=cols[mode], lw=2, label=f"R{mode} group, DISPER80 (analytic U)")
        if db is not None:
            for j in db:
                lab = "disba phase (all mode indices)" if (mode == 0 and j == 0) else None
                ax[1].plot(periods, db[j]["c"], "o", ms=3, mfc="none", color="0.3", label=lab)
                lab = "disba GroupDispersion" if (mode == 0 and j == 0) else None
                ax[2].plot(periods, db[j]["u"], "s", ms=3, mfc="none", color="0.3", label=lab)
                lab = "U from disba phase curve" if (mode == 0 and j == 0) else None
                ax[2].plot(periods, db[j]["u_from_c"], "x", ms=4, color="C2", label=lab)
        if rep is not None and mode in rep:
            dc, du, dufc = rep[mode]
            ax[3].plot(periods, dc * 1e3, "o-", ms=3, color=cols[mode], label=f"R{mode}: disba c - ours")
            ax[3].plot(periods, dufc * 1e3, "x--", ms=4, color=cols[mode], label=f"R{mode}: U(disba phase) - our U")
    if gp is not None:
        T, v, c, u, valid, which = gp
        ax[2 if which == "GROUP" else 1].plot(T, v, ".", color="C1", ms=4, label=f"gpdc .disp ({which.lower()})")
    ax[1].set_title("phase velocity"); ax[2].set_title("group velocity")
    for a in ax[1:3]:
        a.set_xscale("log"); a.set_xlabel("period (s)"); a.set_ylabel("km/s"); a.grid(alpha=0.3); a.legend(fontsize=7)
    ax[3].set_xscale("log"); ax[3].set_xlabel("period (s)"); ax[3].set_ylabel("m/s"); ax[3].grid(alpha=0.3)
    ax[3].set_title("DISPER80 - disba (root matched by velocity)"); ax[3].legend(fontsize=7)
    fig.suptitle("Forward-model check on the Michel (2017) stack, 0.5-6 s", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "fig_forward_check.png"), dpi=150, bbox_inches="tight")
    print("wrote", os.path.join(FIGS, "fig_forward_check.png"))


def main():
    os.makedirs(DATA, exist_ok=True)
    periods = C.PERIODS_FULL
    stack_m = C.load_michel()
    z, dvs, dvpvs = C.truth_nodes(stack_m)
    stack_e = C.engine_stack(z, dvs, dvpvs)
    np.savetxt(os.path.join(DATA, "truth_nodes.csv"), np.column_stack([z, dvs, dvpvs]),
               header="depth_km,dVs_kms,dVpVs", delimiter=",", comments="", fmt="%.6f")
    np.savetxt(os.path.join(DATA, "engine_stack.csv"), stack_e, header="thick_km,rho_gcc,vp_kms,vs_kms",
               delimiter=",", comments="", fmt="%.6f")
    np.savetxt(os.path.join(DATA, "michel_stack.csv"), stack_m, header="thick_km,rho_gcc,vp_kms,vs_kms",
               delimiter=",", comments="", fmt="%.6f")

    with open(os.path.join(DATA, "checks.txt"), "w") as fh:
        log(fh, f"Michel stack: {len(stack_m)} layers, Vs {stack_m[:, 3].min():.3f}-{stack_m[:, 3].max():.3f} km/s, "
                f"interfaces down to {np.cumsum(stack_m[:-1, 0])[-1]:.3f} km")
        log(fh, f"engine stack: {len(stack_e)} layers (17 sampled + {len(C.TAIL_DEPTHS)} tail rows)")
        # --- noise-free truth curves, both stacks ---
        truth_e = forward_truth(stack_e, periods)
        truth_m = forward_truth(stack_m, periods)
        for mode in range(MAXMODE + 1):
            ve, _ = truth_e[(mode, 0)]
            vm, _ = truth_m[(mode, 0)]
            log(fh, f"  R{mode} exists at {ve.sum()}/{len(periods)} periods (engine stack), "
                    f"{vm.sum()}/{len(periods)} (Michel stack)")
        _, c0 = truth_e[(0, 0)]; _, u0 = truth_e[(0, 1)]
        log(fh, f"  R0 engine stack: c {c0[0]:.4f}-{c0[-1]:.4f} km/s, U {u0.min():.4f}-{u0.max():.4f} km/s "
                f"(Airy minimum at T = {periods[np.argmin(u0)]:.2f} s)")
        dm = max(np.max(np.abs(truth_m[(m, g)][1] - truth_e[(m, g)][1])[truth_e[(m, g)][0] & truth_m[(m, g)][0]])
                 for m, g in SLOTS)
        log(fh, f"  model-error size: max |Michel stack - engine stack| curve difference = {dm * 1e3:.1f} m/s "
                f"(density relation only)")
        # --- IMAP identity ---
        log(fh, "\n== engine IMAP vs disp_driver ==")
        imap_identity_check(z, dvs, dvpvs, periods, truth_e, fh)
        # --- disba ---
        log(fh, "\n== DISPER80 vs disba 0.7.0 (dunkin, dc 1e-5, dt 0.005) ==")
        res_m = C.run_disp_driver([stack_m], periods, maxmode=MAXMODE)
        db, rep = disba_check("Michel stack", stack_m, periods, res_m, fh)
        res_e = C.run_disp_driver([stack_e], periods, maxmode=MAXMODE)
        disba_check("engine stack", stack_e, periods, res_e, fh)
        log(fh, "\n== analytic group velocity vs finite difference of our own phase curve ==")
        own_fd_group_check("Michel stack", stack_m, periods, fh)
        own_fd_group_check("engine stack", stack_e, periods, fh)
        # --- gpdc ---
        log(fh, "\n== DISPER80 vs gpdc (Michel2016_gpdc.disp) ==")
        gp = gpdc_check(stack_m, fh)
        # --- noisy data ---
        log(fh, "\n== synthetic data ==")
        rng = np.random.default_rng(C.NOISE_SEED)
        zeta = {s: rng.standard_normal(len(periods)) for s in SLOTS}     # one realisation per slot
        rows = []
        for setname, truth in (("engine", truth_e), ("michel", truth_m)):
            for mode, grp in SLOTS:
                valid, val = truth[(mode, grp)]
                sd = C.NOISE[grp] * val
                obs = val + sd * zeta[(mode, grp)]
                for i in np.where(valid)[0]:
                    rows.append((setname, mode, grp, periods[i], val[i], obs[i], sd[i]))
        with open(os.path.join(DATA, "synthetic_curves.csv"), "w") as f:
            f.write("set,mode,grp,period,true,obs,sd\n")
            for r in rows:
                f.write("%s,%d,%d,%.9f,%.9f,%.9f,%.9f\n" % r)
        log(fh, f"  wrote {len(rows)} rows to synthetic_curves.csv: noise {C.NOISE[0] * 100:.0f} % phase, "
                f"{C.NOISE[1] * 100:.0f} % group, seed {C.NOISE_SEED}")
    figure(stack_m, periods, res_m, db, rep, gp, truth_e)


if __name__ == "__main__":
    main()
