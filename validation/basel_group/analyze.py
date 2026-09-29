#!/usr/bin/env python
"""Diagnostics of the synthetic Basel-1 inversion runs.

For every run directory (runs/<name>/ holding <base>_voro_sample.txt):

  runs/<name>/diagnostics.png   three blocks
    posteriors    Vs(z) per-depth histogram with the 16/50/84 % profile and the
                  TRUE well profile; interface probability with the true
                  interfaces; Vp/Vs(z); p(k) with its Poisson prior; the
                  hierarchical sigma of every slot
    convergence   logL and k per cold chain vs sample index (burn-in marked,
                  truth logL as a reference); running median Vs at three
                  depths per chain; the three L1 stationarity metrics of
                  masw-das/rjmcmc_convergence.py (first vs second half,
                  worst chain pair, 30 % vs 60 % burn-in; <= 0.05 good,
                  >= 0.15 re-run), acceptance, sample counts
    data fits     per slot: observed curve with sd bars, the noise-free truth
                  and the posterior-predictive 16/50/84 % band from
                  forwarding NFWD post-burn-in models with disp_driver (phase
                  AND group from one root search); residuals in sigma units
                  with the reduced chi-square

  results/summary.csv   one row per run (recovery, sigma, chi2, convergence)
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

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402

RUNS = os.path.join(C.HERE, "runs")
DATA = os.path.join(C.HERE, "data")
FIGS = os.path.join(C.HERE, "figures")
RESULTS = os.path.join(C.HERE, "results")
ZGRID = np.arange(0.0, C.HMX + 1e-9, 0.01)       # km
ZEDGES = np.r_[ZGRID - 0.005, ZGRID[-1] + 0.005]   # pcolormesh cell edges
VS_BINS = np.linspace(0.6, 3.6, 151)
VPVS_BINS = np.linspace(1.65, 2.05, 81)
Z_BINS = np.arange(0.0, C.HMX + 1e-9, 0.02)
PROBE_DEPTHS = (0.2, 1.0, 3.0)                   # km, running-median traces


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("runs", nargs="*")
    p.add_argument("--burn", type=float, default=0.3)
    p.add_argument("--nfwd", type=int, default=300)
    p.add_argument("--nprof", type=int, default=4000)
    p.add_argument("--no-overlays", action="store_true")
    return p.parse_args()


# --------------------------------------------------------------------------
def truth_profile():
    t = np.loadtxt(os.path.join(DATA, "truth_nodes.csv"), delimiter=",", skiprows=1)
    z, dvs, dvpvs = t[:, 0], t[:, 1], t[:, 2]
    return z, C.vs_on_grid(z, dvs, ZGRID), C.vpvs_on_grid(z, dvpvs, ZGRID)


def profiles(rows, nlmx, rng, nmax):
    if nmax and len(rows) > nmax:
        rows = rows[rng.choice(len(rows), nmax, replace=False)]
    vs = np.empty((len(rows), len(ZGRID)))
    vpvs = np.empty_like(vs)
    zall = []
    for i, r in enumerate(rows):
        z, dvs, dvpvs = C.sample_nodes(r, nlmx)
        vs[i] = C.vs_on_grid(z, dvs, ZGRID)
        vpvs[i] = C.vpvs_on_grid(z, dvpvs, ZGRID)
        zall.append(z[z > 0])
    return vs, vpvs, np.concatenate(zall) if zall else np.zeros(0)


def quantiles(vs):
    return np.percentile(vs, [16, 50, 84], axis=0)


def l1(qa, qb, ref):
    """Depth-mean |median difference| in units of ref's 16-84 half-width."""
    half = np.maximum((ref[2] - ref[0]) / 2.0, 1e-4)
    return float(np.mean(np.abs(qa[1] - qb[1]) / half))


def probe_vs(rows, nlmx):
    out = np.empty((len(rows), len(PROBE_DEPTHS)))
    zp = np.array(PROBE_DEPTHS)
    for i, r in enumerate(rows):
        z, dvs, _ = C.sample_nodes(r, nlmx)
        out[i] = C.vs_on_grid(z, dvs, zp)
    return out


def running_median(x, npts=80):
    idx = np.unique(np.linspace(1, len(x), npts).astype(int))
    return idx, np.array([np.median(x[:i], axis=0) for i in idx])


def forward_posterior(rows, nlmx, cfg, slot_data, rng, nfwd):
    """Posterior-predictive curves per slot: array (nfwd, npts) with NaN
    where the mode does not exist. One disp_driver call per distinct grid."""
    if len(rows) > nfwd:
        rows = rows[rng.choice(len(rows), nfwd, replace=False)]
    stacks = [C.engine_stack(*C.sample_nodes(r, nlmx)) for r in rows]
    grids = {}
    for i, s in enumerate(slot_data):
        key = tuple(np.round(s["periods"], 9))
        grids.setdefault(key, []).append(i)
    pred = [np.full((len(rows), len(s["periods"])), np.nan) for s in slot_data]
    for key, islots in grids.items():
        periods = np.array(key)
        maxmode = max(cfg["slots"][i][0] for i in islots)
        res = C.run_disp_driver(stacks, periods, maxmode=maxmode)
        for i in islots:
            mode, grp = cfg["slots"][i]
            for im in range(len(stacks)):
                valid, c, u = res[(im, mode)]
                v = u if grp == 1 else c
                pred[i][im] = np.where(valid, v, np.nan)
    return pred


# --------------------------------------------------------------------------
def analyze_run(run_dir, args, truth, rng):
    name = os.path.basename(os.path.normpath(run_dir))
    cfg = C.read_cfg(run_dir)
    lay = C.sample_layout(cfg)
    slot_data = C.read_slot_data(run_dir, cfg)
    dat = C.read_sample(run_dir)
    if dat.shape[1] != lay["ncol"]:
        raise RuntimeError(f"{name}: sample file has {dat.shape[1]} columns, layout expects {lay['ncol']}")
    nlmx, nmode = cfg["NLMX"], cfg["NMODE"]
    isig = lay["isig"]
    truth_logl = float(open(os.path.join(run_dir, "truth_logL.txt")).read()) \
        if os.path.exists(os.path.join(run_dir, "truth_logL.txt")) else np.nan
    tcurves = np.genfromtxt(os.path.join(run_dir, "truth_curves.csv"), delimiter=",", names=True)
    z_true, vs_true, vpvs_true = truth

    post, burn = C.burn_in_split(dat, args.burn)
    chains = np.unique(dat[:, -1]).astype(int)
    print(f"{name}: {len(dat)} rows, {len(chains)} chains, {len(post)} post-burn-in samples, NMODE {nmode}")

    # ---- posteriors ----
    vs, vpvs, znodes = profiles(post, nlmx, rng, args.nprof)
    q = quantiles(vs)
    vs_h = np.array([np.histogram(vs[:, i], VS_BINS)[0] for i in range(len(ZGRID))], float)
    vs_h /= np.maximum(vs_h.max(axis=1, keepdims=True), 1)
    vpvs_h = np.array([np.histogram(vpvs[:, i], VPVS_BINS)[0] for i in range(len(ZGRID))], float)
    vpvs_h /= np.maximum(vpvs_h.max(axis=1, keepdims=True), 1)
    zh, _ = np.histogram(znodes, Z_BINS)
    zh = zh / max(zh.sum() * (Z_BINS[1] - Z_BINS[0]), 1)
    kpost = post[:, 3].astype(int)
    kk = np.arange(C.NLMN, nlmx + 1)
    from math import lgamma
    lam = cfg["lambda"]
    pk = np.exp(-lam + kk * np.log(lam) - np.array([lgamma(k + 1) for k in kk]))
    pk /= pk.sum()
    sig = post[:, isig:isig + nmode]

    # ---- convergence ----
    q_first = quantiles(profiles(post[: len(post) // 2], nlmx, rng, args.nprof)[0])
    q_second = quantiles(profiles(post[len(post) // 2:], nlmx, rng, args.nprof)[0])
    m_stat = l1(q_first, q_second, q)
    qc = [quantiles(profiles(post[post[:, -1] == c], nlmx, rng, args.nprof // 2)[0]) for c in chains]
    m_chain = max((l1(qc[a], qc[b], q) for a in range(len(qc)) for b in range(a + 1, len(qc))), default=0.0)
    post60, _ = C.burn_in_split(dat, 0.6)
    q60 = quantiles(profiles(post60, nlmx, rng, args.nprof)[0])
    q30 = q if abs(args.burn - 0.3) < 1e-9 else quantiles(profiles(C.burn_in_split(dat, 0.3)[0], nlmx, rng, args.nprof)[0])
    m_burn = l1(q30, q60, q)
    # the engine's acceptance column is 0/0 on the master rank; use the fraction
    # of consecutive samples of a chain that differ (a move happened) instead
    moves = []
    for c in chains:
        sub = post[post[:, -1] == c]
        if len(sub) > 1:
            moves.append(np.mean(np.any(sub[1:, 3:isig] != sub[:-1, 3:isig], axis=1)))
    acc = np.array(moves if moves else [np.nan])

    # ---- data fits ----
    pred = forward_posterior(post, nlmx, cfg, slot_data, rng, args.nfwd)
    fit_q, chi2, exist = [], [], []
    for i, s in enumerate(slot_data):
        with np.errstate(all="ignore"):
            fq = np.nanpercentile(pred[i], [16, 50, 84], axis=0)
        fit_q.append(fq)
        exist.append(np.mean(np.isfinite(pred[i]), axis=0))
        sd = s["sd"] if s["sd"] is not None else np.abs(s["obs"]) * np.median(sig[:, i])
        res = (s["obs"] - fq[1]) / sd
        chi2.append(float(np.nanmean(res ** 2)))

    # ---- recovery numbers ----
    def band_err(z0, z1):
        m = (ZGRID >= z0) & (ZGRID < z1)
        half = np.maximum((q[2, m] - q[0, m]) / 2.0, 1e-4)
        return float(np.mean(np.abs(q[1, m] - vs_true[m]) / half)), \
            float(np.mean((vs_true[m] >= q[0, m]) & (vs_true[m] <= q[2, m])))
    rec = {b: band_err(*b) for b in ((0.0, 0.5), (0.5, 2.0), (2.0, 5.0))}
    # interface recovery: a true interface counts as found when the node-depth
    # density within +-3 % of its depth exceeds twice the mean density
    zc = 0.5 * (Z_BINS[1:] + Z_BINS[:-1])
    found = 0
    for zt in z_true[1:]:
        w = np.abs(zc - zt) <= max(0.03 * zt, 0.03)
        if w.any() and zh[w].max() > 2.0 * zh.mean():
            found += 1
    n_int = len(z_true) - 1

    summary = {"run": name, "rows": len(dat), "post": len(post), "chains": len(chains),
               "k_median": float(np.median(kpost)), "move_frac": float(np.mean(acc)),
               "err_0_0.5km": rec[(0.0, 0.5)][0], "cover_0_0.5km": rec[(0.0, 0.5)][1],
               "err_0.5_2km": rec[(0.5, 2.0)][0], "cover_0.5_2km": rec[(0.5, 2.0)][1],
               "err_2_5km": rec[(2.0, 5.0)][0], "cover_2_5km": rec[(2.0, 5.0)][1],
               "interfaces_found": f"{found}/{n_int}",
               "stationarity_L1": m_stat, "between_chain_L1": m_chain, "burn30v60_L1": m_burn,
               "logL_median": float(np.median(post[:, 0])), "truth_logL": truth_logl}
    for i, (mode, grp) in enumerate(cfg["slots"]):
        tag = f"R{mode}{'g' if grp else 'p'}"
        summary[f"sigma_{tag}"] = float(np.median(sig[:, i]))
        summary[f"chi2_{tag}"] = chi2[i]

    # ================= figure =================
    ncol = max(5, nmode)
    fig = plt.figure(figsize=(3.6 * ncol, 15))
    gs = gridspec.GridSpec(4, ncol, figure=fig, height_ratios=[1.6, 1.0, 1.0, 0.7], hspace=0.42, wspace=0.32)
    cmap = plt.get_cmap("bone_r").copy(); cmap.set_under("white")

    ax = fig.add_subplot(gs[0, 0])
    ax.pcolormesh(VS_BINS, ZEDGES, vs_h, cmap=cmap, vmin=1e-3, vmax=1, shading="flat")
    ax.plot(q[1], ZGRID, "w--", lw=1, label="posterior median")
    ax.plot(q[0], ZGRID, "-", color="0.6", lw=0.6); ax.plot(q[2], ZGRID, "-", color="0.6", lw=0.6)
    ax.plot(vs_true, ZGRID, "r-", lw=1.4, label="true (Michel 2017)")
    ax.set_ylim(C.HMX, 0); ax.set_xlim(VS_BINS[0], VS_BINS[-1]); ax.set_xlabel("Vs (km/s)"); ax.set_ylabel("depth (km)")
    ax.set_title("Vs posterior"); ax.legend(fontsize=7, loc="lower left")

    ax = fig.add_subplot(gs[0, 1])
    ax.fill_betweenx(zc, 0, zh, color="0.5", alpha=0.7)
    for zt in z_true[1:]:
        ax.axhline(zt, color="r", lw=0.6, alpha=0.7)
    ax.set_ylim(C.HMX, 0); ax.set_xlabel("node-depth density (1/km)"); ax.set_title(f"interfaces ({found}/{n_int} true found)")

    ax = fig.add_subplot(gs[0, 2])
    ax.pcolormesh(VPVS_BINS, ZEDGES, vpvs_h, cmap=cmap, vmin=1e-3, vmax=1, shading="flat")
    ax.plot(vpvs_true, ZGRID, "r-", lw=1.4)
    ax.set_ylim(C.HMX, 0); ax.set_xlabel("Vp/Vs"); ax.set_title("Vp/Vs posterior")

    ax = fig.add_subplot(gs[0, 3])
    kh = np.array([(kpost == k).mean() for k in kk])
    ax.bar(kk, kh, color="0.4", label="posterior")
    ax.plot(kk, pk, "o-", color="C1", ms=3, label=f"Poisson({lam:g}) prior")
    ax.axvline(len(z_true), color="r", lw=1, label=f"truth k = {len(z_true)}")
    ax.set_xlabel("number of nodes k"); ax.set_title("p(k)"); ax.legend(fontsize=7)

    ax = fig.add_subplot(gs[0, 4])
    for i, (mode, grp) in enumerate(cfg["slots"]):
        tag = f"R{mode} {'group' if grp else 'phase'}"
        if cfg["ICOV_SWD"] == 3:
            x, lo, hi, ref = sig[:, i], cfg["sdmn"][i], cfg["sdmx"][i], 1.0
            xl = "sigma multiplier (sd files carry the noise)"
        else:
            x, lo, hi, ref = 100 * sig[:, i], 100 * cfg["sdmn"][i], 100 * cfg["sdmx"][i], 100 * C.NOISE[grp]
            xl = "sigma (% of datum)"
        h, e = np.histogram(x, bins=np.linspace(lo, min(hi, np.percentile(x, 99.5) * 1.3), 60), density=True)
        ax.step(e[:-1], h, where="post", color=f"C{i}", label=f"{tag}: median {np.median(x):.2f}")
        ax.axvline(ref, color=f"C{i}", ls=":", lw=1)
    ax.set_xlabel(xl); ax.set_title("hierarchical sigma per slot (dotted = expected)"); ax.legend(fontsize=7)

    # ---- row 2: convergence ----
    ax = fig.add_subplot(gs[1, 0:2])
    for j, c in enumerate(chains):
        sub = dat[dat[:, -1] == c]
        ax.plot(np.arange(len(sub)), sub[:, 0], lw=0.5, alpha=0.8, color=plt.get_cmap("tab10")(j % 10), label=f"chain (rank {c})" if j < 6 else None)
        ax.axvline(int(len(sub) * args.burn), color=plt.get_cmap("tab10")(j % 10), lw=0.5, ls="--")
    if np.isfinite(truth_logl):
        ax.axhline(truth_logl, color="r", lw=1, label=f"truth logL = {truth_logl:.1f}")
    lo = np.percentile(post[:, 0], 0.5)
    ax.set_ylim(lo - 0.3 * abs(lo - post[:, 0].max()) - 5, post[:, 0].max() + 5)
    ax.set_xlabel("sample index within chain"); ax.set_ylabel("logL"); ax.set_title("logL per cold chain (dashed = burn-in cut)")
    ax.legend(fontsize=7, ncol=2)

    ax = fig.add_subplot(gs[1, 2])
    for j, c in enumerate(chains):
        sub = dat[dat[:, -1] == c]
        ax.plot(np.arange(len(sub)), sub[:, 3], lw=0.5, alpha=0.8, color=plt.get_cmap("tab10")(j % 10))
    ax.axhline(len(z_true), color="r", lw=1)
    ax.set_xlabel("sample index within chain"); ax.set_ylabel("k"); ax.set_title("k per chain")

    ax = fig.add_subplot(gs[1, 3])
    for j, c in enumerate(chains):
        sub = dat[dat[:, -1] == c]
        pv = probe_vs(sub[:: max(1, len(sub) // 3000)], nlmx)
        idx, rm = running_median(pv)
        for d in range(len(PROBE_DEPTHS)):
            ax.plot(idx * max(1, len(sub) // 3000), rm[:, d], lw=0.8, color=plt.get_cmap("tab10")(j % 10), ls=["-", "--", ":"][d])
    for d, zp in enumerate(PROBE_DEPTHS):
        ax.axhline(np.interp(zp, ZGRID, vs_true), color="r", lw=0.8, ls=["-", "--", ":"][d], label=f"truth at {zp:g} km")
    ax.set_xlabel("sample index within chain"); ax.set_ylabel("running median Vs (km/s)")
    ax.set_title("running median Vs per chain"); ax.legend(fontsize=7)

    ax = fig.add_subplot(gs[1, 4]); ax.axis("off")
    txt = (f"{name}\n\nrows {len(dat)}, chains {len(chains)}\npost-burn-in ({args.burn:.0%}) {len(post)}\n\n"
           f"stationarity (1st vs 2nd half)  {m_stat:.3f}\nbetween-chain (worst pair)      {m_chain:.3f}\n"
           f"burn-in 30 % vs 60 %            {m_burn:.3f}\n(L1 in 68 % half-width units;\n <= 0.05 good, >= 0.15 re-run)\n\n"
           f"moved-sample fraction  {np.mean(acc):.3f}\nk median {np.median(kpost):.0f} (truth {len(z_true)})\n"
           f"logL median {np.median(post[:, 0]):.1f}  truth {truth_logl:.1f}\n")
    for i, (mode, grp) in enumerate(cfg["slots"]):
        txt += f"R{mode}{'g' if grp else 'p'}: sigma {np.median(sig[:, i]):.3f}, chi2_red {chi2[i]:.2f}\n"
    ax.text(0.0, 1.0, txt, va="top", ha="left", family="monospace", fontsize=8.5, transform=ax.transAxes)

    # ---- rows 3-4: data fits ----
    for i, s in enumerate(slot_data):
        mode, grp = cfg["slots"][i]
        tag = f"R{mode} {'group' if grp else 'phase'} velocity"
        tr = tcurves[tcurves["slot"] == i + 1]
        ax = fig.add_subplot(gs[2, i])
        ax.fill_between(s["periods"], fit_q[i][0], fit_q[i][2], color="C0", alpha=0.3, label="posterior predictive 16-84 %")
        ax.plot(s["periods"], fit_q[i][1], "-", color="C0", lw=1.2, label="posterior predictive median")
        ax.plot(tr["period"], tr["true"], "r-", lw=1, label="noise-free truth")
        sd = s["sd"] if s["sd"] is not None else np.abs(s["obs"]) * C.NOISE[grp]
        ax.errorbar(s["periods"], s["obs"], yerr=sd, fmt="k.", ms=4, lw=0.8, capsize=2, label="observed +- sd")
        if np.any(exist[i] < 1):
            ax.plot(s["periods"][exist[i] < 1], s["obs"][exist[i] < 1], "x", color="C3", ms=8, label="mode absent in some models")
        ax.set_xscale("log"); ax.set_xlabel("period (s)"); ax.set_ylabel("km/s"); ax.set_title(tag)
        if i == 0:
            ax.legend(fontsize=7)
        ax = fig.add_subplot(gs[3, i])
        res = (s["obs"] - fit_q[i][1]) / sd
        ax.plot(s["periods"], res, "ko", ms=4)
        for y in (-2, -1, 1, 2):
            ax.axhline(y, color="0.6", lw=0.6, ls="--" if abs(y) == 2 else "-")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xscale("log"); ax.set_xlabel("period (s)"); ax.set_ylabel("(obs - pred) / sd")
        ax.set_title(f"residuals, reduced chi2 = {chi2[i]:.2f}"); ax.set_ylim(-4, 4)

    fig.suptitle(f"{name}: synthetic Basel-1 test, slots " +
                 ", ".join(f"R{m}{'g' if g else 'p'}" for m, g in cfg["slots"]), y=0.995, fontsize=13)
    out = os.path.join(run_dir, "diagnostics.png")
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", out)
    return summary, q


def overlays(results, truth):
    """Vs(z) median + band of several runs over the truth."""
    z_true, vs_true, _ = truth
    groups = {}
    for band in ("full", "short", "long"):
        groups[f"overlay_{band}"] = [n for n in results if n.endswith(f"_{band}") and n.split("_")[0] in ("R0p", "R0g", "R0pg", "R01pg")]
    groups["overlay_variants"] = [n for n in ("R0pg_full", "R0pg_mismatch", "R0pg_modelerr", "R0pg_full_icov1") if n in results]
    os.makedirs(FIGS, exist_ok=True)
    for fname, names in groups.items():
        if len(names) < 2:
            continue
        fig, axes = plt.subplots(1, len(names), figsize=(3.4 * len(names), 7), sharey=True)
        for ax, n in zip(np.atleast_1d(axes), names):
            q = results[n]
            ax.fill_betweenx(ZGRID, q[0], q[2], color="C0", alpha=0.3)
            ax.plot(q[1], ZGRID, "-", color="C0", lw=1.2, label="median, 16-84 %")
            ax.plot(vs_true, ZGRID, "r-", lw=1.2, label="truth")
            ax.set_ylim(C.HMX, 0); ax.set_xlim(0.6, 3.6); ax.set_title(n); ax.set_xlabel("Vs (km/s)"); ax.grid(alpha=0.3)
        np.atleast_1d(axes)[0].set_ylabel("depth (km)"); np.atleast_1d(axes)[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIGS, fname + ".png"), dpi=120)
        plt.close(fig)
        print("wrote", os.path.join(FIGS, fname + ".png"))


def main():
    args = parse_args()
    runs = args.runs or sorted(os.path.join(RUNS, d) for d in os.listdir(RUNS)
                               if os.path.exists(os.path.join(RUNS, d, f"{C.BASE}_voro_sample.txt")))
    if not runs:
        sys.exit("no run directory with a sample file")
    truth = truth_profile()
    rng = np.random.default_rng(0)
    rows, results = [], {}
    for r in runs:
        try:
            s, q = analyze_run(r, args, truth, rng)
        except Exception as e:  # keep going through the matrix
            print(f"  {r}: FAILED: {e}")
            continue
        rows.append(s)
        results[s["run"]] = q
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
    if not args.no_overlays:
        overlays(results, truth)


if __name__ == "__main__":
    main()
