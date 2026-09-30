"""Figures: the rjhist-style diagnostics of one run, data fits, overlays.

rjhist_figure draws three blocks (the layout of the original
rf_plot_rjhist_varpar figures, condensed):
  posteriors    Vs(z) per-depth histogram ('bone' flipped, white zero bin) with
                the 16/50/84 % profile, node-depth (interface) density,
                Vp/Vs(z), p(k) with its Poisson prior, hierarchical sigma per slot
  convergence   logL and k per cold chain vs sample index (burn-in marked),
                running median Vs at three depths per chain, the three L1
                stationarity metrics (<= 0.05 comfortable, >= 0.15 re-run)
  data fits     per slot: observed curve (+- sd when ICOV_SWD = 3), the MAP
                prediction the engine wrote (IMAP), the posterior-predictive
                16/50/84 % band when validation/disp_driver is built, residuals
A `truth` dict adds overlays (synthetic tests): {'nodes': (z, dVs, dVpVs),
'curves': [per slot (period, value) arrays], 'logL': float, 'sigma': [per slot],
'label': str}.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402

from . import forward, io, model, posterior  # noqa: E402

PROBE_FRACTIONS = (0.04, 0.2, 0.6)     # of hmx: depths of the running-median traces


def _truth_profiles(truth, run, zgrid):
    if not truth or "nodes" not in truth:
        return None, None, None
    z, dvs, dvpvs = truth["nodes"]
    return (np.asarray(z), model.vs_on_grid(z, dvs, zgrid, run.cfg, run.vel_ref),
            model.vpvs_on_grid(z, dvpvs, zgrid, run.cfg, run.vel_ref))


def vs_limits(run, zgrid, pad=0.05):
    """The Vs prior box over the depth grid: reference +- dVs."""
    ref = np.array([model.getref(z, run.vel_ref)[0] for z in zgrid])
    lo, hi = ref.min() - run.cfg["dVs"], ref.max() + run.cfg["dVs"]
    return max(lo, 0.05) - pad, hi + pad


def vpvs_limits(run, zgrid, pad=0.02):
    ref = np.array([model.getref(z, run.vel_ref)[1] for z in zgrid])
    if run.cfg["I_VPVS"] != 1:
        return 1.7, 1.8
    return ref.min() - run.cfg["dVpVs"] - pad, ref.max() + run.cfg["dVpVs"] + pad


def sigma_axis(run, i):
    """(scale, label, prior bounds) of the sigma of slot i for display."""
    cfg = run.cfg
    lo, hi = cfg["SDMN_SWD"][i], cfg["SDMX_SWD"][i]
    if cfg["ICOV_SWD"] == 3:
        return 1.0, "sigma multiplier of the sd files", lo, hi
    if cfg["ICOV_SWD"] == 1 and cfg["IMAGSCALE"] == 1:
        return 100.0, "sigma (% of datum)", 100 * lo, 100 * hi
    if cfg["ICOV_SWD"] == 2:
        return 1.0, "sigma scale of Cd", lo, hi
    return 1.0, "sigma (km/s)", lo, hi


def rjhist_figure(run, *, truth=None, burn=0.3, nfwd=300, nprof=4000, zmax=None, out=None,
                  title=None, rng=None, driver=None):
    """Diagnostics figure of one run. Returns (fig, summary dict, extras) with
    extras = {zgrid, q (16/50/84 % Vs), vs, vpvs, znodes, zh, zbins, post, sig,
    chi2, pred} for further analysis (validation/basel_group/analyze.py)."""
    if run.sample is None:
        raise ValueError(f"{run.run_dir}: no sample file")
    rng = rng or np.random.default_rng(0)
    cfg, dat = run.cfg, run.sample
    nmode = cfg["NMODE"]
    zgrid = posterior.depth_grid(run, zmax=zmax)
    zedges = np.r_[zgrid - 0.005, zgrid[-1] + 0.005]
    zbins = np.arange(0.0, zgrid[-1] + 1e-9, 0.02)
    zc = 0.5 * (zbins[1:] + zbins[:-1])
    vlo, vhi = vs_limits(run, zgrid)
    vs_bins = np.linspace(vlo, vhi, 151)
    vpvs_bins = np.linspace(*vpvs_limits(run, zgrid), 81)
    z_true, vs_true, vpvs_true = _truth_profiles(truth, run, zgrid)
    truth_logl = truth.get("logL", np.nan) if truth else np.nan
    truth_label = truth.get("label", "truth") if truth else "truth"

    post, _ = io.burn_in_split(dat, burn)
    chains = np.unique(dat[:, -1]).astype(int)

    # ---- posteriors ----
    vs, vpvs, znodes = posterior.profiles(post, run, zgrid, nprof, rng)
    q = posterior.quantiles(vs)
    vs_h = posterior.per_depth_hist(vs, vs_bins)
    vpvs_h = posterior.per_depth_hist(vpvs, vpvs_bins)
    zh = posterior.interface_density(znodes, zbins)
    kpost = post[:, 3].astype(int)
    kk, pk = posterior.poisson_prior(cfg)
    sig = posterior.sigma_columns(post, run)

    # ---- convergence ----
    metrics = posterior.convergence_metrics(dat, post, q, run, zgrid, burn, nprof, rng)
    move = posterior.move_fraction(post, run)

    # ---- data fits ----
    pred = forward.predict_slots(run, post, nfwd, rng, driver) if forward.available(driver) else None
    fit_q, chi2, exist = [], [], []
    for i, s in enumerate(run.slots):
        if pred is None:
            fit_q.append(None); exist.append(None); chi2.append(np.nan)
            continue
        with np.errstate(all="ignore"):
            fq = np.nanpercentile(pred[i], [16, 50, 84], axis=0)
        fit_q.append(fq)
        exist.append(np.mean(np.isfinite(pred[i]), axis=0))
        sd = s["sd"] if s["sd"] is not None else np.abs(s["obs"]) * np.median(sig[:, i]) if cfg["IMAGSCALE"] == 1 else None
        chi2.append(float(np.nanmean(((s["obs"] - fq[1]) / sd) ** 2)) if sd is not None else np.nan)

    summary = {"run": run.name, "rows": len(dat), "post": len(post), "chains": len(chains),
               "k_median": float(np.median(kpost)), "move_frac": move, **metrics,
               "logL_median": float(np.median(post[:, 0])), "truth_logL": truth_logl}
    for i, (mode, grp) in enumerate(cfg["slots"]):
        tag = io.slot_label(mode, grp, short=True)
        summary[f"sigma_{tag}"] = float(np.median(sig[:, i])) if cfg["ICOV_SWD"] else np.nan
        summary[f"chi2_{tag}"] = chi2[i]

    # ================= figure =================
    ncol = max(5, nmode)
    fig = plt.figure(figsize=(3.6 * ncol, 15))
    gs = gridspec.GridSpec(4, ncol, figure=fig, height_ratios=[1.6, 1.0, 1.0, 0.7], hspace=0.42, wspace=0.32)
    cmap = plt.get_cmap("bone_r").copy(); cmap.set_under("white")
    tab = plt.get_cmap("tab10")

    ax = fig.add_subplot(gs[0, 0])
    ax.pcolormesh(vs_bins, zedges, vs_h, cmap=cmap, vmin=1e-3, vmax=1, shading="flat")
    ax.plot(q[1], zgrid, "w--", lw=1, label="posterior median")
    ax.plot(q[0], zgrid, "-", color="0.6", lw=0.6); ax.plot(q[2], zgrid, "-", color="0.6", lw=0.6)
    if vs_true is not None:
        ax.plot(vs_true, zgrid, "r-", lw=1.4, label=truth_label)
    if run.map_voro is not None and truth is None:
        zm, dm, _ = model.map_nodes(run.map_voro, cfg)
        ax.plot(model.vs_on_grid(zm, dm, zgrid, cfg, run.vel_ref), zgrid, "-", color="C1", lw=0.8, label="map_voro")
    ax.set_ylim(zgrid[-1], 0); ax.set_xlim(vs_bins[0], vs_bins[-1]); ax.set_xlabel("Vs (km/s)"); ax.set_ylabel("depth (km)")
    ax.set_title("Vs posterior"); ax.legend(fontsize=7, loc="lower left")

    ax = fig.add_subplot(gs[0, 1])
    ax.fill_betweenx(zc, 0, zh, color="0.5", alpha=0.7)
    if z_true is not None:
        for zt in z_true[1:]:
            ax.axhline(zt, color="r", lw=0.6, alpha=0.7)
    ax.set_ylim(zgrid[-1], 0); ax.set_xlabel("node-depth density (1/km)"); ax.set_title("interface probability")

    ax = fig.add_subplot(gs[0, 2])
    if cfg["I_VPVS"] == 1:
        ax.pcolormesh(vpvs_bins, zedges, vpvs_h, cmap=cmap, vmin=1e-3, vmax=1, shading="flat")
        if vpvs_true is not None:
            ax.plot(vpvs_true, zgrid, "r-", lw=1.4)
        ax.set_xlabel("Vp/Vs"); ax.set_title("Vp/Vs posterior")
    else:
        ax.text(0.5, 0.5, "Vp/Vs fixed (I_VPVS = -1: Vp = 1.75 Vs)", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Vp/Vs")
    ax.set_ylim(zgrid[-1], 0)

    ax = fig.add_subplot(gs[0, 3])
    kh = np.array([(kpost == k).mean() for k in kk])
    ax.bar(kk, kh, color="0.4", label="posterior")
    if cfg["IPOIPR"] == 1:
        ax.plot(kk, pk, "o-", color="C1", ms=3, label=f"Poisson({cfg['lambda']:g}) prior")
    if z_true is not None:
        ax.axvline(len(z_true), color="r", lw=1, label=f"truth k = {len(z_true)}")
    ax.set_xlabel("number of nodes k"); ax.set_title("p(k)"); ax.legend(fontsize=7)

    ax = fig.add_subplot(gs[0, 4])
    if cfg["ICOV_SWD"] == 0:
        ax.text(0.5, 0.5, "ICOV_SWD = 0: sigma integrated out\n(no hierarchical parameter)", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title("hierarchical sigma"); ax.set_xticks([]); ax.set_yticks([])
    else:
        xl = ""
        for i, (mode, grp) in enumerate(cfg["slots"]):
            scale, xl, lo, hi = sigma_axis(run, i)
            x = scale * sig[:, i]
            edges = np.linspace(lo, min(hi, np.percentile(x, 99.5) * 1.3), 60)
            h, e = np.histogram(x, bins=edges, density=True)
            ax.step(e[:-1], h, where="post", color=f"C{i}", label=f"{io.slot_label(mode, grp)}: median {np.median(x):.3g}")
            if truth and truth.get("sigma") is not None:
                ax.axvline(scale * truth["sigma"][i], color=f"C{i}", ls=":", lw=1)
        ax.set_xlabel(xl); ax.set_title("hierarchical sigma per slot" + (" (dotted = expected)" if truth and truth.get("sigma") is not None else ""))
        ax.legend(fontsize=7)

    # ---- row 2: convergence ----
    ax = fig.add_subplot(gs[1, 0:2])
    for j, c in enumerate(chains):
        sub = dat[dat[:, -1] == c]
        ax.plot(np.arange(len(sub)), sub[:, 0], lw=0.5, alpha=0.8, color=tab(j % 10), label=f"chain (rank {c})" if j < 6 else None)
        ax.axvline(int(len(sub) * burn), color=tab(j % 10), lw=0.5, ls="--")
    if np.isfinite(truth_logl):
        ax.axhline(truth_logl, color="r", lw=1, label=f"{truth_label} logL = {truth_logl:.1f}")
    lo = np.percentile(post[:, 0], 0.5)
    ax.set_ylim(lo - 0.3 * abs(lo - post[:, 0].max()) - 5, post[:, 0].max() + 5)
    ax.set_xlabel("sample index within chain"); ax.set_ylabel("logL"); ax.set_title("logL per cold chain (dashed = burn-in cut)")
    ax.legend(fontsize=7, ncol=2)

    ax = fig.add_subplot(gs[1, 2])
    for j, c in enumerate(chains):
        sub = dat[dat[:, -1] == c]
        ax.plot(np.arange(len(sub)), sub[:, 3], lw=0.5, alpha=0.8, color=tab(j % 10))
    if z_true is not None:
        ax.axhline(len(z_true), color="r", lw=1)
    ax.set_xlabel("sample index within chain"); ax.set_ylabel("k"); ax.set_title("k per chain")

    ax = fig.add_subplot(gs[1, 3])
    probes = [f * zgrid[-1] for f in PROBE_FRACTIONS]
    for j, c in enumerate(chains):
        sub = dat[dat[:, -1] == c]
        step = max(1, len(sub) // 3000)
        pv = posterior.probe_vs(sub[::step], run, probes)
        idx, rm = posterior.running_median(pv)
        for d in range(len(probes)):
            ax.plot(idx * step, rm[:, d], lw=0.8, color=tab(j % 10), ls=["-", "--", ":"][d])
    for d, zp in enumerate(probes):
        if vs_true is not None:
            ax.axhline(np.interp(zp, zgrid, vs_true), color="r", lw=0.8, ls=["-", "--", ":"][d], label=f"{truth_label} at {zp:.2g} km")
        else:
            ax.plot([], [], "k", ls=["-", "--", ":"][d], label=f"{zp:.2g} km")
    ax.set_xlabel("sample index within chain"); ax.set_ylabel("running median Vs (km/s)")
    ax.set_title("running median Vs per chain"); ax.legend(fontsize=7)

    ax = fig.add_subplot(gs[1, 4]); ax.axis("off")
    txt = (f"{run.name}\n\nrows {len(dat)}, chains {len(chains)}\npost-burn-in ({burn:.0%}) {len(post)}\n\n"
           f"stationarity (1st vs 2nd half)  {metrics['stationarity_L1']:.3f}\n"
           f"between-chain (worst pair)      {metrics['between_chain_L1']:.3f}\n"
           f"burn-in 30 % vs 60 %            {metrics['burn30v60_L1']:.3f}\n"
           f"(L1 in 68 % half-width units;\n <= 0.05 good, >= 0.15 re-run)\n\n"
           f"moved-sample fraction  {move:.3f}\nk median {np.median(kpost):.0f}"
           + (f" (truth {len(z_true)})" if z_true is not None else "") + "\n"
           f"logL median {np.median(post[:, 0]):.1f}" + (f"  truth {truth_logl:.1f}" if np.isfinite(truth_logl) else "") + "\n")
    for i, (mode, grp) in enumerate(cfg["slots"]):
        line = f"{io.slot_label(mode, grp, short=True)}: sigma {np.median(sig[:, i]):.3g}" if cfg["ICOV_SWD"] else io.slot_label(mode, grp, short=True)
        if np.isfinite(chi2[i]):
            line += f", chi2_red {chi2[i]:.2f}"
        txt += line + "\n"
    ax.text(0.0, 1.0, txt, va="top", ha="left", family="monospace", fontsize=8.5, transform=ax.transAxes)

    # ---- rows 3-4: data fits ----
    for i, s in enumerate(run.slots):
        mode, grp = cfg["slots"][i]
        ax = fig.add_subplot(gs[2, i])
        if fit_q[i] is not None:
            ax.fill_between(s["periods"], fit_q[i][0], fit_q[i][2], color="C0", alpha=0.3, label="posterior predictive 16-84 %")
            ax.plot(s["periods"], fit_q[i][1], "-", color="C0", lw=1.2, label="posterior predictive median")
            if np.any(exist[i] < 1):
                ax.plot(s["periods"][exist[i] < 1], s["obs"][exist[i] < 1], "x", color="C3", ms=8, label="mode absent in some models")
        if run.mappred is not None:
            ax.plot(s["periods"], run.mappred[i], "-", color="C1", lw=1, label="MAP prediction (IMAP)")
        if truth and truth.get("curves") is not None and truth["curves"][i] is not None:
            tc = truth["curves"][i]
            ax.plot(tc[:, 0], tc[:, 1], "r-", lw=1, label=f"{truth_label} (noise free)")
        if s["sd"] is not None:
            ax.errorbar(s["periods"], s["obs"], yerr=s["sd"], fmt="k.", ms=4, lw=0.8, capsize=2, label="observed +- sd")
        else:
            ax.plot(s["periods"], s["obs"], "k.", ms=4, label="observed")
        ax.set_xscale("log"); ax.set_xlabel("period (s)"); ax.set_ylabel("km/s")
        ax.set_title(f"{io.slot_label(mode, grp)} velocity")
        if i == 0:
            ax.legend(fontsize=7)
        ax = fig.add_subplot(gs[3, i])
        ref = fit_q[i][1] if fit_q[i] is not None else (run.mappred[i] if run.mappred is not None else None)
        if ref is not None:
            if s["sd"] is not None:
                res, yl = (s["obs"] - ref) / s["sd"], "(obs - pred) / sd"
            else:
                res, yl = 100 * (s["obs"] - ref) / s["obs"], "(obs - pred) / obs (%)"
            ax.plot(s["periods"], res, "ko", ms=4)
            if s["sd"] is not None:
                for y in (-2, -1, 1, 2):
                    ax.axhline(y, color="0.6", lw=0.6, ls="--" if abs(y) == 2 else "-")
                ax.set_ylim(-4, 4)
            ax.axhline(0, color="k", lw=0.6); ax.set_ylabel(yl)
            ax.set_title("residuals" + (f", reduced chi2 = {chi2[i]:.2f}" if np.isfinite(chi2[i]) else
                                        (" vs MAP" if fit_q[i] is None else "")))
        else:
            ax.text(0.5, 0.5, "no prediction:\nbuild validation/disp_driver\nor run IMAP", ha="center", va="center", transform=ax.transAxes)
        ax.set_xscale("log"); ax.set_xlabel("period (s)")

    fig.suptitle(title or f"{run.name}: slots " + ", ".join(io.slot_label(m, g, short=True) for m, g in cfg["slots"]),
                 y=0.995, fontsize=13)
    if out:
        fig.savefig(out, dpi=120, bbox_inches="tight")
    extras = {"zgrid": zgrid, "q": q, "vs": vs, "vpvs": vpvs, "znodes": znodes, "zh": zh, "zbins": zbins,
              "post": post, "sig": sig, "chi2": chi2, "pred": pred, "vs_true": vs_true, "z_true": z_true}
    return fig, summary, extras


def data_fit_figure(run, out=None):
    """Observed curves vs the MAP prediction rows the engine wrote in IMAP mode
    (<base>_mappredSWD.dat), residuals per slot."""
    if run.mappred is None:
        raise ValueError(f"{run.run_dir}: no {run.base}_mappredSWD.dat (run the engine with IMAP 1)")
    n = len(run.slots)
    fig, axes = plt.subplots(2, n, figsize=(4.5 * n, 7), squeeze=False)
    for i, s in enumerate(run.slots):
        mode, grp = run.cfg["slots"][i]
        ax = axes[0, i]
        if s["sd"] is not None:
            ax.errorbar(s["periods"], s["obs"], yerr=s["sd"], fmt="k.", ms=4, lw=0.8, capsize=2, label="observed +- sd")
        else:
            ax.plot(s["periods"], s["obs"], "k.", ms=4, label="observed")
        ax.plot(s["periods"], run.mappred[i], "-", color="C1", lw=1.2, label="MAP prediction")
        ax.set_xscale("log"); ax.set_ylabel("km/s"); ax.set_title(f"{io.slot_label(mode, grp)} velocity"); ax.legend(fontsize=8)
        ax = axes[1, i]
        if s["sd"] is not None:
            ax.plot(s["periods"], (s["obs"] - run.mappred[i]) / s["sd"], "ko", ms=4); ax.set_ylabel("(obs - MAP) / sd")
            for y in (-2, -1, 1, 2):
                ax.axhline(y, color="0.6", lw=0.6, ls="--" if abs(y) == 2 else "-")
        else:
            ax.plot(s["periods"], 100 * (s["obs"] - run.mappred[i]) / s["obs"], "ko", ms=4); ax.set_ylabel("(obs - MAP) / obs (%)")
        ax.axhline(0, color="k", lw=0.6); ax.set_xscale("log"); ax.set_xlabel("period (s)")
    fig.suptitle(f"{run.name}: MAP data fit")
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=120)
    return fig


def overlay_figure(entries, zgrid, truth_vs=None, vs_lim=None, out=None, truth_label="truth"):
    """entries: list of (label, q) with q the 16/50/84 % Vs profiles on zgrid."""
    fig, axes = plt.subplots(1, len(entries), figsize=(3.4 * len(entries), 7), sharey=True, squeeze=False)
    for ax, (label, q) in zip(axes[0], entries):
        ax.fill_betweenx(zgrid, q[0], q[2], color="C0", alpha=0.3)
        ax.plot(q[1], zgrid, "-", color="C0", lw=1.2, label="median, 16-84 %")
        if truth_vs is not None:
            ax.plot(truth_vs, zgrid, "r-", lw=1.2, label=truth_label)
        ax.set_ylim(zgrid[-1], 0); ax.set_title(label); ax.set_xlabel("Vs (km/s)"); ax.grid(alpha=0.3)
        if vs_lim is not None:
            ax.set_xlim(*vs_lim)
    axes[0, 0].set_ylabel("depth (km)"); axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=120)
    return fig


def truth_from_map_voro(path, cfg):
    """A truth model written in map_voro layout -> (z, dVs, dVpVs)."""
    return model.map_nodes(io.read_map_voro(path, cfg), cfg)


def truth_curves_from_csv(path, nslots):
    """validation/basel_group truth_curves.csv (slot, mode, grp, period, true, obs, sd)
    -> list per slot of (period, true) arrays."""
    rows = np.genfromtxt(path, delimiter=",", names=True)
    out = []
    for i in range(nslots):
        r = rows[rows["slot"] == i + 1]
        out.append(np.column_stack([r["period"], r["true"]]) if len(r) else None)
    return out


def save_close(fig, out):
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
