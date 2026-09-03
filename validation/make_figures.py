"""Figures for the SWD forward-model validation against disba.

Writes validation/figures/fig1..fig5 (PNG + PDF).

  fig1  dispersion curves, this code vs disba, with residuals
  fig2  summary over all models: root-value agreement vs mode-index agreement
  fig3  accuracy and cost as a function of the c-scan step
  fig4  why the mode index diverges - disba's duplicate roots
  fig5  group velocity: analytic U checked against d/dT of our own phase curve

Colour carries the MODE (categorical slots 1-3 of the reference palette,
validated all-pairs light: worst CVD dE 9.2, normal-vision dE 24.0); line style
carries the CODE.  Every mode curve is also directly labelled, which is the
relief required for the aqua slot's 2.74:1 contrast.
"""
import os
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models as M  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DRIVER = os.path.join(HERE, "disp_driver")
FIGDIR = os.path.join(HERE, "figures")
TMP = os.path.join(HERE, "out", "fig")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

# ---- design tokens -------------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8983"
GRID = "#e6e5e1"
MODE_C = ["#2a78d6", "#eb6834", "#1baf7a"]      # categorical slots 1-3
MODE_L = ["R0 (fundamental)", "R1 (1st overtone)", "R2 (2nd overtone)"]
MODE_S = ["R0", "R1", "R2"]
BAD = "#e34948"                                  # status: failure

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "text.color": INK, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": "#c9c8c3", "axes.linewidth": 0.8,
    "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "figure.dpi": 140,
})


def style(ax, grid=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(True, which="major", axis="both", zorder=0)
        ax.set_axisbelow(True)


# ---- data helpers --------------------------------------------------------
def fort(model, periods, maxmode, cmin, cmax, dc, dcov, iwarm=1, tag="f"):
    fm, fp = os.path.join(TMP, f"m{tag}.txt"), os.path.join(TMP, f"p{tag}.txt")
    M.write_model_file(fm, model)
    np.savetxt(fp, periods, fmt="%.9f")
    t0 = time.time()
    r = subprocess.run(
        [DRIVER, fm, fp, str(maxmode), f"{cmin:.9f}", f"{cmax:.9f}",
         f"{dc:.9f}", f"{dcov:.9f}", str(iwarm)],
        capture_output=True, text=True, check=True)
    dt = time.time() - t0
    c = np.full((maxmode + 1, len(periods)), np.nan)
    u = np.full((maxmode + 1, len(periods)), np.nan)
    for ln in r.stdout.strip().split("\n")[1:]:
        p = ln.split(",")
        if int(p[4]):
            c[int(p[1]), int(p[2]) - 1] = float(p[5])
            u[int(p[1]), int(p[2]) - 1] = float(p[6])
    return c, u, dt


def disba(model, periods, nm=3, dc=1e-5, want_u=False, dt=0.005):
    from disba import GroupDispersion, PhaseDispersion
    th = model[:, 0].copy(); th[-1] = max(th[-1], 1.0)
    args = (th, model[:, 2], model[:, 3], model[:, 1])
    pd = PhaseDispersion(*args, algorithm="dunkin", dc=dc)
    c = np.full((nm, len(periods)), np.nan)
    u = np.full((nm, len(periods)), np.nan)
    solvers = [(pd, c)]
    if want_u:
        solvers.append((GroupDispersion(*args, algorithm="dunkin", dc=dc, dt=dt), u))
    for solver, arr in solvers:
        for m in range(nm):
            try:
                cur = solver(periods, mode=m, wave="rayleigh")
                if len(cur.period):
                    arr[m, np.searchsorted(periods, cur.period)] = cur.velocity
            except Exception:
                pass
    return (c, u) if want_u else c


DC_FINE = 0.0001          # converged scan step, km/s
PER = np.logspace(np.log10(0.01), np.log10(0.6), 60)
FREQ = 1.0 / PER


def plot_profile(ax, model, zmax=None, color=INK2, lw=1.8, fill=True):
    """Vs(z) as a step profile, depth downwards. Last entry is the half-space."""
    h = model[:, 0] * 1000.0                      # m
    vs = model[:, 3] * 1000.0                     # m/s
    ztop = np.concatenate([[0.0], np.cumsum(h[:-1])])
    base = ztop[-1]
    if zmax is None:
        zmax = base * 1.35 if base > 0 else 100.0
    zz, vv = [], []
    for i in range(len(vs)):
        top = ztop[i]
        bot = ztop[i + 1] if i + 1 < len(ztop) else zmax
        if i == len(vs) - 1:
            bot = zmax
        zz += [top, bot]
        vv += [vs[i], vs[i]]
    if fill:
        ax.fill_betweenx(zz, 0, vv, color=color, alpha=0.08, lw=0, zorder=2)
    ax.plot(vv, zz, color=color, lw=lw, zorder=4, solid_joinstyle="miter")
    ax.axhline(base, color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=3)
    ax.annotate("half-space", (0.96, base), xycoords=("axes fraction", "data"),
                ha="right", va="bottom", fontsize=7, color=MUTED)
    ax.set_ylim(zmax, 0)
    ax.set_xlim(0, max(vv) * 1.12)
    return zmax


def label_curve(ax, x, y, text, color, frac=0.62, dy=6):
    """Direct label on a curve (relief for the low-contrast slot)."""
    ok = ~np.isnan(y)
    if ok.sum() < 3:
        return
    xs, ys = x[ok], y[ok]
    i = int(len(xs) * frac)
    ax.annotate(text, (xs[i], ys[i]), textcoords="offset points",
                xytext=(0, dy), ha="center", fontsize=8, color=color,
                fontweight="bold", zorder=6)


# =========================================================================
# fig 1 - dispersion curves + residuals
# =========================================================================
def fig1():
    mods = M.build_models()
    show = ["tailings_dry", "grad_sqrt_fine", "lvz_strong"]
    titles = ["Tailings profile (5 layers, 120-350 m/s)",
              "Compaction gradient (50 layers)",
              "Strong low-velocity zone (180 m/s at 15-40 m)"]

    fig, axes = plt.subplots(3, 3, figsize=(11.5, 8.6),
                             gridspec_kw=dict(height_ratios=[1.35, 2.3, 0.95],
                                              hspace=0.30, wspace=0.24))
    for j, (name, title) in enumerate(zip(show, titles)):
        model = mods[name]
        cmin, cmax = M.scan_window(model)
        cf, _, _ = fort(model, PER, 2, cmin, cmax, DC_FINE, DC_FINE)
        cd = disba(model, PER)
        prof, top, bot = axes[0, j], axes[1, j], axes[2, j]

        plot_profile(prof, model)
        style(prof)
        prof.set_xlabel("$V_S$  (m/s)", labelpad=2)
        if j == 0:
            prof.set_ylabel("depth  (m)")
        prof.set_title(title, color=INK, pad=8, fontsize=9.5)

        for m in range(3):
            col = MODE_C[m]
            top.plot(FREQ, cd[m] * 1000, lw=4.5, color=col, alpha=0.30,
                     solid_capstyle="round", zorder=2)
            top.plot(FREQ, cf[m] * 1000, lw=1.5, color=col, zorder=4)
            sub = slice(2, None, 7)
            top.plot(FREQ[sub], cf[m][sub] * 1000, "o", ms=3.6, mfc=SURFACE,
                     mec=col, mew=1.1, zorder=5)
            label_curve(top, FREQ, cf[m] * 1000, MODE_S[m], col,
                        frac=[0.62, 0.46, 0.34][m], dy=7)
            d = (cf[m] - cd[m]) * 1000
            bot.plot(FREQ, d, lw=1.3, color=col, zorder=4)

        top.set_xscale("log"); bot.set_xscale("log")
        top.tick_params(labelbottom=False)
        style(top); style(bot)
        if j == 0:
            top.set_ylabel("phase velocity  (m/s)")
            bot.set_ylabel("this code $-$ disba\n(m/s)")
        bot.set_xlabel("frequency  (Hz)")
        bot.axhspan(-0.02, 0.02, color=MUTED, alpha=0.16, lw=0, zorder=1)
        bot.axhline(0, color=MUTED, lw=0.7, zorder=2)
        bot.set_ylim(-0.055, 0.055)
        if j == 2:
            bot.annotate("root-finder tolerance\n($\\pm$0.02 m/s)", (0.97, 0.06),
                         xycoords="axes fraction", ha="right", va="bottom",
                         fontsize=7.5, color=MUTED)

    h = [plt.Line2D([], [], lw=4.5, color=MUTED, alpha=0.45,
                    solid_capstyle="round", label="disba 0.7.0 (reference)"),
         plt.Line2D([], [], lw=1.5, color=INK2, marker="o", ms=3.6,
                    mfc=SURFACE, mec=INK2, label="this code (DISPER-80)")]
    h += [plt.Line2D([], [], lw=2.4, color=MODE_C[m], label=MODE_L[m])
          for m in range(3)]
    fig.legend(handles=h, ncol=5, loc="lower center",
               bbox_to_anchor=(0.5, 0.005), columnspacing=2.0)
    fig.suptitle("Rayleigh dispersion: this code vs disba, 0-100 m models, "
                 f"scan step {DC_FINE*1000:.1f} m/s",
                 y=0.985, fontsize=11, color=INK)
    save(fig, "fig1_dispersion_validation")


# =========================================================================
# fig 2 - summary across all models
# =========================================================================
def fig2():
    mods = M.build_models()
    names, root_err, idx_ok, nvalid = [], [], [], []
    for name, model in mods.items():
        cmin, cmax = M.scan_window(model)
        cf, _, _ = fort(model, PER, 2, cmin, cmax, DC_FINE, DC_FINE)
        cd10 = disba(model, PER, nm=10)
        for m in range(3):
            errs, hits, n = [], 0, 0
            for i in range(len(PER)):
                v = cf[m, i]
                if np.isnan(v):
                    continue
                col = cd10[:, i]
                if np.all(np.isnan(col)):
                    continue
                k = int(np.nanargmin(np.abs(col - v)))
                errs.append(abs(col[k] - v) * 1000)
                hits += (k == m)
                n += 1
            if n:
                names.append(f"{name}  {MODE_S[m]}")
                root_err.append(max(errs))
                idx_ok.append(hits / n)
                nvalid.append(n)

    order = np.arange(len(names))[::-1]
    y = np.arange(len(names))
    fig, (a1, a2) = plt.subplots(
        1, 2, figsize=(11.5, 8.6), sharey=True,
        gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.06))

    cols = [MODE_C[int(n.split()[-1][1])] for n in names]
    a1.barh(y, np.array(root_err)[order], height=0.62,
            color=[cols[i] for i in order], zorder=3)
    a1.set_xscale("log")
    a1.set_xlim(1e-5, 2.0)
    a1.axvline(0.02, color=MUTED, lw=1.0, ls="--", zorder=4)
    a1.annotate("root-finder\ntolerance", (0.02, len(names) - 0.15),
                fontsize=7.5, color=MUTED, ha="center", va="bottom",
                annotation_clip=False)
    a1.set_xlabel("worst |Δc| to the nearest disba root  (m/s)")
    a1.set_title("Are the roots right?", color=INK, pad=10, loc="left")
    n_ok = sum(1 for e in root_err if e <= 0.02)
    n_big = sum(1 for e in root_err if e > 0.1)
    a1.annotate(f"{n_ok} of {len(root_err)} curves within the root-finder\n"
                f"tolerance; only {n_big} exceed 0.1 m/s",
                (0.985, 0.015), xycoords="axes fraction", ha="right",
                va="bottom", fontsize=7.8, color=INK2)

    ok = np.array(idx_ok)[order]
    a2.barh(y, ok * 100, height=0.62,
            color=[cols[i] if ok[k] > 0.999 else BAD
                   for k, i in enumerate(order)], zorder=3)
    a2.set_xlim(0, 104)
    a2.set_xlabel("periods where the mode index matches disba  (%)")
    a2.set_title("Is the mode numbering the same?", color=INK, pad=10, loc="left")
    for k, i in enumerate(order):
        if ok[k] <= 0.999:
            a2.annotate("disba duplicates roots here", (ok[k] * 100 + 2, k),
                        va="center", fontsize=7.5, color=BAD)

    a1.set_yticks(y)
    a1.set_yticklabels([names[i] for i in order], fontsize=7.6)
    a1.set_ylim(-0.8, len(names) - 0.2)
    for ax in (a1, a2):
        style(ax)
        ax.grid(axis="y", visible=False)
    fig.suptitle("All 16 models, modes R0-R2   —   left: do we find the same "
                 "root VALUES?    right: do we NUMBER them the same way?",
                 y=0.965, fontsize=10.5, color=INK)
    save(fig, "fig2_residual_summary", tight_rect=(0, 0, 1, 0.95))


# =========================================================================
# fig 3 - accuracy and cost vs the scan step
# =========================================================================
def fig3():
    mods = M.build_models()
    usable = [n for n in mods if n not in ("2layer_normal", "near_homogeneous")]
    steps = [0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005]
    worst, med, cost = [], [], []
    ref = {}
    for name in usable:
        model = mods[name]
        cmin, cmax = M.scan_window(model)
        ref[name] = (disba(model, PER), cmin, cmax,
                     fort(model, PER, 2, cmin, cmax, 2e-5, 2e-5)[0])
    for dcv in steps:
        w, alld, tt = 0.0, [], 0.0
        for name in usable:
            model = mods[name]
            cd, cmin, cmax, cfine = ref[name]
            c, _, dt = fort(model, PER, 2, cmin, cmax, dcv, dcv)
            tt += dt
            ok = (~np.isnan(c) & ~np.isnan(cd) & ~np.isnan(cfine)
                  & (np.abs(cfine - cd) * 1000 < 1.0))
            if ok.sum():
                d = np.abs(c[ok] - cd[ok]) * 1000
                w = max(w, d.max()); alld.extend(d.tolist())
        worst.append(w); med.append(np.median(alld))
        cost.append(tt / len(usable) * 1000)

    x = np.array(steps) * 1000
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.4, 4.2),
                                 gridspec_kw=dict(wspace=0.26))
    a1.plot(x, worst, lw=2, color=BAD, marker="o", ms=6, mfc=SURFACE,
            mew=1.6, mec=BAD, zorder=4)
    a1.plot(x, med, lw=2, color=MODE_C[0], marker="s", ms=5.5, mfc=SURFACE,
            mew=1.6, mec=MODE_C[0], zorder=4)
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.invert_xaxis()
    a1.axvspan(5.0, 1.0, color=BAD, alpha=0.07, lw=0, zorder=1)
    a1.annotate("values suggested\nfor near-surface\nwork (commit dad67af)",
                (2.4, 4e-3), fontsize=7.6, color=BAD, ha="center")
    a1.annotate("worst case", (0.35, 3.0), color=BAD, fontsize=8.5,
                fontweight="bold")
    a1.annotate("median", (0.35, 3e-4), color=MODE_C[0], fontsize=8.5,
                fontweight="bold")
    a1.set_xlabel("c-scan step  dc = dc_over  (m/s)   $\\longrightarrow$ finer")
    a1.set_ylabel("|Δc| vs disba  (m/s)")
    a1.set_title("Accuracy is set by the scan step", color=INK, loc="left", pad=8)

    a2.plot(x, cost, lw=2, color=MODE_C[2], marker="o", ms=6, mfc=SURFACE,
            mew=1.6, mec=MODE_C[2], zorder=4)
    a2.set_xscale("log"); a2.invert_xaxis()
    a2.set_ylim(0, max(cost) * 1.18)
    a2.axvspan(5.0, 1.0, color=BAD, alpha=0.07, lw=0, zorder=1)
    a2.set_xlabel("c-scan step  dc = dc_over  (m/s)   $\\longrightarrow$ finer")
    a2.set_ylabel("forward model  (ms per model)")
    a2.set_title("and the cost is modest", color=INK, loc="left", pad=8)
    for xx, yy in zip(x, cost):
        a2.annotate(f"{yy:.0f}", (xx, yy), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=7.4, color=INK2)
    for ax in (a1, a2):
        style(ax)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{v:g}" for v in x])
        ax.minorticks_off()
    fig.suptitle("3 modes x 60 periods, 14 models, scan window fitted per model",
                 y=1.0, fontsize=9.5, color=MUTED)
    save(fig, "fig3_scan_step")


# =========================================================================
# fig 4 - the mode-index divergence is disba's duplicate roots
# =========================================================================
def fig4():
    mods = M.build_models()
    name, T = "lid_over_channel", 0.050
    model = mods[name]
    cmin, cmax = M.scan_window(model)
    per = np.array([T])
    cf, _, _ = fort(model, per, 9, cmin, cmax, DC_FINE, DC_FINE)
    cd = disba(model, per, nm=10)
    fl = [v * 1000 for v in cf[:, 0] if not np.isnan(v)]
    dl = [v * 1000 for v in cd[:, 0] if not np.isnan(v)]

    fig, ax = plt.subplots(figsize=(11.0, 3.9))
    yF, yD = 1.0, 0.0
    ax.plot([min(fl + dl) - 6, max(fl + dl) + 6], [yF, yF], color=GRID, lw=1.2, zorder=1)
    ax.plot([min(fl + dl) - 6, max(fl + dl) + 6], [yD, yD], color=GRID, lw=1.2, zorder=1)

    # connect equal root values between the two lanes
    for v in fl:
        m = [w for w in dl if abs(w - v) < 1.0]
        if m:
            ax.plot([v, m[0]], [yF, yD], color=MUTED, lw=0.7, ls=":", zorder=2)

    def cluster(vals, tol=1.0):
        """Group roots that are the same value; keep the indices in each group."""
        out = []
        for i, v in enumerate(vals):
            if out and abs(v - out[-1][0]) < tol:
                out[-1][1].append(i)
            else:
                out.append([v, [i]])
        return out

    for v, idxs in cluster(fl):
        ax.plot(v, yF, "o", ms=8, color=MODE_C[0], zorder=4)
        ax.annotate(",".join(map(str, idxs)), (v, yF),
                    textcoords="offset points", xytext=(0, 11), ha="center",
                    fontsize=8, color=MODE_C[0], fontweight="bold")
    for v, idxs in cluster(dl):
        dup = len(idxs) > 1
        col = BAD if dup else MODE_C[1]
        ax.plot(v, yD, "o", ms=8, color=col, zorder=4)
        ax.annotate(",".join(map(str, idxs)), (v, yD),
                    textcoords="offset points", xytext=(0, -17), ha="center",
                    fontsize=8, color=col, fontweight="bold")
        if dup:
            ax.annotate("one root,\ntwo mode slots", (v, yD),
                        textcoords="offset points", xytext=(0, -31),
                        ha="center", va="top", fontsize=7, color=BAD)

    ax.set_yticks([yD, yF])
    ax.set_yticklabels(["disba 0.7.0", "this code"], fontsize=9.5)
    ax.set_ylim(-0.75, 1.55)
    ax.set_xlabel("Rayleigh phase velocity  (m/s)")
    ax.set_title(f"Every trapped root at {1/T:.0f} Hz, model "
                 f"'{name}' (400 / 180 / 600 m/s)   —   numbers are mode indices",
                 color=INK, loc="left", pad=12)
    style(ax)
    ax.grid(axis="y", visible=False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.annotate("Both codes find the same root VALUES.  disba spends two mode slots\n"
                "on one root (red), which shifts every index above it by one — and\n"
                "its list runs out before ours, so the top roots have no counterpart.",
                (0.995, 0.99), xycoords="axes fraction", ha="right", va="top",
                fontsize=8, color=INK2)
    save(fig, "fig4_mode_index")


# =========================================================================
# fig 5 - group velocity
# =========================================================================
def fig5():
    mods = M.build_models()
    h = 2e-3
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4),
                             gridspec_kw=dict(wspace=0.22))
    for ax, name, title in zip(
            axes, ["tailings_dry", "lvz_weak"],
            ["Tailings profile", "Low-velocity zone (200 m/s at 20-40 m)"]):
        model = mods[name]
        cmin, cmax = M.scan_window(model)
        grid = np.sort(np.concatenate([PER * (1 - h), PER, PER * (1 + h)]))
        cg, ug, _ = fort(model, grid, 2, cmin, cmax, DC_FINE, DC_FINE, tag="g")
        idx = {round(float(t), 12): j for j, t in enumerate(grid)}
        _, ud = disba(model, PER, nm=3, want_u=True)

        umax = 0.0
        for m in range(3):
            col = MODE_C[m]
            ua, ufd = np.full(len(PER), np.nan), np.full(len(PER), np.nan)
            for i, T in enumerate(PER):
                j0 = idx[round(float(T), 12)]
                jm = idx[round(float(T * (1 - h)), 12)]
                jp = idx[round(float(T * (1 + h)), 12)]
                cm, c0, cp = cg[m, jm], cg[m, j0], cg[m, jp]
                ua[i] = ug[m, j0] * 1000
                if np.isnan(cm) or np.isnan(c0) or np.isnan(cp):
                    continue
                if abs(cp - cm) > 0.05 * c0:
                    continue
                den = 1.0 + (T / c0) * (cp - cm) / (T * 2 * h)
                if abs(den) > 1e-6:
                    ufd[i] = c0 / den * 1000
            ax.plot(FREQ, ufd, lw=4.5, color=col, alpha=0.30,
                    solid_capstyle="round", zorder=2)
            ax.plot(FREQ, ua, lw=1.5, color=col, zorder=4)
            # drawn on top so it is visible where it agrees, not just where it fails
            ax.plot(FREQ, ud[m] * 1000, lw=1.4, ls=(0, (2.5, 2.5)), color=BAD,
                    zorder=7)
            umax = max(umax, np.nanmax(ua))
            label_curve(ax, FREQ, ua, MODE_S[m], col,
                        frac=0.90, dy=8)
        ax.set_xscale("log")
        style(ax)
        ax.set_xlabel("frequency  (Hz)")
        ax.set_title(title, color=INK, loc="left", pad=8)
        # scale to our own curves; disba's spikes are allowed to run off
        ax.set_ylim(0, umax * 1.20)
        ins = ax.inset_axes([0.68, 0.60, 0.30, 0.37])
        plot_profile(ins, model, lw=1.3)
        ins.tick_params(labelsize=6.5, length=2, pad=1.5)
        ins.set_xlabel("$V_S$ (m/s)", fontsize=6.5, labelpad=1)
        ins.set_ylabel("depth (m)", fontsize=6.5, labelpad=1)
        for s in ("top", "right"):
            ins.spines[s].set_visible(False)
        for t in ins.texts:
            t.set_visible(False)
    axes[0].set_ylabel("group velocity  (m/s)")
    axes[0].annotate("all three curves agree", (0.035, 0.97),
                     xycoords="axes fraction", ha="left", va="top",
                     fontsize=8.5, color=INK2)
    axes[1].annotate("disba's GroupDispersion collapses to 10-30 m/s here,\n"
                     "which is not physical — and it does so exactly where\n"
                     "the two PHASE curves agree to 1e-6 km/s",
                     xy=(6.0, 13.0), xycoords="data",
                     xytext=(0.055, 0.20), textcoords="axes fraction",
                     ha="left", va="top", fontsize=7.6, color=BAD,
                     arrowprops=dict(arrowstyle="->", color=BAD, lw=1.0,
                                     shrinkB=3, connectionstyle="arc3,rad=-0.25"))

    hh = [plt.Line2D([], [], lw=4.5, color=MUTED, alpha=0.45,
                     solid_capstyle="round",
                     label="$U = c\\,/\\,(1+(T/c)\\,dc/dT)$ from OUR own phase curve"),
          plt.Line2D([], [], lw=1.5, color=INK2,
                     label="our analytic $U$ (DISPER-80 eigenfunctions)"),
          plt.Line2D([], [], lw=1.4, ls=(0, (2.5, 2.5)), color=BAD,
                     label="disba GroupDispersion (agrees left, fails right)")]
    fig.legend(handles=hh, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, -0.04), columnspacing=2.0)
    fig.suptitle("Group velocity: the analytic U agrees with a derivative of our "
                 "own phase curve to 0.006-0.09 m/s on smooth models",
                 y=1.0, fontsize=10.5, color=INK)
    save(fig, "fig5_group_velocity")


def save(fig, stem, tight_rect=None):
    if tight_rect:
        fig.tight_layout(rect=tight_rect)
    else:
        fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGDIR, f"{stem}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figures/{stem}.png / .pdf")


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3", "4", "5"]
    for n in which:
        print(f"fig{n} ...")
        globals()[f"fig{n}"]()
