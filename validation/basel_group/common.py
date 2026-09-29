"""Shared pieces of the Basel-1 / Otterbach-2 group-velocity synthetic test
harness (validation/basel_group).

Everything here mirrors the ENGINE's own conventions (src/loglhood.f90
MAKE_CURMOD + the deep tail of LOGLHOOD_SWD, src/read_input.f90 file naming,
the sample-file layout of SAVESAMPLE), so that the synthetic truth lies inside
the model space the sampler explores and the posterior can be forwarded with
exactly the stack the engine would build.

Layer stacks are arrays of shape (n, 4): thickness km, density g/cc, Vp km/s,
Vs km/s, half-space last (thickness 0) -- the disp_driver model block.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
ENGINE = os.path.join(REPO, "src", "bin", "prjmh_temper_rf")
DISP_DRIVER = os.path.join(REPO, "validation", "disp_driver")

WELL_DIR = "/Volumes/T7blue/riehen-data/well-data"
MICHEL_MODEL = os.path.join(WELL_DIR, "Michel2016_gpdc.model")        # thick m, Vp, Vs m/s, rho kg/m3
MICHEL_CSV = os.path.join(WELL_DIR, "Vsmodel_well_Basel1_Otterbach_Michel2016.csv")
MICHEL_DISP = os.path.join(WELL_DIR, "Michel2016_gpdc.disp")          # gpdc: frequency, slowness

BASE = "BAS"                 # file prefix of every run
NPL = 3                      # parameters per node: depth, dVs, dVpVs (I_VPVS = 1)

# ---- prior / reference (constant reference => uniform absolute prior) ----
REF_VS, REF_VPVS = 2.1, 1.85           # km/s, -
DVS, DVPVS = 1.4, 0.15                 # one-sided widths: Vs 0.7-3.5, Vp/Vs 1.70-2.00
HMX, HMIN = 5.0, 0.02                  # km
NLMN, NLMX = 2, 30
TAIL_DEPTHS = (5.0, 8.0, 12.0)         # km, the "PREM" rows of vel_ref; first = HMX
SCAN = (0.6, 3.6, 0.005, 0.001)        # SWD_SCAN cmin cmax dc dc_over (km/s)
LAMBDA_K = 10.0                        # Poisson prior rate on the number of nodes

# ---- data ----
PERIODS_FULL = np.logspace(np.log10(0.5), np.log10(6.0), 30)
BANDS = {"full": (0.5, 6.0), "short": (0.5, 2.0), "long": (2.0, 6.0), "mismatch_g": (1.0, 4.0)}
NOISE = {0: 0.02, 1: 0.04}             # relative Gaussian noise by velocity type (0 phase, 1 group)
NOISE_SEED = 20260929


def legacy_rho(vp_kms):
    """The engine's density for the sampled layers (MAKE_CURMOD): g/cc."""
    vp = np.asarray(vp_kms, float)
    return 2.35 + 0.036 * (vp - 3.0) ** 2


TAIL_RHO = float(legacy_rho(3.3661 * 1.75))   # tail density column of vel_ref (fixed by the engine)


# --------------------------------------------------------------------------
# Michel et al. (2017) Basel-1 / Otterbach-2 model
# --------------------------------------------------------------------------
def load_michel(path=MICHEL_MODEL):
    """Literal Michel stack: (n,4) thick km, rho g/cc, Vp km/s, Vs km/s.

    The geopsy .model file holds thickness, Vp, Vs (m/s) and density (kg/m3;
    Gardner 310 Vp^0.25). Falls back to the Vs-only CSV with Vp = 1.9 Vs and
    Gardner density when the .model file is absent.
    """
    if os.path.exists(path):
        rows = np.loadtxt(path, skiprows=1)
        th, vp, vs, rho = rows[:, 0] / 1e3, rows[:, 1] / 1e3, rows[:, 2] / 1e3, rows[:, 3] / 1e3
    else:
        csv = np.loadtxt(MICHEL_CSV, delimiter=",", skiprows=1)
        vs_m, z_m = csv[:, 0], csv[:, 1]
        # rows come in pairs (top, bottom) of each layer
        vs = vs_m[0::2] / 1e3
        z = z_m[0::2] / 1e3
        th = np.r_[np.diff(z), 0.0]
        vp = 1.9 * vs
        rho = 0.31 * (vp * 1e3) ** 0.25
    th = th.copy()
    th[-1] = 0.0                                   # half-space
    return np.column_stack([th, rho, vp, vs])


def truth_nodes(stack=None):
    """The Michel model as ENGINE nodes: depth km, dVs, dVpVs (relative to the
    constant reference), one node per layer top (the first at 0)."""
    stack = load_michel() if stack is None else stack
    th, vp, vs = stack[:, 0], stack[:, 2], stack[:, 3]
    z = np.r_[0.0, np.cumsum(th[:-1])]
    dvs = vs - REF_VS
    dvpvs = vp / vs - REF_VPVS
    if np.any(np.abs(dvs) > DVS) or np.any(np.abs(dvpvs) > DVPVS) or z[-1] > HMX:
        raise ValueError("the Michel model lies outside the prior box; widen DVS/DVPVS/HMX")
    return z, dvs, dvpvs


# --------------------------------------------------------------------------
# Engine layer stack (MAKE_CURMOD + the tail of LOGLHOOD_SWD)
# --------------------------------------------------------------------------
def engine_stack(z, dvs, dvpvs, tail_depths=TAIL_DEPTHS):
    """Stack the engine builds from nodes (z ascending, all active):
    k sampled layers (the last runs to the first tail depth), then the tail
    rows shifted by the deepest node's dVs / dVpVs, half-space last.
    Vp = Vs * (REF_VPVS + dVpVs); rho = legacy relation; tail rho fixed."""
    z = np.asarray(z, float)
    order = np.argsort(z)
    z, dvs, dvpvs = z[order], np.asarray(dvs, float)[order], np.asarray(dvpvs, float)[order]
    vs = REF_VS + dvs
    vp = vs * (REF_VPVS + dvpvs)
    th = np.r_[np.diff(z), tail_depths[0] - z[-1]]
    rows = [th, legacy_rho(vp), vp, vs]
    # tail: vel_prem Vs (= REF_VS) + factvs, VpVs (= REF_VPVS) + factvpvs
    tvs = REF_VS + dvs[-1]
    tvp = tvs * (REF_VPVS + dvpvs[-1])
    tth = np.r_[np.diff(tail_depths), 0.0]
    tail = [tth, np.full(len(tail_depths), TAIL_RHO), np.full(len(tail_depths), tvp),
            np.full(len(tail_depths), tvs)]
    return np.column_stack([np.r_[a, b] for a, b in zip(rows, tail)])


def sample_nodes(row, nlmx=NLMX):
    """Nodes of one sample row: depth, dVs, dVpVs with inactive (-100)
    parameters carried down from the node above (INTERPLAYER: an inactive
    parameter merges the layer with the one above it)."""
    k = int(round(row[3]))
    v = row[4:4 + nlmx * NPL].reshape(nlmx, NPL)[:k].copy()
    v = v[v[:, 0] > -99.0]
    v = v[np.argsort(v[:, 0])]
    for j in (1, 2):
        for i in range(len(v)):
            if v[i, j] < -99.0:
                v[i, j] = v[i - 1, j] if i > 0 else 0.0
    return v[:, 0], v[:, 1], v[:, 2]


def vs_on_grid(z, dvs, zgrid):
    """Absolute Vs (km/s) of a node model on a depth grid (km)."""
    order = np.argsort(z)
    z, vs = np.asarray(z)[order], REF_VS + np.asarray(dvs)[order]
    lay = np.clip(np.searchsorted(z, zgrid, "right") - 1, 0, len(z) - 1)
    return vs[lay]


def vpvs_on_grid(z, dvpvs, zgrid):
    order = np.argsort(z)
    z, r = np.asarray(z)[order], REF_VPVS + np.asarray(dvpvs)[order]
    lay = np.clip(np.searchsorted(z, zgrid, "right") - 1, 0, len(z) - 1)
    return r[lay]


# --------------------------------------------------------------------------
# disp_driver (the production forward, phase AND group from one root search)
# --------------------------------------------------------------------------
def run_disp_driver(stacks, periods, maxmode=1, scan=SCAN, iwarm=0, driver=DISP_DRIVER):
    """Forward every stack at the (ascending) periods.

    Returns dict[(imodel, mode)] -> (valid bool array, c km/s, u km/s),
    imodel 0-based. Invalid periods carry c = u = 0.
    """
    periods = np.asarray(periods, float)
    if np.any(np.diff(periods) <= 0):
        raise ValueError("periods must be ascending")
    if not os.path.exists(driver):
        raise FileNotFoundError(f"{driver}: build it with the command in validation/README.md")
    with tempfile.TemporaryDirectory() as tmp:
        fmod = os.path.join(tmp, "models.txt")
        fper = os.path.join(tmp, "periods.txt")
        with open(fmod, "w") as fh:
            for st in stacks:
                fh.write(f"{len(st)}\n")
                for th, rho, vp, vs in st:
                    fh.write(f"{th:.8f} {rho:.6f} {vp:.7f} {vs:.7f}\n")
        np.savetxt(fper, periods, fmt="%.9f")
        cmd = [driver, fmod, fper, str(maxmode)] + [f"{x:g}" for x in scan] + [str(iwarm)]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines()[1:]:
        im, mode, ip, per, iv, c, u = line.split(",")
        key = (int(im) - 1, int(mode))
        if key not in res:
            n = len(periods)
            res[key] = (np.zeros(n, bool), np.zeros(n), np.zeros(n))
        i = int(ip) - 1
        res[key][0][i] = int(iv) == 1
        res[key][1][i] = float(c)
        res[key][2][i] = float(u)
    return res


# --------------------------------------------------------------------------
# Run-directory files
# --------------------------------------------------------------------------
def slot_filename(mode, grp, sd=False, base=BASE):
    """Engine naming (SWD_SLOT_FILE): _SWD / _SWDG [+ _M<m>], _sd... for sd files."""
    stem = ("_sdSWD" if sd else "_SWD") + ("G" if grp == 1 else "")
    return f"{base}{stem}.dat" if mode == 0 else f"{base}{stem}_M{mode}.dat"


def write_filebase(run_dir, base=BASE):
    with open(os.path.join(run_dir, "filebase.txt"), "w") as fh:
        fh.write(f"{len(base)}\n{base}\n")


def write_vel_ref(run_dir, base=BASE):
    """Constant reference (uniform absolute prior) + the tail rows."""
    ref = [(0.0, REF_VS, REF_VPVS, TAIL_RHO), (HMX, REF_VS, REF_VPVS, TAIL_RHO)]
    tail = [(d, REF_VS, REF_VPVS, TAIL_RHO) for d in TAIL_DEPTHS]
    with open(os.path.join(run_dir, f"{base}_vel_ref.txt"), "w") as fh:
        fh.write(f"{len(ref)}  {len(ref) + len(tail)}\n")
        for row in ref + tail:
            fh.write("   %.7e   %.7e   %.7e   %.7e\n" % row)


def write_map_voro(path, z, dvs, dvpvs, nmode, sdpar_swd=1.0, nmode_ell=1, nlmx=NLMX):
    """k, voro(NLMX*NPL) (unused slots 0), sdparSWD(NMODE), sdparELL, arparSWD(NMODE), arparELL."""
    k = len(z)
    voro = np.zeros((nlmx, NPL))
    order = np.argsort(z)
    voro[:k, 0] = np.asarray(z)[order]
    voro[:k, 1] = np.asarray(dvs)[order]
    voro[:k, 2] = np.asarray(dvpvs)[order]
    sd = np.full(nmode, float(sdpar_swd)) if np.isscalar(sdpar_swd) else np.asarray(sdpar_swd, float)
    row = np.concatenate([[k], voro.ravel(), sd, [1e-3] * nmode_ell, np.zeros(nmode), np.zeros(nmode_ell)])
    np.savetxt(path, row[None, :], fmt="%.10e")


def parameter_lines(nmode, mode_of, grp_of, ndat_swd, *, imap=0, icov_swd=3, imagscale=0,
                    isetseed=1, ismpprior=0, sdmn_swd=None, sdmx_swd=None, isd_swd=1,
                    sdmn=(0.3, 1e-3), sdmx=(3.0, 1e-1), nptchains1=6, nkeep=100, lam=LAMBDA_K,
                    extra_keywords=()):
    """The 48 positional lines + keyword tail of <base>_parameter.dat."""
    pos = [
        (imap, "IMAP"), (imagscale, "IMAGSCALE"), (0, "ENOS"), (1, "IPOIPR"), (0, "IAR"),
        (1, "I_VARPAR"), (0, "IBD_SINGLE"), (1, "I_SWD"), (0, "I_ELL"), (1, "I_VREF"),
        (1, "I_VPVS"), (ismpprior, "ISMPPRIOR"), (isetseed, "ISETSEED"), (1, "IEXCHANGE"),
        (ndat_swd, "NDAT_SWD"), (nmode, "NMODE"), (30, "NDAT_ELL"), (1, "NMODE_ELL"),
        (NLMN, "NLMN"), (NLMX, "NLMX"), (1, "ICHAINTHIN"), (nkeep, "NKEEP"), (nptchains1, "NPTCHAINS1"),
        (1.15, "dTlog"), (lam, "lambda"), (HMX, "hmx"), (HMIN, "hmin"), (0.1, "armxSWD"),
        (0.1, "armxELL"), (1000, "TCHCKPT"), (DVS, "dVs"), (DVPVS, "dVpVs"),
        (f"{sdmn[0]:g} {sdmn[1]:g}", "sdmn"), (f"{sdmx[0]:g} {sdmx[1]:g}", "sdmx"),
        (isd_swd, "ISD_SWD"), (0, "ISD_ELL"), (icov_swd, "ICOV_SWD"), (1, "ICOV_ELL"),
        (0, "ELL_verbose"), ("1.e-5", "ELL_prec"), (1, "I_ABS_ELL"), (0, "I_LOG10_ELL"),
        (0, "I_SAMPLING_TYPE_ELL"), (0, "I_SET_STEP_ELL"), (1.5, "STEP_SIZE_ELL"),
        (0, "I_SET_COUNT_ELL"), (145, "COUNT_ELL"), (1, "I_SET_RANGE_ELL"),
    ]
    assert len(pos) == 48
    lines = [f"{str(v):<18} !! {i + 1:2d} {name}" for i, (v, name) in enumerate(pos)]
    lines.append("!! ---- keyword tail ----")
    lines.append("MODE_OF " + " ".join(str(m) for m in mode_of))
    lines.append("GRP_OF  " + " ".join(str(g) for g in grp_of))
    if sdmn_swd is not None:
        lines.append("SDMN_SWD " + " ".join(f"{s:g}" for s in sdmn_swd))
    if sdmx_swd is not None:
        lines.append("SDMX_SWD " + " ".join(f"{s:g}" for s in sdmx_swd))
    lines.append("SWD_SCAN " + " ".join(f"{x:g}" for x in SCAN))
    lines.append("SWD_WARM 0     !! the Michel model has LVZs: no DVSCON/DVSMONO, cold scans")
    lines += list(extra_keywords)
    return "\n".join(lines) + "\n"


def write_covparameter(run_dir, base=BASE):
    """Covariance iteration OFF (single-pass rjMcMC); other lines are inert."""
    lines = [
        (0, "Icov_iterUpdate_SWD"), (0, "Icov_iterUpdate_ELL"), (2000, "covIter_zero_nsamples"),
        (2000, "covIter_period"), (100, "MAXcovIter"), (0, "ICOVest"),
        (5, "CHAINTHIN_COVest_period_zeroIter"), (4, "CHAINTHIN_COVest_period_nonzeroIter"),
        (1, "ISD_SWD_covIter"), (0, "ISD_ELL_covIter"), ("1. 1.", "sdmn_covIter (SWD ELL)"),
        ("10. 10.", "sdmx_covIter (SWD ELL)"), ("2. 2.", "sdpar_covIter (SWD ELL)"),
        (100, "NKEEP_covIter"), (100, "NKEEP_covIter_res"), (1, "iSAVEsample_covIter"),
        (0, "iSAVEsample_only_zeroIter"), (0, "iMAP_calc"), (1, "iconverge_criterion"),
        (2, "iconverge_criterion_SWD"), (2, "iconverge_criterion_ELL"),
        ("1.e-4", "converge_threshold_SWD"), ("1.e-4", "converge_threshold_ELL"),
        (10, "nfrac_SWD"), (40, "MAX_NAVE_SWD"), (1, "inonstat_SWD"), (0, "iunbiased_SWD"),
        (0, "imr_SWD"), (1.2, "damp_power_SWD"), (10, "nfrac_ELL"), (40, "MAX_NAVE_ELL"),
        (1, "inonstat_ELL"), (0, "iunbiased_ELL"), (0, "imr_ELL"), (1.2, "damp_power_ELL"),
    ]
    assert len(lines) == 35
    with open(os.path.join(run_dir, f"{base}_covparameter.dat"), "w") as fh:
        for v, name in lines:
            fh.write(f"{str(v):<14} !!{name}\n")


def write_curve(path, periods, values):
    np.savetxt(path, np.column_stack([periods, values]), fmt="%.9e")


def read_cfg(run_dir, base=BASE):
    """Settings of a run dir needed by the analysis (positional + keyword)."""
    cfg = {}
    with open(os.path.join(run_dir, f"{base}_parameter.dat")) as fh:
        lines = [l.rstrip("\n") for l in fh]
    val = lambda i: lines[i - 1].split("!!")[0].split()
    cfg["IMAGSCALE"] = int(val(2)[0]); cfg["NDAT_SWD"] = int(val(15)[0]); cfg["NMODE"] = int(val(16)[0])
    cfg["NMODE_ELL"] = int(val(18)[0]); cfg["NLMX"] = int(val(20)[0]); cfg["NPTCHAINS1"] = int(val(23)[0])
    cfg["lambda"] = float(val(25)[0]); cfg["hmx"] = float(val(26)[0]); cfg["ICOV_SWD"] = int(val(37)[0])
    cfg["sdmn"] = [float(val(33)[0])] * cfg["NMODE"]; cfg["sdmx"] = [float(val(34)[0])] * cfg["NMODE"]
    cfg["MODE_OF"] = list(range(cfg["NMODE"])); cfg["GRP_OF"] = [0] * cfg["NMODE"]
    for l in lines[48:]:
        l = l.split("!")[0].split()
        if not l:
            continue
        kw, rest = l[0].upper(), l[1:]
        if kw == "MODE_OF":
            cfg["MODE_OF"] = [int(x) for x in rest]
        elif kw == "GRP_OF":
            cfg["GRP_OF"] = [int(x) for x in rest]
        elif kw == "IGRP":
            cfg["GRP_OF"] = [int(rest[0])] * cfg["NMODE"]
        elif kw == "SDMN_SWD":
            cfg["sdmn"] = [float(x) for x in rest]
        elif kw == "SDMX_SWD":
            cfg["sdmx"] = [float(x) for x in rest]
    cfg["slots"] = list(zip(cfg["MODE_OF"], cfg["GRP_OF"]))
    return cfg


def read_slot_data(run_dir, cfg, base=BASE):
    """Observed curves and per-point sd of every slot: list of dicts."""
    out = []
    for mode, grp in cfg["slots"]:
        d = np.loadtxt(os.path.join(run_dir, slot_filename(mode, grp, base=base)))
        sd_path = os.path.join(run_dir, slot_filename(mode, grp, sd=True, base=base))
        sd = np.loadtxt(sd_path) if os.path.exists(sd_path) else None
        out.append({"mode": mode, "grp": grp, "periods": d[:, 0], "obs": d[:, 1], "sd": sd})
    return out


def read_sample(run_dir, base=BASE):
    """Sample matrix. Rows whose logL is -HUGE (a start state that could not
    predict an observed mode) are written by the engine with a 3-digit
    exponent that drops the 'E' ("-1.79769313+308"); they are unreadable for
    numpy and useless for the posterior, so they are removed."""
    path = os.path.join(run_dir, f"{base}_voro_sample.txt")
    try:
        dat = np.loadtxt(path)
    except ValueError:
        import re
        bad = re.compile(r"\d[+-]\d{3}\b")
        with open(path) as fh:
            lines = [l for l in fh if not bad.search(l)]
        dat = np.loadtxt(lines)
    dat = dat[None, :] if dat.ndim == 1 else dat
    return dat[dat[:, 0] > -1e300]


def sample_layout(cfg):
    """Column offsets of the sample row: logL, logPr, tcmp, k, voro, sdparSWD,
    sdparELL, arparSWD, arparELL, acc, iacc_bd, irej_bd, iacc_bds, chain, source."""
    isig = 4 + cfg["NLMX"] * NPL
    return {"isig": isig, "ncol": isig + 2 * cfg["NMODE"] + 2 * cfg["NMODE_ELL"] + 6}


def burn_in_split(dat, burn_frac):
    """Drop the first burn_frac of every chain (grouped by the source column)."""
    ranks = dat[:, -1].astype(int)
    keep = np.zeros(len(dat), bool)
    for r in np.unique(ranks):
        idx = np.where(ranks == r)[0]
        keep[idx[int(len(idx) * burn_frac):]] = True
    return dat[keep], dat[~keep]


def run_engine_imap(run_dir, np_ranks=4, engine=ENGINE, timeout=600):
    """Run the engine in IMAP mode in run_dir (its parameter file must have
    IMAP 1). Returns (logL, stdout)."""
    out = subprocess.run(["mpirun", "-np", str(np_ranks), engine], cwd=run_dir,
                         capture_output=True, text=True, timeout=timeout)
    logl = None
    for line in out.stdout.splitlines():
        if line.strip().startswith("logL =") or line.strip().startswith("logL ="):
            logl = float(line.split("=")[1])
    if logl is None:
        raise RuntimeError(f"IMAP run in {run_dir} printed no logL:\n{out.stdout[-3000:]}\n{out.stderr[-2000:]}")
    return logl, out.stdout


def read_mappred(run_dir, base=BASE):
    """Rows of <base>_mappredSWD.dat (one per slot, NDAT_SWD wide, zero padded)."""
    rows = np.loadtxt(os.path.join(run_dir, f"{base}_mappredSWD.dat"), ndmin=2)
    return rows


def copy_run_inputs(src, dst, base=BASE):
    os.makedirs(dst, exist_ok=True)
    for f in os.listdir(src):
        if f == "filebase.txt" or (f.startswith(base + "_") and not f.startswith(base + "_voro_sample")
                                   and "mappred" not in f and "_obs" not in f and "_mapar" not in f):
            shutil.copy(os.path.join(src, f), os.path.join(dst, f))
