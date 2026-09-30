"""Synthetic-test-specific pieces of the Basel-1 / Otterbach-2 group-velocity
harness (validation/basel_group): the prior / reference constants, the Michel
model, the run-directory writers. Everything that reads engine files or
rebuilds the engine's layer stack comes from tools/swdell (the general
post-processing package); the thin aliases below keep the harness scripts
short.
"""
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "tools"))
from swdell import forward, io, model  # noqa: E402

ENGINE = os.path.join(REPO, "src", "bin", "prjmh_temper_rf")
DISP_DRIVER = forward.DISP_DRIVER

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

legacy_rho = model.legacy_rho
TAIL_RHO = float(legacy_rho(3.3661 * 1.75))   # tail density column of vel_ref (fixed by the engine)

# the engine configuration every staged run shares (what parameter_lines writes)
CFG = {"NLMX": NLMX, "NLMN": NLMN, "NPL": NPL, "I_VPVS": 1, "VP_BROCHER": 0, "RHO_BROCHER": 0, "hmx": HMX}
VEL_REF = np.array([(0.0, REF_VS, REF_VPVS, TAIL_RHO), (HMX, REF_VS, REF_VPVS, TAIL_RHO)])
VEL_PREM = np.array([(d, REF_VS, REF_VPVS, TAIL_RHO) for d in TAIL_DEPTHS])


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
        vs = vs_m[0::2] / 1e3                      # rows come in pairs (top, bottom) of each layer
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


def engine_stack(z, dvs, dvpvs):
    """Stack the engine builds from nodes under the harness configuration
    (swdell.model.engine_stack with the constant reference and the tail rows)."""
    return model.engine_stack(z, dvs, dvpvs, CFG, VEL_REF, VEL_PREM)


def run_disp_driver(stacks, periods, maxmode=1, scan=SCAN, iwarm=0, driver=None):
    return forward.run_disp_driver(stacks, periods, maxmode=maxmode, scan=scan, iwarm=iwarm, driver=driver)


# --------------------------------------------------------------------------
# Run-directory files
# --------------------------------------------------------------------------
def slot_filename(mode, grp, sd=False, base=BASE):
    return io.slot_filename(base, mode, grp, sd=sd)


def write_filebase(run_dir, base=BASE):
    with open(os.path.join(run_dir, "filebase.txt"), "w") as fh:
        fh.write(f"{len(base)}\n{base}\n")


def write_vel_ref(run_dir, base=BASE):
    """Constant reference (uniform absolute prior) + the tail rows."""
    with open(os.path.join(run_dir, f"{base}_vel_ref.txt"), "w") as fh:
        fh.write(f"{len(VEL_REF)}  {len(VEL_REF) + len(VEL_PREM)}\n")
        for row in np.vstack([VEL_REF, VEL_PREM]):
            fh.write("   %.7e   %.7e   %.7e   %.7e\n" % tuple(row))


def write_map_voro(path, z, dvs, dvpvs, nmode, sdpar_swd=1.0, nmode_ell=1, nlmx=NLMX):
    """k, voro(NLMX*NPL) (unused slots 0), sdparSWD(NMODE), sdparELL, arparSWD(NMODE), arparELL."""
    cfg = {"NLMX": nlmx, "NPL": NPL, "NMODE": nmode, "NMODE_ELL": nmode_ell}
    io.write_map_voro(path, cfg, z, dvs, dvpvs, sdpar_swd=sdpar_swd, sdpar_ell=1e-3)


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


def run_engine_imap(run_dir, np_ranks=4, engine=ENGINE, timeout=600):
    """Run the engine in IMAP mode in run_dir (its parameter file must have
    IMAP 1). Returns (logL, stdout)."""
    out = subprocess.run(["mpirun", "-np", str(np_ranks), engine], cwd=run_dir,
                         capture_output=True, text=True, timeout=timeout)
    logl = None
    for line in out.stdout.splitlines():
        if line.strip().startswith("logL ="):
            logl = float(line.split("=")[1])
    if logl is None:
        raise RuntimeError(f"IMAP run in {run_dir} printed no logL:\n{out.stdout[-3000:]}\n{out.stderr[-2000:]}")
    return logl, out.stdout


def read_mappred(run_dir, base=BASE):
    """Rows of <base>_mappredSWD.dat (one per slot, NDAT_SWD wide, zero padded)."""
    return np.loadtxt(os.path.join(run_dir, f"{base}_mappredSWD.dat"), ndmin=2)


def copy_run_inputs(src, dst, base=BASE):
    os.makedirs(dst, exist_ok=True)
    for f in os.listdir(src):
        if f == "filebase.txt" or (f.startswith(base + "_") and not f.startswith(base + "_voro_sample")
                                   and "mappred" not in f and "_obs" not in f and "_mapar" not in f):
            shutil.copy(os.path.join(src, f), os.path.join(dst, f))
