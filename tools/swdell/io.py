"""Readers and writers for the engine's input and output files.

Everything here follows src/read_input.f90 and src/prjmh_temper_rf.f90:
- <base>_parameter.dat: 48 positional lines, then KEYWORD lines (unknown
  lines ignored; GRP_OF wins over IGRP whatever their order; SDMN_SWD /
  SDMX_SWD default to the scalar sdmn(1) / sdmx(1)).
- curve slots identified by (MODE_OF, GRP_OF); file names by SWD_SLOT_FILE.
- <base>_voro_sample.txt rows: logL, logPr, tcmp, k, voro(NLMX*NPL),
  sdparSWD(NMODE), sdparELL(NMODE_ELL), arparSWD(NMODE), arparELL(NMODE_ELL),
  acc, iaccept_bd, ireject_bd, iaccept_bds, chain, source_rank.
- <base>_map_voro.dat: k, voro(NLMX*NPL), sdparSWD, sdparELL, arparSWD, arparELL.
- <base>_mappredSWD.dat / _obsSWD.dat / _maparSWD.dat: one row per slot,
  NDAT_SWD wide, zero padded (SAVEREPLICA, IMAP mode).
- <base>_vel_ref.txt: "NVELREF NTOT", then rows depth Vs VpVs rho; the rows
  beyond NVELREF are the fixed deep tail ("PREM").
"""
import os
import re
from dataclasses import dataclass, field

import numpy as np

KEYWORDS = ("DVSCON", "DVSMONO", "VP_BROCHER", "RHO_BROCHER", "MODE_OF", "GRP_OF", "IGRP",
            "SDMN_SWD", "SDMX_SWD", "SWD_SCAN", "SWD_WARM")

# the 48 positional lines of READPARFILE, in order
POSITIONAL = ("IMAP", "IMAGSCALE", "ENOS", "IPOIPR", "IAR", "I_VARPAR", "IBD_SINGLE", "I_SWD", "I_ELL",
              "I_VREF", "I_VPVS", "ISMPPRIOR", "ISETSEED", "IEXCHANGE", "NDAT_SWD", "NMODE", "NDAT_ELL",
              "NMODE_ELL", "NLMN", "NLMX", "ICHAINTHIN", "NKEEP", "NPTCHAINS1", "dTlog", "lambda", "hmx",
              "hmin", "armxSWD", "armxELL", "TCHCKPT", "dVs", "dVpVs", "sdmn", "sdmx", "ISD_SWD", "ISD_ELL",
              "ICOV_SWD", "ICOV_ELL", "ELL_verbose", "ELL_prec", "I_ABS_ELL", "I_LOG10_ELL",
              "I_SAMPLING_TYPE_ELL", "I_SET_STEP_ELL", "STEP_SIZE_ELL", "I_SET_COUNT_ELL", "COUNT_ELL",
              "I_SET_RANGE_ELL")
INT_FIELDS = {"IMAP", "IMAGSCALE", "ENOS", "IPOIPR", "IAR", "I_VARPAR", "IBD_SINGLE", "I_SWD", "I_ELL",
              "I_VREF", "I_VPVS", "ISMPPRIOR", "ISETSEED", "IEXCHANGE", "NDAT_SWD", "NMODE", "NDAT_ELL",
              "NMODE_ELL", "NLMN", "NLMX", "ICHAINTHIN", "NKEEP", "NPTCHAINS1", "ISD_SWD", "ISD_ELL",
              "ICOV_SWD", "ICOV_ELL", "ELL_verbose", "I_ABS_ELL", "I_LOG10_ELL", "I_SAMPLING_TYPE_ELL",
              "I_SET_STEP_ELL", "I_SET_COUNT_ELL", "COUNT_ELL", "I_SET_RANGE_ELL"}


def find_base(run_dir):
    """The file prefix from filebase.txt (length line, then the name)."""
    with open(os.path.join(run_dir, "filebase.txt")) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    return lines[-1]


def read_parfile(path):
    """<base>_parameter.dat -> dict of the 48 positional values + keywords +
    derived fields NPL, slots (list of (mode, grp)), sdmn_swd, sdmx_swd."""
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    if len(lines) < 48:
        raise ValueError(f"{path}: {len(lines)} lines, the parameter file has 48 positional lines")
    cfg = {}
    for i, name in enumerate(POSITIONAL):
        tok = lines[i].split("!!")[0].split()
        if name in ("sdmn", "sdmx"):
            cfg[name] = [float(tok[0]), float(tok[1])]
        elif name in INT_FIELDS:
            cfg[name] = int(float(tok[0]))
        else:
            cfg[name] = float(tok[0])
    nmode = cfg["NMODE"]
    cfg.update({"DVSCON": -1.0, "DVSMONO": -1.0, "VP_BROCHER": 0, "RHO_BROCHER": 0, "IGRP": 0,
                "SWD_SCAN": [2.0, 6.5, 0.05, -1.0], "SWD_WARM": -1,
                "MODE_OF": list(range(nmode)), "GRP_OF": None, "SDMN_SWD": None, "SDMX_SWD": None})
    for l in lines[48:]:
        tok = l.split("!")[0].split()
        if not tok:
            continue
        kw, rest = tok[0].upper(), tok[1:]
        if kw in ("DVSCON", "DVSMONO"):
            cfg[kw] = float(rest[0])
        elif kw in ("VP_BROCHER", "RHO_BROCHER", "IGRP", "SWD_WARM"):
            cfg[kw] = int(rest[0])
        elif kw in ("MODE_OF", "GRP_OF"):
            cfg[kw] = [int(x) for x in rest[:nmode]]
        elif kw in ("SDMN_SWD", "SDMX_SWD"):
            cfg[kw] = [float(x) for x in rest[:nmode]]
        elif kw == "SWD_SCAN":
            v = [float(x) for x in rest[:4]]
            cfg[kw] = v + [-1.0] * (4 - len(v))
    # resolve the defaults exactly as READPARFILE does after its keyword loop
    if cfg["GRP_OF"] is None:
        cfg["GRP_OF"] = [cfg["IGRP"]] * nmode
    if cfg["SDMN_SWD"] is None:
        cfg["SDMN_SWD"] = [cfg["sdmn"][0]] * nmode
    if cfg["SDMX_SWD"] is None:
        cfg["SDMX_SWD"] = [cfg["sdmx"][0]] * nmode
    cfg["NPL"] = 3 if cfg["I_VPVS"] == 1 else 2
    cfg["slots"] = list(zip(cfg["MODE_OF"], cfg["GRP_OF"]))
    return cfg


def slot_filename(base, mode, grp, sd=False):
    """SWD_SLOT_FILE: _SWD / _SWDG (+ _M<m>), _sd... for the sd files."""
    stem = ("_sdSWD" if sd else "_SWD") + ("G" if grp == 1 else "")
    return f"{base}{stem}.dat" if mode == 0 else f"{base}{stem}_M{mode}.dat"


def slot_label(mode, grp, short=False):
    return f"R{mode}{'g' if grp else 'p'}" if short else f"R{mode} {'group' if grp else 'phase'}"


def read_vel_ref(path):
    """-> (vel_ref (NVELREF,4), vel_prem (NPREM,4)); columns depth km, Vs, VpVs, rho."""
    with open(path) as fh:
        head = fh.readline().split()
        nref, ntot = int(head[0]), int(head[1])
        rows = np.array([[float(x) for x in fh.readline().split()] for _ in range(ntot)])
    return rows[:nref], rows[nref:]


def read_slots(run_dir, cfg, base):
    """Observed curves per slot: list of dicts (mode, grp, periods, obs, sd or None)."""
    out = []
    for mode, grp in cfg["slots"]:
        d = np.loadtxt(os.path.join(run_dir, slot_filename(base, mode, grp)), ndmin=2)
        sd = None
        if cfg["ICOV_SWD"] == 3:
            sd = np.loadtxt(os.path.join(run_dir, slot_filename(base, mode, grp, sd=True)))
        out.append({"mode": mode, "grp": grp, "periods": d[:, 0], "obs": d[:, 1], "sd": sd})
    return out


_BAD_TOKEN = re.compile(r"\d[+-]\d{3}\b")   # "-1.79769313+308": -HUGE printed without the E


def read_sample(path):
    """Sample matrix. Rows with logL = -HUGE (a start state that could not
    predict an observed mode; written with a 3-digit exponent that drops the
    'E') are unreadable and useless: they are removed."""
    try:
        dat = np.loadtxt(path)
    except ValueError:
        with open(path) as fh:
            dat = np.loadtxt([l for l in fh if not _BAD_TOKEN.search(l)])
    dat = dat[None, :] if dat.ndim == 1 else dat
    return dat[dat[:, 0] > -1e300]


def sample_layout(cfg):
    """Column offsets of a sample row."""
    isig = 4 + cfg["NLMX"] * cfg["NPL"]
    nm, ne = cfg["NMODE"], cfg["NMODE_ELL"]
    return {"k": 3, "voro": (4, isig), "sdparSWD": (isig, isig + nm), "sdparELL": (isig + nm, isig + nm + ne),
            "arparSWD": (isig + nm + ne, isig + 2 * nm + ne), "arparELL": (isig + 2 * nm + ne, isig + 2 * nm + 2 * ne),
            "acc": isig + 2 * nm + 2 * ne, "chain": -2, "source": -1, "ncol": isig + 2 * nm + 2 * ne + 6}


def burn_in_split(dat, burn_frac):
    """Drop the first burn_frac of every chain (grouped by the source column)."""
    ranks = dat[:, -1].astype(int)
    keep = np.zeros(len(dat), bool)
    for r in np.unique(ranks):
        idx = np.where(ranks == r)[0]
        keep[idx[int(len(idx) * burn_frac):]] = True
    return dat[keep], dat[~keep]


def read_map_voro(path, cfg):
    """-> dict k, voro (NLMX, NPL), sdparSWD, sdparELL, arparSWD, arparELL."""
    v = np.loadtxt(path)
    n = cfg["NLMX"] * cfg["NPL"]
    nm, ne = cfg["NMODE"], cfg["NMODE_ELL"]
    i = 1 + n
    return {"k": int(round(v[0])), "voro": v[1:i].reshape(cfg["NLMX"], cfg["NPL"]),
            "sdparSWD": v[i:i + nm], "sdparELL": v[i + nm:i + nm + ne],
            "arparSWD": v[i + nm + ne:i + 2 * nm + ne], "arparELL": v[i + 2 * nm + ne:i + 2 * nm + 2 * ne]}


def write_map_voro(path, cfg, z, dvs, dvpvs=None, sdpar_swd=1.0, sdpar_ell=1e-3, arpar_swd=0.0, arpar_ell=0.0):
    """k, voro(NLMX*NPL) (unused slots 0), sdparSWD(NMODE), sdparELL, arparSWD(NMODE), arparELL."""
    nlmx, npl, nm, ne = cfg["NLMX"], cfg["NPL"], cfg["NMODE"], cfg["NMODE_ELL"]
    k = len(z)
    voro = np.zeros((nlmx, npl))
    order = np.argsort(z)
    voro[:k, 0] = np.asarray(z)[order]
    voro[:k, 1] = np.asarray(dvs)[order]
    if npl == 3:
        voro[:k, 2] = np.asarray(dvpvs if dvpvs is not None else np.zeros(k))[order]
    expand = lambda v, n: np.full(n, float(v)) if np.isscalar(v) else np.asarray(v, float)
    row = np.concatenate([[k], voro.ravel(), expand(sdpar_swd, nm), expand(sdpar_ell, ne),
                          expand(arpar_swd, nm), expand(arpar_ell, ne)])
    np.savetxt(path, row[None, :], fmt="%.10e")


def sample_row_to_map(row, cfg):
    """A sample row -> the map_voro vector (k, voro, sdpar/arpar), as UPDATE_MAPfile
    does (columns 4 .. end-6 of the sample row)."""
    lay = sample_layout(cfg)
    return row[lay["k"]:lay["acc"]]


def read_pred_rows(path, cfg, slots):
    """Rows of <base>_mappredSWD.dat / _obsSWD.dat / _maparSWD.dat trimmed to
    each slot's point count -> list of arrays (or None if the file is absent)."""
    if not os.path.exists(path):
        return None
    rows = np.loadtxt(path, ndmin=2)
    return [rows[i, :len(s["periods"])] for i, s in enumerate(slots)]


@dataclass
class Run:
    run_dir: str
    base: str
    cfg: dict
    vel_ref: np.ndarray
    vel_prem: np.ndarray
    slots: list
    sample: np.ndarray = None
    layout: dict = field(default_factory=dict)
    map_voro: dict = None
    mappred: list = None      # MAP prediction rows per slot (IMAP output), if present
    mapobs: list = None

    @property
    def name(self):
        return os.path.basename(os.path.abspath(self.run_dir))

    def path(self, suffix):
        return os.path.join(self.run_dir, self.base + suffix)


def read_run(run_dir, with_sample=True):
    """Everything the plots need from a run directory."""
    base = find_base(run_dir)
    cfg = read_parfile(os.path.join(run_dir, f"{base}_parameter.dat"))
    if cfg["I_SWD"] != 1:
        raise ValueError(f"{run_dir}: I_SWD = {cfg['I_SWD']}; only SWD runs are supported")
    if cfg["I_VPVS"] not in (1, -1):
        raise ValueError(f"{run_dir}: I_VPVS = {cfg['I_VPVS']} (0 = sample Vp) is not supported")
    if cfg["I_VREF"] != 1:
        raise ValueError(f"{run_dir}: I_VREF = 0 (no vel_ref) is not supported")
    vel_ref, vel_prem = read_vel_ref(os.path.join(run_dir, f"{base}_vel_ref.txt"))
    slots = read_slots(run_dir, cfg, base)
    run = Run(run_dir=run_dir, base=base, cfg=cfg, vel_ref=vel_ref, vel_prem=vel_prem, slots=slots,
              layout=sample_layout(cfg))
    mp = os.path.join(run_dir, f"{base}_map_voro.dat")
    if os.path.exists(mp):
        run.map_voro = read_map_voro(mp, cfg)
    run.mappred = read_pred_rows(os.path.join(run_dir, f"{base}_mappredSWD.dat"), cfg, slots)
    run.mapobs = read_pred_rows(os.path.join(run_dir, f"{base}_obsSWD.dat"), cfg, slots)
    sp = os.path.join(run_dir, f"{base}_voro_sample.txt")
    if with_sample and os.path.exists(sp):
        run.sample = read_sample(sp)
        if run.sample.shape[1] != run.layout["ncol"]:
            raise ValueError(f"{sp}: {run.sample.shape[1]} columns, the layout expects {run.layout['ncol']} "
                             f"(NLMX {cfg['NLMX']}, NPL {cfg['NPL']}, NMODE {cfg['NMODE']}, NMODE_ELL {cfg['NMODE_ELL']})")
    return run
