"""Forward modelling through validation/disp_driver (the production DISPER80
code, dispersion_cu: phase AND group velocity from one root search)."""
import os
import subprocess
import tempfile

import numpy as np

from . import model

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DISP_DRIVER = os.environ.get("SWDELL_DISP_DRIVER", os.path.join(REPO, "validation", "disp_driver"))


def available(driver=None):
    return os.path.exists(driver or DISP_DRIVER)


def run_disp_driver(stacks, periods, maxmode=0, scan=(2.0, 6.5, 0.05, -1.0), iwarm=0, driver=None):
    """Forward every stack at the ascending periods.
    Returns dict[(imodel, mode)] -> (valid bool array, c km/s, u km/s)."""
    driver = driver or DISP_DRIVER
    periods = np.asarray(periods, float)
    if np.any(np.diff(periods) <= 0):
        raise ValueError("periods must be ascending")
    if not os.path.exists(driver):
        raise FileNotFoundError(f"{driver}: build it with the command in validation/README.md")
    cmin, cmax, dc, dc_over = scan
    if dc_over <= 0:
        dc_over = dc / 5.0
    with tempfile.TemporaryDirectory() as tmp:
        fmod = os.path.join(tmp, "models.txt")
        fper = os.path.join(tmp, "periods.txt")
        with open(fmod, "w") as fh:
            for st in stacks:
                fh.write(f"{len(st)}\n")
                for th, rho, vp, vs in st:
                    fh.write(f"{th:.8f} {rho:.6f} {vp:.7f} {vs:.7f}\n")
        np.savetxt(fper, periods, fmt="%.9f")
        cmd = [driver, fmod, fper, str(maxmode), f"{cmin:g}", f"{cmax:g}", f"{dc:g}", f"{dc_over:g}", str(iwarm)]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    res = {}
    n = len(periods)
    for line in out.splitlines()[1:]:
        im, mode, ip, per, iv, c, u = line.split(",")
        key = (int(im) - 1, int(mode))
        if key not in res:
            res[key] = (np.zeros(n, bool), np.zeros(n), np.zeros(n))
        i = int(ip) - 1
        res[key][0][i] = int(iv) == 1
        res[key][1][i] = float(c)
        res[key][2][i] = float(u)
    return res


def predict_slots(run, rows, nfwd=300, rng=None, driver=None):
    """Posterior-predictive curves of every slot for a subsample of rows:
    list of arrays (nfwd, npts), NaN where the mode does not exist.
    One disp_driver call per distinct period grid, with the run's SWD_SCAN."""
    rng = rng or np.random.default_rng(0)
    if len(rows) > nfwd:
        rows = rows[rng.choice(len(rows), nfwd, replace=False)]
    cfg = run.cfg
    stacks = [model.engine_stack(*model.sample_nodes(r, cfg), cfg, run.vel_ref, run.vel_prem) for r in rows]
    grids = {}
    for i, s in enumerate(run.slots):
        grids.setdefault(tuple(np.round(s["periods"], 9)), []).append(i)
    pred = [np.full((len(rows), len(s["periods"])), np.nan) for s in run.slots]
    for key, islots in grids.items():
        periods = np.array(key)
        maxmode = max(cfg["slots"][i][0] for i in islots)
        res = run_disp_driver(stacks, periods, maxmode=maxmode, scan=cfg["SWD_SCAN"], driver=driver)
        for i in islots:
            mode, grp = cfg["slots"][i]
            for im in range(len(stacks)):
                valid, c, u = res[(im, mode)]
                pred[i][im] = np.where(valid, u if grp == 1 else c, np.nan)
    return pred
