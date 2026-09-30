"""Posterior summaries: depth profiles, node densities, convergence metrics."""
from math import lgamma

import numpy as np

from . import model


def depth_grid(run, dz=0.01, zmax=None):
    zmax = run.cfg["hmx"] if zmax is None else zmax
    return np.arange(0.0, zmax + 1e-9, dz)


def profiles(rows, run, zgrid, nmax=4000, rng=None):
    """Vs(z), VpVs(z) per sample on zgrid and all node depths (> 0)."""
    rng = rng or np.random.default_rng(0)
    if nmax and len(rows) > nmax:
        rows = rows[rng.choice(len(rows), nmax, replace=False)]
    cfg, vr = run.cfg, run.vel_ref
    vs = np.empty((len(rows), len(zgrid)))
    vpvs = np.empty_like(vs)
    zall = []
    for i, r in enumerate(rows):
        z, dvs, dvpvs = model.sample_nodes(r, cfg)
        vs[i] = model.vs_on_grid(z, dvs, zgrid, cfg, vr)
        vpvs[i] = model.vpvs_on_grid(z, dvpvs, zgrid, cfg, vr)
        zall.append(z[z > 0])
    return vs, vpvs, (np.concatenate(zall) if zall else np.zeros(0))


def quantiles(a, q=(16, 50, 84)):
    return np.percentile(a, list(q), axis=0)


def per_depth_hist(vals, bins):
    """Normalised (per depth, max = 1) histogram, rows = depths."""
    h = np.array([np.histogram(vals[:, i], bins)[0] for i in range(vals.shape[1])], float)
    return h / np.maximum(h.max(axis=1, keepdims=True), 1)


def interface_density(znodes, zbins):
    h, _ = np.histogram(znodes, zbins)
    return h / max(h.sum() * (zbins[1] - zbins[0]), 1)


def poisson_prior(cfg):
    kk = np.arange(cfg["NLMN"], cfg["NLMX"] + 1)
    lam = cfg["lambda"]
    pk = np.exp(-lam + kk * np.log(lam) - np.array([lgamma(k + 1) for k in kk]))
    return kk, pk / pk.sum()


def l1_halfwidth(qa, qb, ref):
    """Depth-mean |median difference| in units of ref's 16-84 half-width."""
    half = np.maximum((ref[2] - ref[0]) / 2.0, 1e-4)
    return float(np.mean(np.abs(qa[1] - qb[1]) / half))


def probe_vs(rows, run, depths):
    zp = np.asarray(depths, float)
    out = np.empty((len(rows), len(zp)))
    for i, r in enumerate(rows):
        z, dvs, _ = model.sample_nodes(r, run.cfg)
        out[i] = model.vs_on_grid(z, dvs, zp, run.cfg, run.vel_ref)
    return out


def running_median(x, npts=80):
    idx = np.unique(np.linspace(1, len(x), npts).astype(int))
    return idx, np.array([np.median(x[:i], axis=0) for i in idx])


def convergence_metrics(dat, post, q, run, zgrid, burn, nprof=4000, rng=None):
    """The three L1 metrics of masw-das/rjmcmc_convergence.py (in 68 %
    half-width units of the full post-burn-in posterior q):
    stationarity (first vs second half), between-chain (worst pair of source
    ranks), burn-in 30 % vs 60 %. <= 0.05 comfortable, >= 0.15 re-run."""
    from .io import burn_in_split
    rng = rng or np.random.default_rng(0)
    chains = np.unique(dat[:, -1]).astype(int)
    q1 = quantiles(profiles(post[: len(post) // 2], run, zgrid, nprof, rng)[0])
    q2 = quantiles(profiles(post[len(post) // 2:], run, zgrid, nprof, rng)[0])
    m_stat = l1_halfwidth(q1, q2, q)
    qc = [quantiles(profiles(post[post[:, -1] == c], run, zgrid, nprof // 2, rng)[0]) for c in chains]
    m_chain = max((l1_halfwidth(qc[a], qc[b], q) for a in range(len(qc)) for b in range(a + 1, len(qc))), default=0.0)
    q60 = quantiles(profiles(burn_in_split(dat, 0.6)[0], run, zgrid, nprof, rng)[0])
    q30 = q if abs(burn - 0.3) < 1e-9 else quantiles(profiles(burn_in_split(dat, 0.3)[0], run, zgrid, nprof, rng)[0])
    m_burn = l1_halfwidth(q30, q60, q)
    return {"stationarity_L1": m_stat, "between_chain_L1": m_chain, "burn30v60_L1": m_burn}


def move_fraction(post, run):
    """Fraction of consecutive samples of a chain that differ in the model
    columns (the engine's own acceptance column is 0/0 on the master rank)."""
    isig = run.layout["sdparSWD"][0]
    moves = []
    for c in np.unique(post[:, -1]):
        sub = post[post[:, -1] == c]
        if len(sub) > 1:
            moves.append(np.mean(np.any(sub[1:, 3:isig] != sub[:-1, 3:isig], axis=1)))
    return float(np.mean(moves)) if moves else float("nan")


def sigma_columns(post, run):
    a, b = run.layout["sdparSWD"]
    return post[:, a:b]


def map_row(dat, k=None):
    """The sample row of maximum logL (optionally among rows with a given k)."""
    rows = dat if k is None else dat[dat[:, 3].astype(int) == k]
    if len(rows) == 0:
        raise ValueError(f"no sample with k = {k}")
    return rows[np.argmax(rows[:, 0])]
