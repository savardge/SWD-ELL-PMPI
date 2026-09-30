"""Node model -> layer stack, exactly as the engine builds it.

GETREF (loglhood.f90): piecewise-linear interpolation of the vel_ref rows
(Vs and Vp/Vs), constant beyond the last row; INTERPLAYER: a node with an
inactive (-100) parameter merges with the layer above, i.e. the value is
carried down; MAKE_CURMOD: Vp = (Vp/Vs) * Vs when I_VPVS = 1, Vp = 1.75 Vs
when I_VPVS = -1, density 2.35 + 0.036 (Vp - 3)^2 g/cc, Brocher (2005)
replacements with VP_BROCHER / RHO_BROCHER; LOGLHOOD_SWD: the last sampled
layer runs to the first tail depth, the tail rows take vel_prem Vs and Vp/Vs
shifted by the deepest node's dVs / dVpVs, their own density, half-space last.

Stacks are (n, 4) arrays: thickness km, density g/cc, Vp km/s, Vs km/s.
"""
import numpy as np


def getref(z, vel_ref):
    """(Vs_ref, VpVs_ref) at depth z km, the engine's GETREF."""
    zr, vs, r = vel_ref[:, 0], vel_ref[:, 1], vel_ref[:, 2]
    if z >= zr[-1]:
        return vs[-1], r[-1]
    if z == 0.0:
        return vs[0], r[0]
    i = int(np.searchsorted(zr, z, side="left")) - 1     # first row with zr >= z, minus one
    i = max(i, 0)
    f = (z - zr[i]) / (zr[i + 1] - zr[i])
    return vs[i] + f * (vs[i + 1] - vs[i]), r[i] + f * (r[i + 1] - r[i])


def legacy_rho(vp_kms):
    return 2.35 + 0.036 * (np.asarray(vp_kms, float) - 3.0) ** 2


def brocher_vp(vs_kms):
    """Brocher (2005, BSSA 95, eq. 9), km/s, valid 0 < Vs < 4.5."""
    v = np.asarray(vs_kms, float)
    return 0.9409 + 2.0947 * v - 0.8206 * v ** 2 + 0.2683 * v ** 3 - 0.0251 * v ** 4


def brocher_rho(vp_kms):
    """Brocher (2005, eq. 1), g/cc, Vp clamped to 1.5-8.5 km/s as APPLY_BROCHER does."""
    v = np.clip(np.asarray(vp_kms, float), 1.5, 8.5)
    return 1.6612 * v - 0.4721 * v ** 2 + 0.0671 * v ** 3 - 0.0043 * v ** 4 + 0.000106 * v ** 5


def sample_nodes(row, cfg):
    """Nodes of one sample row (or a map_voro vector prefixed by k at index 3),
    sorted by depth: depth, dVs, dVpVs (zeros when NPL = 2). An INACTIVE
    parameter (-100 in the file) is returned as NaN: the engine gives that
    node the ABSOLUTE value of the node above (INTERPLAYER), which
    absolute_nodes() applies once the reference is known."""
    nlmx, npl = cfg["NLMX"], cfg["NPL"]
    k = int(round(row[3]))
    v = np.array(row[4:4 + nlmx * npl], float).reshape(nlmx, npl)[:k]
    v = v[v[:, 0] > -99.0]
    v = v[np.argsort(v[:, 0], kind="stable")]
    v[:, 1:][v[:, 1:] < -99.0] = np.nan
    dvpvs = v[:, 2] if npl == 3 else np.zeros(len(v))
    return v[:, 0], v[:, 1], dvpvs


def map_nodes(mv, cfg):
    """Nodes of a map_voro dict (io.read_map_voro)."""
    row = np.concatenate([[0, 0, 0, mv["k"]], mv["voro"].ravel()])
    return sample_nodes(row, cfg)


def absolute_nodes(z, dvs, dvpvs, cfg, vel_ref):
    """Absolute Vs and Vp/Vs at the (depth-sorted) nodes: reference at the
    node depth + perturbation for an active parameter, the value of the node
    above for an inactive one (NaN), as INTERPLAYER does."""
    vs = np.empty(len(z))
    vpvs = np.empty(len(z))
    for i, zi in enumerate(z):
        vr, rr = getref(zi, vel_ref)
        if np.isfinite(dvs[i]):
            vs[i] = vr + dvs[i]
        else:
            vs[i] = vs[i - 1] if i > 0 else vr
        if cfg["I_VPVS"] != 1:
            vpvs[i] = 1.75
        elif np.isfinite(dvpvs[i]):
            vpvs[i] = rr + dvpvs[i]
        else:
            vpvs[i] = vpvs[i - 1] if i > 0 else rr
    return vs, vpvs


def vs_on_grid(z, dvs, zgrid, cfg, vel_ref):
    """Absolute Vs (km/s) on a depth grid; layer i spans [z_i, z_{i+1})."""
    order = np.argsort(z, kind="stable")
    z, dvs = np.asarray(z)[order], np.asarray(dvs)[order]
    vs, _ = absolute_nodes(z, dvs, np.zeros(len(z)), cfg, vel_ref)
    lay = np.clip(np.searchsorted(z, zgrid, "right") - 1, 0, len(z) - 1)
    return vs[lay]


def vpvs_on_grid(z, dvpvs, zgrid, cfg, vel_ref):
    order = np.argsort(z, kind="stable")
    z, dvpvs = np.asarray(z)[order], np.asarray(dvpvs)[order]
    _, r = absolute_nodes(z, np.zeros(len(z)), dvpvs, cfg, vel_ref)
    lay = np.clip(np.searchsorted(z, zgrid, "right") - 1, 0, len(z) - 1)
    return r[lay]


def engine_stack(z, dvs, dvpvs, cfg, vel_ref, vel_prem):
    """The layer stack LOGLHOOD_SWD forwards (curmod2), km / g/cc / km/s."""
    order = np.argsort(z, kind="stable")
    z, dvs, dvpvs = np.asarray(z, float)[order], np.asarray(dvs, float)[order], np.asarray(dvpvs, float)[order]
    vs, vpvs = absolute_nodes(z, dvs, dvpvs, cfg, vel_ref)
    vp = vpvs * vs
    rho = legacy_rho(vp)
    if cfg["VP_BROCHER"] == 1:
        vp = brocher_vp(vs)
    if cfg["RHO_BROCHER"] == 1:
        rho = brocher_rho(vp)
    th = np.r_[np.diff(z), vel_prem[0, 0] - z[-1]]
    # tail: vel_prem Vs / VpVs shifted by the deepest ACTIVE node's perturbations
    factvs = dvs[np.isfinite(dvs)][-1] if np.any(np.isfinite(dvs)) else 0.0
    factvpvs = (dvpvs[np.isfinite(dvpvs)][-1] if np.any(np.isfinite(dvpvs)) else 0.0) if cfg["I_VPVS"] == 1 else 0.0
    tvs = vel_prem[:, 1] + factvs
    tvp = tvs * (vel_prem[:, 2] + factvpvs)
    trho = vel_prem[:, 3].copy()
    if cfg["VP_BROCHER"] == 1:
        tvp = brocher_vp(tvs)
    if cfg["RHO_BROCHER"] == 1:
        trho = brocher_rho(tvp)
    tth = np.r_[np.diff(vel_prem[:, 0]), 0.0]
    return np.column_stack([np.r_[th, tth], np.r_[rho, trho], np.r_[vp, tvp], np.r_[vs, tvs]])
