"""Shallow (0-100 m) test models for the SWD forward-model validation.

Each model is a layer stack (thickness km, rho g/cm3, vp km/s, vs km/s) over a
half-space, chosen to exercise the cases a tailings-dam inversion actually
visits: normal dispersion, smooth and steep gradients, low-velocity zones of
increasing severity, a fast lid over a slow channel (inverse dispersion), and
saturated profiles where Vp/Vs is large.

The last entry of every stack is the half-space; its thickness is ignored by
both codes.
"""

import numpy as np


def _rho_from_vs(vs):
    """Plausible near-surface bulk density [g/cm3]; monotone in Vs, 1.6-2.2."""
    return np.clip(1.6 + 0.9 * vs, 1.6, 2.2)


def _stack(layers, vp_mode="dry"):
    """layers = [(thickness_m, vs_ms), ...]; last entry is the half-space.

    vp_mode 'dry'       : Vp = 1.9 * Vs        (unsaturated, Poisson ~ 0.31)
            'saturated' : Vp = max(1.5 km/s, 1.9 Vs)  (below the water table,
                          Poisson -> 0.49; a stress test for the Rayleigh
                          secular function)
    """
    h = np.array([l[0] for l in layers], dtype=float) / 1000.0  # m -> km
    vs = np.array([l[1] for l in layers], dtype=float) / 1000.0  # m/s -> km/s
    if vp_mode == "dry":
        vp = 1.9 * vs
    elif vp_mode == "saturated":
        vp = np.maximum(1.5, 1.9 * vs)
    else:
        raise ValueError(vp_mode)
    rho = _rho_from_vs(vs)
    return np.column_stack([h, rho, vp, vs])


def _gradient(z_top, z_bot, vs_top, vs_bot, nlay, power=1.0):
    """Discretise a Vs gradient into nlay equal-thickness layers."""
    dz = (z_bot - z_top) / nlay
    out = []
    for i in range(nlay):
        f = (i + 0.5) / nlay
        vs = vs_top + (vs_bot - vs_top) * f**power
        out.append((dz, vs))
    return out


def build_models():
    """Return an ordered dict {name: (nlayer, 4) array}."""
    m = {}

    # --- normal dispersion, few layers ---------------------------------
    m["2layer_normal"] = _stack([(50, 200), (0, 600)])
    m["3layer_normal"] = _stack([(20, 150), (40, 300), (0, 700)])

    # --- gradients ------------------------------------------------------
    m["grad_smooth"] = _stack(_gradient(0, 100, 120, 500, 20) + [(0, 600)])
    m["grad_steep"] = _stack(_gradient(0, 100, 100, 700, 20) + [(0, 800)])
    # sqrt-like compaction gradient, fine discretisation (2 m layers)
    m["grad_sqrt_fine"] = _stack(
        _gradient(0, 100, 110, 480, 50, power=0.5) + [(0, 650)]
    )

    # --- low-velocity zones ---------------------------------------------
    m["lvz_weak"] = _stack([(20, 250), (20, 200), (60, 400), (0, 700)])
    m["lvz_strong"] = _stack([(15, 350), (25, 180), (60, 500), (0, 800)])
    m["lvz_deep"] = _stack([(40, 250), (30, 170), (30, 450), (0, 750)])
    # fast lid over a slow channel: inverse dispersion, the case the warm
    # start explicitly does not trust (LOGLHOOD_SWD turns it off)
    m["lid_over_channel"] = _stack([(10, 400), (40, 180), (0, 600)])
    # LVZ inside a gradient
    m["grad_with_lvz"] = _stack(
        _gradient(0, 40, 130, 300, 8)
        + [(20, 200)]
        + _gradient(60, 100, 380, 520, 8)
        + [(0, 700)]
    )

    # --- realistic tailings profile --------------------------------------
    m["tailings_dry"] = _stack(
        [(5, 120), (15, 180), (30, 250), (50, 350), (0, 600)]
    )
    m["tailings_saturated"] = _stack(
        [(5, 120), (15, 180), (30, 250), (50, 350), (0, 600)],
        vp_mode="saturated",
    )
    m["lvz_saturated"] = _stack(
        [(15, 350), (25, 180), (60, 500), (0, 800)], vp_mode="saturated"
    )

    # --- near-homogeneous (higher modes should mostly not exist) ---------
    m["near_homogeneous"] = _stack([(100, 300), (0, 320)])

    # The two models above have only 2 entries in the stack, which DISPER80
    # refuses (RAYDSPN: IF L.LE.2 -> IER = -1).  These are the identical
    # physical models with the top layer split in two, to confirm that the
    # rejection is the layer count and nothing else.
    m["2layer_normal_split"] = _stack([(25, 200), (25, 200), (0, 600)])
    m["near_homogeneous_split"] = _stack([(50, 300), (50, 300), (0, 320)])

    return m


def write_model_file(path, model):
    """One model per file, in the driver's format."""
    with open(path, "w") as f:
        f.write(f"{len(model)}\n")
        for h, rho, vp, vs in model:
            f.write(f"{h:.9f} {rho:.9f} {vp:.9f} {vs:.9f}\n")


def scan_window(model, pad_lo=0.80, pad_hi=0.999):
    """Root-scan window [km/s] for this model.

    Rayleigh phase velocity is bounded below by ~0.86*min(Vs) and above by the
    half-space Vs (roots above it are leaky and DISPER80 legitimately refuses
    them), so the window is set from the model itself.
    """
    vs = model[:, 3]
    return pad_lo * vs.min(), pad_hi * vs[-1]
