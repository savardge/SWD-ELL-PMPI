#!/usr/bin/env python3
"""Export layer stacks of real posterior models for swd/test_warm_driver.f90.

Rebuilds each sampled model exactly as MAKE_CURMOD + LOGLHOOD_SWD do
(k Voronoi nodes -> k layers, the last extended to the first PREM depth,
then the PREM continuation perturbed by the deepest node's dVs/dVpVs) and
writes them as  nl / (thick[km] rho vp[km/s] vs[km/s]) x nl  records.

Usage: export_test_models.py SAMPLE_FILE VEL_REF NLMX NMODE N_MODELS OUT [--fmt receiver|swdell]
"""
import sys
import numpy as np

samp, vref_file, nlmx, nmode, nmodels, out = sys.argv[1:7]
nlmx, nmode, nmodels = int(nlmx), int(nmode), int(nmodels)
fmt = sys.argv[8] if len(sys.argv) > 8 else "receiver"

with open(vref_file) as f:
    nref = int(f.readline().split()[0])
    rows = np.array([[float(x) for x in f.readline().split()] for _ in range(nref + 2)])
vel_ref, vel_prem = rows[:nref], rows[nref:]
zref, vsref = vel_ref[:, 0], vel_ref[:, 1]

dat = np.loadtxt(samp)
rng = np.random.default_rng(0)
sel = dat[rng.choice(len(dat), min(nmodels, len(dat)), replace=False)]

with open(out, "w") as fh:
    for row in sel:
        k = int(row[3])
        v = row[4:4 + nlmx * 3].reshape(nlmx, 3)[:k]
        v = v[v[:, 0] > -99.0]
        v = v[np.argsort(v[:, 0])]
        z, dvs, dvpvs = v[:, 0], v[:, 1], v[:, 2]
        vs = np.interp(z, zref, vsref) + dvs          # GETREF + perturbation, km/s
        vpvs = 2.0 + dvpvs
        factvs, factvpvs = dvs[-1], dvpvs[-1]         # deepest active node
        th = np.empty(len(z))
        th[:-1] = np.diff(z)
        th[-1] = vel_prem[0, 0] - z[-1]               # last layer runs to the first PREM depth
        # PREM continuation (thicknesses between PREM depths, half-space last)
        th_p = np.append(np.diff(vel_prem[:, 0]), 0.0)
        vs_p = vel_prem[:, 1] + factvs
        vpvs_p = vel_prem[:, 2] + factvpvs
        TH = np.concatenate([th, th_p])
        VS = np.concatenate([vs, vs_p])
        VP = np.concatenate([vs * vpvs, vs_p * vpvs_p])
        RHO = 2.35 + 0.036 * (VP - 3.0) ** 2          # Birch-type law used in MAKE_CURMOD
        RHO[len(th):] = vel_prem[:, 3]                # PREM rows carry their own density
        fh.write(f"{len(TH)}\n")
        for t, r, p, s in zip(TH, RHO, VP, VS):
            fh.write(f"{t:.8f} {r:.8f} {p:.8f} {s:.8f}\n")
print(f"wrote {len(sel)} models ({out})")
