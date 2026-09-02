#!/usr/bin/env python3
"""Convert run inputs from the legacy RF-SWD-ELL-MT-PMPI format (81-line
parameter file, 54-line covparameter file) or from the
receiver_rjmcmc_varpar_sourceinv_joint format (44/46-line parameter file) to
the SWD-ELL-PMPI format:

  <base>_parameter.dat     48 positional lines + optional keyword tail
                           (DVSCON, MODE_OF, IGRP, SWD_SCAN)
  <base>_covparameter.dat  35 lines (SWD and ELL datasets only)
  <base>_map_voro.dat      k, voro(NLMX*NPL), sdparSWD(NMODE), sdparELL(NMODE_ELL),
                           arparSWD(NMODE), arparELL(NMODE_ELL)

Usage:  convert_legacy_inputs.py {rfswdellmt|receiver} SRC_DIR DST_DIR
"""
import os, re, shutil, sys

fmt, src, dst = sys.argv[1:4]
os.makedirs(dst, exist_ok=True)
fb = open(os.path.join(src, 'filebase.txt')).read().split()[-1]
L = open(os.path.join(src, f'{fb}_parameter.dat')).read().splitlines()
v = lambda i: L[i-1].split('!!')[0].strip()          # 1-based line -> value text

if fmt == 'rfswdellmt':
    g = dict(IMAP=v(1), IMAGSCALE=v(3), ENOS=v(4), IPOIPR=v(5), IAR=v(6), I_VARPAR=v(7), IBD_SINGLE=v(8),
             I_SWD=v(11), I_ELL=v(12), I_VREF=v(15), I_VPVS=v(16), ISMPPRIOR=v(17), ISETSEED=v(18), IEXCHANGE=v(19),
             NDAT_SWD=v(21), NMODE=v(22), NDAT_ELL=v(23), NMODE_ELL=v(24), NLMN=v(28), NLMX=v(29), ICHAINTHIN=v(30),
             NKEEP=v(31), NPTCHAINS1=v(32), dTlog=v(33), lambda_=v(34), hmx=v(35), hmin=v(36), armxSWD=v(39), armxELL=v(40),
             TCHCKPT=v(41), dVs=v(46), dVpVs=v(47), sdmn=' '.join(v(50).split()[1:3]), sdmx=' '.join(v(51).split()[1:3]),
             ISD_SWD=v(61), ISD_ELL=v(62), ICOV_SWD=v(65), ICOV_ELL=v(66), ELL_verbose=v(68), ELL_prec=v(69),
             I_ABS_ELL=v(70), I_LOG10_ELL=v(71), I_SAMPLING_TYPE_ELL=v(72), I_SET_STEP_ELL=v(73), STEP_SIZE_ELL=v(74),
             I_SET_COUNT_ELL=v(75), COUNT_ELL=v(76), I_SET_RANGE_ELL=v(77))
    tail = [l for l in L[77:] if re.match(r'^\s*[A-Za-z]', l)]      # existing keyword lines, if any
    nmode, nmode_ell = int(g['NMODE']), int(g['NMODE_ELL'])
    swd_scan = None
elif fmt == 'receiver':
    sd1, sd2 = v(40).split(), v(41).split()
    g = dict(IMAP=v(1), IMAGSCALE='0', ENOS=v(3), IPOIPR=v(4), IAR=v(5), I_VARPAR=v(6), IBD_SINGLE=v(7),
             I_SWD=v(10), I_ELL='0', I_VREF=v(11), I_VPVS=v(12), ISMPPRIOR=v(13), ISETSEED=v(14), IEXCHANGE=v(15),
             NDAT_SWD=v(17), NMODE=v(18), NDAT_ELL='1', NMODE_ELL='1', NLMN=v(21), NLMX=v(22), ICHAINTHIN=v(23),
             NKEEP=v(24), NPTCHAINS1=v(25), dTlog=v(26), lambda_=v(27), hmx=v(28), hmin=v(29), armxSWD=v(32), armxELL='0.1',
             TCHCKPT=v(33), dVs=v(38), dVpVs=v(39), sdmn=f'{sd1[1]} 1.0e-3', sdmx=f'{sd2[1]} 1.0e-1',
             ISD_SWD='1', ISD_ELL='0', ICOV_SWD='1', ICOV_ELL='1', ELL_verbose='0', ELL_prec='1.e-5',
             I_ABS_ELL='1', I_LOG10_ELL='0', I_SAMPLING_TYPE_ELL='0', I_SET_STEP_ELL='0', STEP_SIZE_ELL='1.5',
             I_SET_COUNT_ELL='0', COUNT_ELL='145', I_SET_RANGE_ELL='1')
    tail = [f'IGRP {v(44)}']
    if len(L) >= 45 and re.match(r'^\s*[-0-9.]', L[44]): tail.append(f'DVSCON {v(45)}')
    if len(L) >= 46 and re.match(r'^\s*[0-9]', L[45]): tail.append(f'MODE_OF {v(46)}')
    tail.append('SWD_SCAN 0.08 1.6 0.005 0.001')
    nmode, nmode_ell = int(g['NMODE']), 1
else:
    sys.exit('fmt must be rfswdellmt or receiver')

order = ['IMAP','IMAGSCALE','ENOS','IPOIPR','IAR','I_VARPAR','IBD_SINGLE','I_SWD','I_ELL','I_VREF','I_VPVS','ISMPPRIOR',
         'ISETSEED','IEXCHANGE','NDAT_SWD','NMODE','NDAT_ELL','NMODE_ELL','NLMN','NLMX','ICHAINTHIN','NKEEP','NPTCHAINS1',
         'dTlog','lambda_','hmx','hmin','armxSWD','armxELL','TCHCKPT','dVs','dVpVs','sdmn','sdmx','ISD_SWD','ISD_ELL',
         'ICOV_SWD','ICOV_ELL','ELL_verbose','ELL_prec','I_ABS_ELL','I_LOG10_ELL','I_SAMPLING_TYPE_ELL','I_SET_STEP_ELL',
         'STEP_SIZE_ELL','I_SET_COUNT_ELL','COUNT_ELL','I_SET_RANGE_ELL']
assert len(order) == 48
out = [f"{g[k]:<18} !! {i+1:2d} {k.rstrip('_')}" for i, k in enumerate(order)]
out.append('!! ---- optional keyword lines: DVSCON x | MODE_OF m1 m2 .. | IGRP 0|1 | SWD_SCAN cmin cmax dc [dc_over] ----')
out += tail
open(os.path.join(dst, f'{fb}_parameter.dat'), 'w').write('\n'.join(out) + '\n')

# ---- covparameter (SWD, ELL) ----
if fmt == 'rfswdellmt' and os.path.exists(os.path.join(src, f'{fb}_covparameter.dat')):
    C = [l.split('!!')[0].strip() for l in open(os.path.join(src, f'{fb}_covparameter.dat')).read().splitlines()]
    c = lambda i: C[i-1]
    two = lambda s: ' '.join(s.split()[1:3])
    cov = [c(2), c(3), c(5), c(6), c(7), c(8), c(9), c(10), c(12), c(13), two(c(15)), two(c(16)), two(c(17)),
           c(18), c(19), c(20), c(21), c(22), c(23), c(25), c(26), c(29), c(30)] + C[37:43] + C[43:49]
else:
    cov = ['0','0','2000','2000','100','0','5','4','0','0','1. 1.','10. 10.','2. 2.','100','100','0','0','0','1','2','2',
           '1.e-4','1.e-4','10','40','1','0','0','1.2','10','40','1','0','0','1.2']
names = ['Icov_iterUpdate_SWD','Icov_iterUpdate_ELL','covIter_zero_nsamples','covIter_period','MAXcovIter','ICOVest',
         'CHAINTHIN_COVest_period_zeroIter','CHAINTHIN_COVest_period_nonzeroIter','ISD_SWD_covIter','ISD_ELL_covIter',
         'sdmn_covIter (SWD ELL)','sdmx_covIter (SWD ELL)','sdpar_covIter (SWD ELL)','NKEEP_covIter','NKEEP_covIter_res',
         'iSAVEsample_covIter','iSAVEsample_only_zeroIter','iMAP_calc','iconverge_criterion','iconverge_criterion_SWD',
         'iconverge_criterion_ELL','converge_threshold_SWD','converge_threshold_ELL','nfrac_SWD','MAX_NAVE_SWD',
         'inonstat_SWD','iunbiased_SWD','imr_SWD','damp_power_SWD','nfrac_ELL','MAX_NAVE_ELL','inonstat_ELL',
         'iunbiased_ELL','imr_ELL','damp_power_ELL']
assert len(cov) == 35 == len(names), (len(cov), len(names))
open(os.path.join(dst, f'{fb}_covparameter.dat'), 'w').write(''.join(f'{a:<14} !!{b}\n' for a, b in zip(cov, names)))

# ---- map_voro ----
nlmx = int(g['NLMX']); npl = 3 if g['I_VPVS'] == '1' else 2; nv = nlmx * npl
m = open(os.path.join(src, f'{fb}_map_voro.dat')).read().split()
k, voro, rest = m[0], m[1:1+nv], m[1+nv:]
if fmt == 'receiver':   # rest: sdR,sdV,sdT, sdSWD(nmode), ar(3), arSWD(nmode)
    sdswd, arswd = rest[3:3+nmode], rest[6+nmode:6+2*nmode]; sdell, arell = ['0.05'], ['0.0']
else:                   # rest: sdR,sdV,sdT(NRF1=1 each), sdSWD, sdELL, sdMT, ar(3), arSWD, arELL
    sdswd = rest[3:3+nmode]; sdell = rest[3+nmode:3+nmode+nmode_ell]
    o = 3+nmode+nmode_ell+1+3; arswd = rest[o:o+nmode]; arell = rest[o+nmode:o+nmode+nmode_ell]
    if len(arell) < nmode_ell: arell = ['0.0'] * nmode_ell
    if len(arswd) < nmode: arswd = ['0.0'] * nmode
newmap = [k] + voro + sdswd + sdell + arswd + arell
assert len(newmap) == 1 + nv + 2*nmode + 2*nmode_ell, len(newmap)
open(os.path.join(dst, f'{fb}_map_voro.dat'), 'w').write(' '.join(newmap) + '\n')
for f in os.listdir(src):
    if f.startswith(f'{fb}_SWD') or f.startswith(f'{fb}_ELL') or f == f'{fb}_vel_ref.txt' or f == 'filebase.txt':
        shutil.copy(os.path.join(src, f), os.path.join(dst, f))
print(f'converted [{fmt}] {src} -> {dst}: NMODE={nmode}, map {len(newmap)} numbers, tail {tail}')
