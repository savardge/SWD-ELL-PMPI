!! Validation + benchmark for the warm-started DISPER80 root scan.
!!
!! Reads a model list (tools/export_test_models.py: layer stack + PREM
!! continuation exactly as LOGLHOOD_SWD builds it) and one mode's observed
!! periods, then computes the modal phase velocities twice:
!!   warm : dispersion()      -- period-marching warm start (production path)
!!   cold : raydspn(clow = 0) -- the original scan from cmin at every period
!! and reports every difference, the validity pattern, the period-to-period
!! step distribution of the reference curve, and the wall time of each path.
!!
!! Build: make -f Makefile.testwarm   Run: ./test_warm models.txt periods.txt MODE
PROGRAM test_warm
  IMPLICIT NONE
  INTEGER, PARAMETER :: NTMAX = 200, NLMAX = 200, NMODMAX = 40000
  REAL    :: peri(NTMAX), cw(NTMAX), cc(NTMAX)
  INTEGER :: ivalid(NTMAX), ivc(NTMAX)
  REAL    :: thick(NLMAX), rho(NLMAX), vp(NLMAX), vs(NLMAX), ap(NLMAX), ae(NLMAX)
  REAL    :: y0r(3), yij(15), ekd, u1, c1, dc, dcm, cmin, cmax, tol
  REAL    :: dmax, w, pi2, dstep, smin, smax, thr(4)
  INTEGER :: nt, nl, imod, nmod, iper, ier, mode, itr, ia, nbad, nmiss, i
  INTEGER :: nstep, nb(4), ib
  REAL    :: t0, t1, twarm, tcold
  CHARACTER(LEN=256) :: fmod, fper, arg
  EXTERNAL raymrx

  CALL GET_COMMAND_ARGUMENT(1, fmod)
  CALL GET_COMMAND_ARGUMENT(2, fper)
  CALL GET_COMMAND_ARGUMENT(3, arg); READ(arg,*) mode

  cmin = 0.08; cmax = 1.6; tol = 0.0001; itr = 10; ia = 0
  dc = 0.005; dcm = 0.001
  IF (mode > 0) dc = dcm
  pi2 = 6.283185; ap = 0.; ae = 0.
  thr = (/ 0.02, 0.05, 0.10, 0.20 /)

  OPEN(10, FILE=fper, STATUS='OLD', ACTION='READ')
  nt = 0
  DO
    READ(10,*,IOSTAT=ier) peri(nt+1)
    IF (ier /= 0) EXIT
    nt = nt + 1
  ENDDO
  CLOSE(10)

  OPEN(11, FILE=fmod, STATUS='OLD', ACTION='READ')
  nmod = 0; nbad = 0; nmiss = 0; dmax = 0.
  nstep = 0; nb = 0; smin = 1.e30; smax = -1.e30
  twarm = 0.; tcold = 0.
  DO imod = 1, NMODMAX
    READ(11,*,IOSTAT=ier) nl
    IF (ier /= 0) EXIT
    DO i = 1, nl
      READ(11,*) thick(i), rho(i), vp(i), vs(i)
    ENDDO
    nmod = nmod + 1

    CALL CPU_TIME(t0)
    CALL dispersion(nl, rho, vp, vs, thick, cw, peri, nt, 0, ier, mode, ivalid, &
                    cmin, cmax, 0.005, dcm, 1)
    CALL CPU_TIME(t1); twarm = twarm + (t1 - t0)

    CALL CPU_TIME(t0)
    DO iper = 1, nt
      w = pi2 / peri(iper)
      c1 = 0.; u1 = 0.
      CALL raydspn(raymrx, thick, rho, vp, vs, ap, ae, nl, w, cmin, cmax, dc, &
                   tol, itr, ia, mode, 0.0, c1, u1, ekd, y0r, yij, ier)
      cc(iper) = c1
      ivc(iper) = 0
      IF (ier == 0) ivc(iper) = 1
    ENDDO
    CALL CPU_TIME(t1); tcold = tcold + (t1 - t0)

    IF (SUM(ivalid(1:nt)) /= SUM(ivc(1:nt))) nmiss = nmiss + 1
    DO iper = 1, nt
      IF (ivalid(iper) == 1 .AND. ivc(iper) == 1 .AND. cw(iper) /= cc(iper)) THEN
        nbad = nbad + 1
        dmax = MAX(dmax, ABS(cw(iper) - cc(iper)))
      ENDIF
    ENDDO
    DO iper = 2, nt
      IF (ivc(iper) == 1 .AND. ivc(iper-1) == 1) THEN
        dstep = cc(iper) - cc(iper-1)
        smin = MIN(smin, dstep); smax = MAX(smax, dstep); nstep = nstep + 1
        DO ib = 1, 4
          IF (ABS(dstep) > thr(ib)) nb(ib) = nb(ib) + 1
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  CLOSE(11)

  WRITE(*,'(A,I6,A,I3,A,I2)') 'models ', nmod, '  periods ', nt, '  mode ', mode
  WRITE(*,'(A,I8,A,ES11.3)')  'periods warm /= cold : ', nbad, '   max |diff| km/s ', dmax
  WRITE(*,'(A,I6)')           'models with different validity pattern : ', nmiss
  WRITE(*,'(A,I8,A,F8.4,A,F8.4)') 'cold steps ', nstep, '  min ', smin, '  max ', smax
  WRITE(*,'(A,4I9)')          '|step| > .02 .05 .10 .20 km/s : ', nb
  WRITE(*,'(A,F8.3,A,F8.3,A,F7.2,A)') 'cpu warm ', twarm, ' s  cold ', tcold, &
        ' s  speedup ', tcold/MAX(twarm,1.e-9), 'x'
END PROGRAM test_warm
