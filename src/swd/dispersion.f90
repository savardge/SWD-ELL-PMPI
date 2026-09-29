subroutine dispersion_cu(nlyrs,rho,alpha,beta,thick,velc,velu,peri,NTMAX,ier,nmode_in,ivalid,cmin_in,cmax_in,dc_in,dc_over_in,iwarm_in)
!! Returns BOTH the phase velocity velc and the group velocity velu of mode
!! nmode_in at every period: DISPER80 evaluates the analytic group velocity
!! (energy integrals, RAYMRX with IG=3) at every converged phase root, so
!! both come out of one root search. dispersion() below is the original
!! single-output entry point (selects one of the two by IGRP).
! C ----------------------------------------------------------------------
! C Hrvoje Tkalcic, February 25, 2005, LLNL
! C The only input file is model.0              
! C Format - more or less self-explanatory, first line being number of layers
! C Uses the same logics as programs by Julia and Ammon but modified so that
! C it is useful for a grid search. This program can be executed from a 
! C shell script, which will form models in a grid search manner before each run.
! C The output will be 2 files, for Love and Rayleigh wave dispersion with
! C first column representing period, second - group and third - phase velocity.
! C Subroutines: lovdsp     -> DISPER80 (Saito)
! C              lovmrx     -> DISPER80 (Saito)
! C              raydsp     -> DISPER80 (Saito)
! C              raymrx     -> DISPER80 (Saito)
! C Libraries:   sac.a      -> SAC10.6
! C Sample run:  program 'dispersion' is invoked as follows
! C              dispersion > disp.out
! C ----------------------------------------------------------------------
      include 'params.h'
      PARAMETER (NLMAX=50,NCUMAX=2)
!C
      INTEGER   nlyrs, ilay, idum1, itr, iper, NTMAX,ier,iwarm,ilast
      REAL      alpha(nlyrs), beta(nlyrs),rho(nlyrs),&
               thick(nlyrs), smooth(nlyrs), beta_ap(nlyrs),&
               weight(nlyrs), cmin, cmax,& 
               dc, tol, pi2, w, ax(nlyrs), c(NTMAX), u(NTMAX), &
               velc(NTMAX), velu(NTMAX), ekd, y0l(6),&
               vp(2*nlyrs), vs(2*nlyrs), z(2*nlyrs), za,&
               cwarm, pwarm, cwback, cwfwd, clow,&
               y0r(3), yij(15), ap(nlyrs), ae(nlyrs), peri (NTMAX)
      CHARACTER*20 title
!C
      EXTERNAL  lovmrx, raymrx
!! Mode index: 0 = fundamental (= original DISPER80 behaviour), 1 = first
!! higher mode, ... Resolved by RAYDSPN's sign-change counting. REQUIRED, not
!! optional: dispersion() is a bare external subroutine with no explicit
!! interface at its call site, and OPTIONAL/PRESENT is invalid without one.
      INTEGER, INTENT(IN) :: nmode_in
!! Per-period validity: 1 = root found, 0 = not found (mode below cut-off, or
!! leaky i.e. c > half-space Vs, which DISPER80 legitimately refuses).
!! NEEDED because ier below is overwritten on every period, so on its own it
!! only ever reports the status of the LAST period - a latent bug that becomes
!! important for overtones, which routinely fail near their cut-off.
      INTEGER, INTENT(OUT) :: ivalid(NTMAX)
!! Root-scan window and step [km/s] (keyword SWD_SCAN in the parameter file;
!! defaults 2.0 6.5 0.05 = the original crustal values; near-surface work
!! needs e.g. 0.08 1.6 0.005).
      REAL, INTENT(IN) :: cmin_in, cmax_in, dc_in, dc_over_in
!! Warm-start switch: 1 = march the root from the previous period (fast),
!! 0 = original cold scan from cmin at every period.  See the block below.
      INTEGER, INTENT(IN) :: iwarm_in
!C
      DATA      pi2/6.283185/, ia/0/
      cmin=cmin_in
      cmax=cmax_in
      dc=dc_in
      tol=0.0001
      itr=10
! CCCCCCCCCCCCCCCCCCCCCCCC
! C Reads in Earth model C
! CCCCCCCCCCCCCCCCCCCCCCCC

    za=0.
  !    OPEN (1,file='model.0',status='Old')
  !    READ (1,*) nlyrs, title

!  DO ilay=1,nlyrs
!  !        READ (1,*) idum1,alpha(ilay),beta(ilay),rho(ilay),thick(ilay),
!   !  &         smooth(ilay),beta_ap(ilay),weight(ilay)
!        
!        alpha(ilay)=vpvs*beta(ilay)
!        !rho(ilay)=3
!        rho(ilay)= (2.35+0.036*(beta(ilay)*vpvs-3.0)**2)
!  END DO

DO ilay=1,nlyrs
  z(2*ilay-1)=za
  vs(2*ilay-1)=beta(ilay)
  vp(2*ilay-1)=alpha(ilay)
  z(2*ilay)=z(2*ilay-1)-thick(ilay)
  vs(2*ilay)=beta(ilay)
  vp(2*ilay)=alpha(ilay)
  za=z(2*ilay)
ENDDO

!!write (*,*)alpha,beta,rho,thick

!! Overtones need a finer scan: R1/R2 approach to ~12 m/s near 15.5 Hz on the
!! HVC dam model, i.e. 2 steps at dc=0.005. At dc=0.001 mode-2 agreement with
!! disba improves from 1.59 to 0.01 m/s. The fundamental is unaffected
!! (identical to 0.00 m/s at either step), so keep it cheap.
IF (nmode_in > 0) dc = dc_over_in

ivalid = 0

!! ---------------------------------------------------------------------
!! Warm-started root search (2026-09-02).
!! DISPER80 finds the MODE-th root by scanning phase velocity from cmin in
!! steps of dc and counting sign changes of the secular function; at the
!! overtone step (dc = 0.001 km/s) that is ~850 propagator evaluations per
!! period, and it is repeated from cmin for every period and every mode.
!! Rayleigh phase velocity of a given mode increases with period, so when
!! the periods are supplied in ASCENDING order the previous period's root is
!! a lower bound for the current one: the scan then starts just below it
!! instead of at cmin.  Measured on the HVC posteriors: 8-fold fewer
!! propagator calls for the overtones, 3-fold for the fundamental.
!!
!! Exactness: the scan GRID is unchanged (the start is the grid point
!! cmin + K*dc just below the bound), so the bracket -- and therefore the
!! refined root -- are bit-identical to the cold scan whenever no root lies
!! between cmin and the start. Two guards make that self-correcting:
!!   (a) if no crossing is found above the warm start, and
!!   (b) if the root found lies further than CWFWD above the previous one
!!       (a jump that large means the warm bracket cannot be trusted),
!! the period is recomputed with a full cold scan, which is authoritative.
!! RAYDSPN adds a sign-parity certificate on the skipped interval, and the
!! whole curve is recomputed cold if it is not monotone in period.
!! LIMIT: strongly inverse-dispersive models (a fast lid over a slow
!! channel) can move a root DOWN by more than the warm window and skip a
!! PAIR of roots, which parity cannot see.  Such models only occur when the
!! adjacent-layer contrast is unconstrained, so the caller switches the warm
!! start off there (iwarm_in = 0); see LOGLHOOD_SWD.
!! ---------------------------------------------------------------------
iwarm  = iwarm_in    !! 1 = warm start allowed on this pass, 0 = cold pass
 900 CONTINUE            !! re-entry point for the cold verification pass
cwarm  = -1.
pwarm  = -1.
cwback = 10.*dc      !! safety margin below the previous root [km/s]
cwfwd  = 0.30        !! max trusted rise in c between adjacent periods [km/s]

DO iper=1,NTMAX
  w=pi2/peri(iper)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C Call DISPER80 for Love dispersion C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALL lovdsp (lovmrx,thick,rho,beta,ax,nlyrs,w,
! &                         cmin,cmax,dc,
! &                         tol,itr,ia,c(iper),
! &                         u(iper),ekd,y0l,ier)
! write (8,77)t(iper),u(iper),c(iper)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C Call DISPER80 for Rayleigh dispersion C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!write(*,*)'------------------------'
  c(iper)=100.
  u(iper)=100.

  !write (*,*) 'thick',thick
  !write (*,*) 'rho',rho
  !write (*,*) 'alpha',alpha
  !write (*,*) 'beta',beta
  !write (*,*) 'misc1',nlyrs,cmin,cmax,dc,tol
  !write (*,*) 'itr',itr
  !write (*,*) 'ia',ia
  !write (*,*) 'ekd',ekd
  !write (*,*) 'y0r',y0r
  !write (*,*) 'yij',yij

  c(iper) = 0.
  u(iper) = 0.
  !! Warm start only when marching to a LONGER period (c increases); any
  !! other ordering falls back to the cold scan.
  clow = 0.
  IF (iwarm == 1 .AND. cwarm > 0. .AND. peri(iper) > pwarm) clow = cwarm - cwback
  CALL raydspn(raymrx,thick,rho,alpha,beta,ap,ae,nlyrs,&
               w,cmin,cmax,dc,tol,itr,ia,nmode_in,clow,c(iper),u(iper),&
               ekd,y0r,yij,ier)
  IF (clow > 0.) THEN
    !! Guards (a) and (b): recompute cold when the warm bracket is not
    !! trustworthy. The cold result then replaces it unconditionally.
    IF (ier /= 0 .OR. c(iper) > cwarm + cwfwd) THEN
      c(iper) = 0.
      u(iper) = 0.
      CALL raydspn(raymrx,thick,rho,alpha,beta,ap,ae,nlyrs,&
                   w,cmin,cmax,dc,tol,itr,ia,nmode_in,0.,c(iper),u(iper),&
                   ekd,y0r,yij,ier)
    ENDIF
  ENDIF
  IF(ier < 0)THEN
    !! Hard input error: abort the whole call (ier propagates to the caller).
    ivalid(iper:NTMAX) = 0
    velc = 0.
    velu = 0.
    RETURN
  ENDIF
  !! ier = 1 (slow convergence) and 2 (root not found) both mean "no usable
  !! datum here"; this matches the original all-or-nothing treatment of ier/=0,
  !! but now resolved PER PERIOD instead of for the whole curve.
  IF(ier == 0)THEN
    ivalid(iper) = 1
    cwarm = c(iper)
    pwarm = peri(iper)
  ELSE
    ivalid(iper) = 0
  ENDIF
  !write (*,*) w,u(iper),c(iper)
  !write(*,*)'------------------------'
ENDDO
!write (*,*) 'c',c
!write (*,*) 'u',u
!! JD
!! ---------------------------------------------------------------------
!! Monotonicity audit.  The warm start is exact while the mode's phase
!! velocity increases with period, which is the normal case (and was true
!! for every model of the constrained HVC posteriors: 7000 models, minimum
!! step +0.5 m/s).  Strongly inverse-dispersive models -- e.g. a fast lid
!! over a slow channel, which only exist when the adjacent-layer contrast
!! is unconstrained -- can make c DECREASE with period by more than the
!! warm window, and then the warm bracket may skip a pair of roots without
!! tripping the parity certificate inside RAYDSPN.  Any decrease in the
!! computed curve is that signature, so the whole curve is recomputed with
!! cold scans, which are authoritative.  Only pathological models pay.
!! ---------------------------------------------------------------------
IF (iwarm == 1) THEN
  ilast = 0
  DO iper = 1,NTMAX
    IF (ivalid(iper) == 1) THEN
      IF (ilast > 0) THEN
        IF (c(iper) < c(ilast) .AND. peri(iper) > peri(ilast)) THEN
          iwarm = 0
          ivalid = 0
          GO TO 900
        ENDIF
      ENDIF
      ilast = iper
    ENDIF
  ENDDO
ENDIF

ier = 0
velc = c
velu = u
RETURN
END SUBROUTINE dispersion_cu

!! ---------------------------------------------------------------------
!! Original entry point: one output array, phase (IGRP = 0) or group
!! (IGRP = 1) velocity. Thin wrapper around dispersion_cu, so every caller
!! that only needs one of the two is unchanged and bit-identical.
!! ---------------------------------------------------------------------
subroutine dispersion(nlyrs,rho,alpha,beta,thick,vel,peri,NTMAX,IGRP,ier,nmode_in,ivalid,cmin_in,cmax_in,dc_in,dc_over_in,iwarm_in)
      INTEGER   nlyrs, NTMAX, IGRP, ier
      INTEGER, INTENT(IN)  :: nmode_in, iwarm_in
      INTEGER, INTENT(OUT) :: ivalid(NTMAX)
      REAL      alpha(nlyrs), beta(nlyrs), rho(nlyrs), thick(nlyrs)
      REAL      vel(NTMAX), peri(NTMAX)
      REAL, INTENT(IN) :: cmin_in, cmax_in, dc_in, dc_over_in
      REAL      velc(NTMAX), velu(NTMAX)
CALL dispersion_cu(nlyrs,rho,alpha,beta,thick,velc,velu,peri,NTMAX,ier,nmode_in,ivalid,&
                   cmin_in,cmax_in,dc_in,dc_over_in,iwarm_in)
IF(IGRP == 1)THEN
  vel = velu
ELSE
  vel = velc
ENDIF
RETURN
END SUBROUTINE dispersion
