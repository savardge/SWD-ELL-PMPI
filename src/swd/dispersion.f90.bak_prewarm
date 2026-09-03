subroutine dispersion(nlyrs,rho,alpha,beta,thick,vel,peri,NTMAX,IGRP,ier,nmode_in,ivalid,cmin_in,cmax_in,dc_in,dc_over_in)
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
      INTEGER   nlyrs, ilay, idum1, itr, iper, NTMAX,ier
      REAL      alpha(nlyrs), beta(nlyrs),rho(nlyrs),&
               thick(nlyrs), smooth(nlyrs), beta_ap(nlyrs),&
               weight(nlyrs), cmin, cmax,& 
               dc, tol, pi2, w, ax(nlyrs), c(NTMAX), u(NTMAX), &
               vel(NTMAX), ekd, y0l(6),&
               vp(2*nlyrs), vs(2*nlyrs), z(2*nlyrs), za,&
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
  CALL raydspn(raymrx,thick,rho,alpha,beta,ap,ae,nlyrs,&
               w,cmin,cmax,dc,tol,itr,ia,nmode_in,c(iper),u(iper),&
               ekd,y0r,yij,ier)
  IF(ier < 0)THEN
    !! Hard input error: abort the whole call (ier propagates to the caller).
    ivalid(iper:NTMAX) = 0
    vel = 0.
    RETURN
  ENDIF
  !! ier = 1 (slow convergence) and 2 (root not found) both mean "no usable
  !! datum here"; this matches the original all-or-nothing treatment of ier/=0,
  !! but now resolved PER PERIOD instead of for the whole curve.
  IF(ier == 0)THEN
    ivalid(iper) = 1
  ELSE
    ivalid(iper) = 0
  ENDIF
  !write (*,*) w,u(iper),c(iper)
  !write(*,*)'------------------------'
ENDDO
!write (*,*) 'c',c
!write (*,*) 'u',u
!! JD
!! Select either group (u) or phase (c) velocity
ier = 0
IF(IGRP == 1)THEN
  vel = u
ELSE
  vel = c
ENDIF
RETURN
END SUBROUTINE
