 ! DO ivo = 1,obj%k
  !  CALL GETREF(vref,vpvsref,obj%voro(ivo,1))
   ! partmp(ivo,2) = vref + obj%voro(ivo,2)
   ! partmp(ivo,3) = vpvsref + obj%voro(ivo,3)
   ! !PRINT*,ivo,obj%voro(ivo,1),obj%voro(ivo,2),vref
 ! ENDDO 
!!=======================================================================

SUBROUTINE LOGLHOOD(obj,ipred)
!!=======================================================================
!! Joint log likelihood of the surface-wave dispersion (SWD) and Rayleigh
!! ellipticity (ELL) data sets (independent errors).
USE RJMCMC_COM
IMPLICIT NONE
INCLUDE 'params.h'
INTEGER(KIND=IB):: ipred
TYPE (objstruc) :: obj
REAL(KIND=RP)   :: logL_SWD, logL_ELL
REAL(KIND=SP),DIMENSION(maxlay,10) :: curmod
!!
!! Layer model for the forward codes: curmod = (/ thick, rho, alph, beta, ... /)
!!
curmod = 0._SP
CALL MAKE_CURMOD(obj,curmod)
IF(I_SWD == 1)THEN
  CALL LOGLHOOD_SWD(obj,curmod,ipred,logL_SWD)
ELSE
  logL_SWD = 0._RP
ENDIF
IF(I_ELL == 1)THEN
  CALL LOGLHOOD_ELL(obj,curmod,ipred,logL_ELL)
ELSE
  logL_ELL = 0._RP
ENDIF
obj%logL = logL_SWD + logL_ELL
RETURN
END SUBROUTINE LOGLHOOD
!!=======================================================================


!!=======================================================================


!!=======================================================================


!=======================================================================


!!=======================================================================


!!=======================================================================


!!=======================================================================

SUBROUTINE LOGLHOOD_SWD(obj,curmod,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted SWD data and compute logL_SWD
!!
USE RJMCMC_COM
USE ieee_arithmetic
IMPLICIT NONE
INCLUDE 'params.h'

INTEGER(KIND=IB)                      :: ipred,ierr_swd,imod,ibadlogL,ilay
TYPE (objstruc)                       :: obj
REAL(KIND=SP),DIMENSION(maxlay,10)    :: curmod
REAL(KIND=SP),DIMENSION(maxlay+NPREM,10):: curmod2
REAL(KIND=SP),DIMENSION(NDAT_SWD)     :: periods,DpredSWD
INTEGER(KIND=IB),DIMENSION(NDAT_SWD)  :: ivalidSWD
INTEGER(KIND=IB)                      :: imode,nswd_m
REAL(KIND=RP)                         :: dc_over
REAL(KIND=RP),DIMENSION(NMODE)        :: EtmpSWD
REAL(KIND=RP)                         :: logL,factvs,factvpvs
REAL(KIND=RP)                         :: tstart, tend, tcmp   ! Overall time 
LOGICAL :: ISNAN
REAL(KIND=SP), ALLOCATABLE, DIMENSION(:,:) :: curmodtest

IF(IMAP == 1)THEN
  PRINT*,'CURMOD IN SWD'
  DO ilay=1,obj%nunique+1
    WRITE(*,206)ilay,curmod(ilay,1:10)
  ENDDO
  206   FORMAT(I3,10F12.4)
ENDIF

!!
!!  curmod = thick  rho  alph  beta  %P  %S  tr  pl  st  dip 
!!

!! Find lowest valid entry for perturbations
!factvs   = obj%par(obj%nunique*NPL+1)  !! Vs is in km/s here
!factvpvs = obj%par(obj%nunique*NPL+2)
DO ilay=obj%k,1,-1
  IF(obj%voroidx(ilay,2) == 1)THEN
    factvs   = obj%voro(ilay,2)  !! Vs is in km/s here
    EXIT
  ENDIF
ENDDO
IF(I_VPVS==1) THEN
DO ilay=obj%k,1,-1
  IF(obj%voroidx(ilay,3) == 1)THEN
    factvpvs = obj%voro(ilay,3)
    EXIT
  ENDIF
ENDDO
ELSE
  factvpvs = 0.
END IF
!WRITE(*,*) 'vel_prem = ', SHAPE(vel_prem), vel_prem
curmod2  = 0.
curmod2(1:obj%nunique+1,:) = curmod(1:obj%nunique+1,:)
!! last layer thickness
!curmod2(obj%nunique+1,1)   = (vel_prem(1,1)*1000.)-obj%par((obj%nunique-1)*NPL+1)*1000.
curmod2(obj%nunique+1,1)   = (vel_prem(1,1)*1000.)-obj%voro(obj%k,1)*1000.

!PRINT*,'obj voro'
!DO ilay=1,obj%k+1
!  PRINT*,obj%voro(ilay,:)
!ENDDO
!PRINT*,'obj par',obj%par(1:obj%NFP+2)
!PRINT*,'fact',factvs,factvpvs,obj%par((obj%nunique-1)*NPL+1)

curmod2(obj%nunique+2:obj%nunique+NPREM,1)   = (vel_prem(1,2:NPREM)-vel_prem(1,1:NPREM-1)) * 1000. !! thickness in km
curmod2(obj%nunique+1+NPREM,1)               = 0.                                                  !! HS thickness is 0 
curmod2(obj%nunique+2:obj%nunique+1+NPREM,2) = vel_prem(4,1:NPREM) * 1000.            !! Density
curmod2(obj%nunique+2:obj%nunique+1+NPREM,4) = (vel_prem(2,1:NPREM) + factvs) * 1000. !! Vs
curmod2(obj%nunique+2:obj%nunique+1+NPREM,3) = (vel_prem(2,1:NPREM)+factvs)*(vel_prem(3,1:NPREM)+factvpvs)*1000. !! Vp
IF(IMAP == 1)THEN
  WRITE(*,*) 'curmod2 (including PREM)'
  DO ilay=1,obj%nunique+1+NPREM+1
    WRITE(*,206)ilay,curmod2(ilay,1:10)
  ENDDO
  !206   FORMAT(I3,10F12.4)
ENDIF

!!
!! Adjacent-layer Vs contrast constraint (indicator prior, ABSOLUTE km/s like
!! BayHunter's lvz/hvz). Checked on the final layer stack BEFORE any forward
!! call, so violating proposals cost nothing. Motivation: Kennett (2023, 2026,
!! Seismica) - strong internal contrasts in compaction-gradient structures
!! create spurious internal waveguides; without this the multi-mode posterior
!! can be dominated by a data-fitting but geologically absurd fast-lid family.
!!
IF(DVSCON > 0._RP)THEN
  DO ilay = 1,obj%nunique
    IF(ABS(curmod2(ilay+1,4)-curmod2(ilay,4)) > DVSCON*1000._RP)THEN
      logL = -HUGE(1._RP)
      RETURN
    ENDIF
  ENDDO
ENDIF

!!
!! Forward-model every curve slot on its OWN period grid and point count, as
!! the Rayleigh mode MODE_OF(imode) (0 = fundamental; RAYDSPN counts roots).
!!
dc_over = SWD_DC_OVER
IF(dc_over <= 0._RP) dc_over = SWD_DC/5._RP
DO imode = 1,NMODE
  nswd_m = NDAT_MODE(imode)
  IF(nswd_m <= 0) CYCLE
  periods = 0._SP
  periods(1:nswd_m) = REAL(obj%periods(imode,1:nswd_m),SP)
  !!
  !!  Need to append PREM perturbed by half-space perturbation here to
  !!  ensure that long period SWD can be properly modelled.
  !!
  CALL dispersion(obj%nunique+1+NPREM,curmod2(1:obj%nunique+1+NPREM,2)/1000., &
       curmod2(1:obj%nunique+1+NPREM,3)/1000.,curmod2(1:obj%nunique+1+NPREM,4)/1000.,&
       curmod2(1:obj%nunique+1+NPREM,1)/1000.,DpredSWD,&
       periods,nswd_m,IGRP,ierr_swd,MODE_OF(imode),ivalidSWD,&
       REAL(SWD_CMIN,SP),REAL(SWD_CMAX,SP),REAL(SWD_DC,SP),REAL(dc_over,SP))

  IF(ierr_swd < 0)THEN
    !! Hard input error in the propagator: reject.
    logL = -HUGE(1._RP)
    RETURN
  ENDIF
  !!
  !! Every OBSERVED datum must be predictable. A model that cannot produce a
  !! mode at a period where we actually measured it is inconsistent with the
  !! data, so it is rejected. Only observed periods are ever tested, so a
  !! model is never rejected for failing where there is no datum. Dropping
  !! failed points instead would be WRONG: the number of data would then vary
  !! between models and the likelihood would reward vanishing modes.
  !!
  IF(SUM(ivalidSWD(1:nswd_m)) < nswd_m)THEN
    logL = -HUGE(1._RP)
    RETURN
  ENDIF

  obj%DpredSWD(imode,1:nswd_m) = REAL(DpredSWD(1:nswd_m),RP)
  obj%DresSWD(imode,1:nswd_m)  = obj%DobsSWD(imode,1:nswd_m)-obj%DpredSWD(imode,1:nswd_m)
ENDDO

ibadlogL = 0
IF(IAR == 1)THEN
  obj%DarSWD  = 0._RP
  !!
  !!  Compute autoregressive model
  !!
  DO imode = 1,NMODE
    IF(NDAT_MODE(imode) > 0) CALL ARPRED_SWD(obj,imode,1,NDAT_MODE(imode))
  ENDDO
  !! Recompute predicted data as ith autoregressive model
  obj%DresSWD = obj%DresSWD-obj%DarSWD

  !! Check if predicted AR model data are outside max allowed bounds
  CALL CHECKBOUNDS_ARMXSWD(obj,ibadlogL)
ENDIF

!!
!!  Compute log likelihood
!!
IF(ibadlogL == 0)THEN
  EtmpSWD = 0._RP
  IF(ICOV_SWD == 0)THEN
  !!
  !! implicitly sample over sigma
  !!
  DO imod = 1,NMODE
    EtmpSWD(imod) = -REAL(NDAT_MODE(imod),RP)/2._RP * LOG( SUM(obj%DresSWD(imod,1:NDAT_MODE(imod))**2._RP) / REAL(NDAT_MODE(imod),RP) )
  ENDDO
  ELSEIF(ICOV_SWD == 1)THEN
    !!
    !! Sample over sigma (one per mode)
    !!
    IF (IMAGSCALE==0) THEN
      DO imod = 1,NMODE
        EtmpSWD(imod) = LOG(1._RP/(2._RP*PI2)**(REAL(NDAT_MODE(imod),RP)/2._RP)) &
                      -(SUM(obj%DresSWD(imod,1:NDAT_MODE(imod))**2._RP)/(2._RP*obj%sdparSWD(imod)**2._RP)&
                      +REAL(NDAT_MODE(imod),RP)*LOG(obj%sdparSWD(imod)))
      ENDDO
    ELSE !!IMAGSCALE
      DO imod = 1,NMODE
        EtmpSWD(imod) = -(REAL(NDAT_MODE(imod),RP)/2._RP)*LOG(2._RP*PI2) &
                        -SUM(obj%DresSWD(imod,1:NDAT_MODE(imod))**2._RP/ &
                            (2._RP*(obj%DobsSWD(imod,1:NDAT_MODE(imod))*obj%sdparSWD(imod))**2._RP)) &
                        -REAL(NDAT_MODE(imod),RP)*LOG(obj%sdparSWD(imod)) - &
                         SUM(LOG(ABS(obj%DobsSWD(imod,1:NDAT_MODE(imod)))))
      END DO
    END IF !!IMAGSCALE
  ELSEIF(ICOV_SWD == 2)THEN
    !! Empirical cov mat estimate and magnitude scaling (scaling is done via std
    !! dev parameter... (see Dettmer et al. 2014 GJI)  scaling: Cd = sC'd not Cd = (s^2)C'd 
    DO imod = 1,NMODE
      EtmpSWD(imod) = -DOT_PRODUCT(MATMUL(obj%DresSWD(imod,1:NDAT_MODE(imod)),CdiSWD(1:NDAT_MODE(imod),1:NDAT_MODE(imod))),obj%DresSWD(imod,1:NDAT_MODE(imod))) !TRANSPOSE(obj%DresSWD(imod,1:NDAT_MODE(imod))))
      EtmpSWD(imod) = EtmpSWD(imod)/(2._RP*obj%sdparSWD(imod)) &
           -REAL(NDAT_MODE(imod),RP)/2._RP*LOG(obj%sdparSWD(imod))
       !EtmpSWD(imod) = LOG(1._RP/(2._RP*PI2)**(REAL(NDAT_MODE(imod),RP)/2._RP)) &
       !              -(SUM(obj%DresSWD(imod,1:NDAT_MODE(imod))**2._RP)/(2._RP*obj%sdparSWD(imod)**2._RP)&
       !              +REAL(NDAT_MODE(imod),RP)*LOG(obj%sdparSWD(imod)))
    ENDDO
  ELSEIF(ICOV_SWD ==3)THEN
    DO imod = 1,NMODE
      EtmpSWD(imod) = -REAL(NDAT_MODE(imod),RP)/2._RP*LOG(2._RP*PI2) &
                      -REAL(NDAT_MODE(imod),RP)*LOG(obj%sdparSWD(imod)) - SUM(LOG(sdSWD(imod,1:NDAT_MODE(imod)))) &
                      -SUM((obj%DresSWD(imod,1:NDAT_MODE(imod))/sdSWD(imod,1:NDAT_MODE(imod)))**2)/(2._RP*obj%sdparSWD(imod)**2)
    END DO
  ENDIF
  logL = SUM(EtmpSWD)
  IF(ieee_is_nan(logL))THEN
    logL = -HUGE(1._RP)
  ENDIF
ELSE
   logL = -HUGE(1._RP)
   ibadlogL = 0
ENDIF
IF (IMAP == 1) WRITE(*,*) 'logL_SWD = ', logL
RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD_SWD

!!=======================================================================

SUBROUTINE LOGLHOOD_ELL(obj,curmod,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted ELL data and compute logL_ELL
!!
USE RJMCMC_COM
USE ieee_arithmetic
IMPLICIT NONE
INCLUDE 'params.h'

INTEGER(KIND=IB)                      :: ipred,ierr_ELL,imod,ibadlogL,ilay, iidx
TYPE (objstruc)                       :: obj
REAL(KIND=SP),DIMENSION(maxlay,10)    :: curmod
REAL(KIND=SP),DIMENSION(maxlay+NPREM,10):: curmod2
REAL(KIND=RP),DIMENSION(NDAT_ELL)     :: periods,angFres,DpredELL
REAL(KIND=RP),DIMENSION(NMODE_ELL)    :: EtmpELL
REAL(KIND=RP)                         :: logL,factvs,factvpvs
REAL(KIND=RP)                         :: tstart, tend, tcmp   ! Overall time 
LOGICAL :: ISNAN
REAL(KIND=RP), ALLOCATABLE, DIMENSION(:,:) :: curmodtest
INTEGER(KIND=IB)                      :: nlayers, NDAT_ELL_tmp
REAL(KIND=RP), DIMENSION(maxlay)      :: thicknesses, vp, vs, rho

IF(IMAP == 1)THEN
  PRINT*,'CURMOD IN ELL'
  DO ilay=1,obj%nunique+1
    WRITE(*,206)ilay,curmod(ilay,1:10)
  ENDDO
  206   FORMAT(I3,10F12.4)
ENDIF

!!
!!  curmod = thick  rho  alph  beta  %P  %S  tr  pl  st  dip 
!!

!! Find lowest valid entry for perturbations
!factvs   = obj%par(obj%nunique*NPL+1)  !! Vs is in km/s here
!factvpvs = obj%par(obj%nunique*NPL+2)
DO ilay=obj%k,1,-1
  IF(obj%voroidx(ilay,2) == 1)THEN
    factvs   = obj%voro(ilay,2)  !! Vs is in km/s here
    EXIT
  ENDIF
ENDDO
IF (I_VPVS==1) THEN
DO ilay=obj%k,1,-1
  IF(obj%voroidx(ilay,3) == 1)THEN
    factvpvs = obj%voro(ilay,3)
    EXIT
  ENDIF
ENDDO
ELSE
factvpvs = 0.
END IF
!WRITE(*,*) 'vel_prem = ', SHAPE(vel_prem), vel_prem
curmod2  = 0.
curmod2(1:obj%nunique+1,:) = curmod(1:obj%nunique+1,:)
!! last layer thickness
!curmod2(obj%nunique+1,1)   = (vel_prem(1,1)*1000.)-obj%par((obj%nunique-1)*NPL+1)*1000.
curmod2(obj%nunique+1,1)   = (vel_prem(1,1)*1000.)-obj%voro(obj%k,1)*1000.

!PRINT*,'obj voro'
!DO ilay=1,obj%k+1
!  PRINT*,obj%voro(ilay,:)
!ENDDO
!PRINT*,'obj par',obj%par(1:obj%NFP+2)
!PRINT*,'fact',factvs,factvpvs,obj%par((obj%nunique-1)*NPL+1)

curmod2(obj%nunique+2:obj%nunique+NPREM,1)   = (vel_prem(1,2:NPREM)-vel_prem(1,1:NPREM-1)) * 1000. !! thickness in km
curmod2(obj%nunique+1+NPREM,1)               = 0.                                                  !! HS thickness is 0 
curmod2(obj%nunique+2:obj%nunique+1+NPREM,2) = vel_prem(4,1:NPREM) * 1000.            !! Density
curmod2(obj%nunique+2:obj%nunique+1+NPREM,4) = (vel_prem(2,1:NPREM) + factvs) * 1000. !! Vs
curmod2(obj%nunique+2:obj%nunique+1+NPREM,3) = (vel_prem(2,1:NPREM)+factvs)*(vel_prem(3,1:NPREM)+factvpvs)*1000. !! Vp
IF(IMAP == 1)THEN
  WRITE(*,*) 'curmod2 (including PREM)'
  DO ilay=1,obj%nunique+1+NPREM+1
    WRITE(*,206)ilay,curmod2(ilay,1:10)
  ENDDO
  !206   FORMAT(I3,10F12.4)
ENDIF

IF(IMAP == 1)THEN
  PRINT*, 'PAR:'
  WRITE(*,*) obj%par(1:NPL*obj%k)
END IF

periods = REAL(obj%periods_ELL(1,1:NDAT_ELL), RP)
!WRITE(*,*) SHAPE(periods)
!OPEN(UNIT=10, FILE='periods', STATUS='REPLACE', ACTION='WRITE')
!WRITE(10,*) periods
!CLOSE(10)
nlayers = obj%nunique+1  !! if use curmod
!WRITE(*,*) 'nlayers ', nlayers
thicknesses(1:nlayers)  = REAL(curmod(1:nlayers,1), RP)
!WRITE(*,*) 'h = ', thicknesses(1:nlayers)
vp(1:nlayers) = REAL(curmod(1:nlayers,3), RP)
!WRITE(*,*) 'vp = ', vp(1:nlayers)
vs(1:nlayers) = REAL(curmod(1:nlayers,4), RP)
!WRITE(*,*) 'vs = ', vs(1:nlayers)
rho(1:nlayers) =REAL(curmod(1:nlayers,2), RP)
!WRITE(*,*) 'rho = ', rho(1:nlayers)
!WRITE(*,*) 'prec = ', ELL_prec
!!
!!  Need to append PREM perturbed by half-space perturbation here to 
!!  ensure that long period SWD can be properly modelled. 
!!
!nlayers = obj%nunique+1+NPREM
!thicknesses(1:nlayers)  =  REAL(curmod2(1:nlayers,1), RP)
!vp(1:nlayers) = REAL(curmod2(1:nlayers,3), RP)
!vs(1:nlayers) = REAL(curmod2(1:nlayers,4), RP)
!rho(1:nlayers) = REAL(curmod2(1:nlayers,2), RP)

!PRINT*, 'before calling gpell'

!CALL ellipticity_gpell(nlayers,thicknesses(1:nlayers),&
!                       vp(1:nlayers),vs(1:nlayers),rho(1:nlayers),&
!                       NDAT_ELL,periods,NMODE_ELL,DpredELL,I_ABS_ELL,I_LOG10_ELL,&
!                       ELL_prec,I_SAMPLING_TYPE_ELL, I_SET_STEP_ELL, STEP_SIZE_ELL,&
!                       I_SET_COUNT_ELL, COUNT_ELL,&
!                       I_SET_RANGE_ELL, ierr_ELL)


IF(ierr_ELL /= 0)THEN
!!  WRITE(*,*)'THIS MODEL IS WEIRD, Cannot compute ellipticity'
!!  WRITE(*,*)'ierr_ELL = ', ierr_ELL
  
!  !stop
!  !CALL PRINTPAR(obj)
  logL = -HUGE(1._RP)
  IF (IMAP == 0) WRITE(77,*) obj%par
  RETURN
ENDIF

!DO iidx = 1, NDAT_ELL
!    IF (DpredELL(iidx) == 0._RP) THEN
!        WRITE(*,*)'THIS MODEL IS WEIRD, Cannot compute ellipticity'
!        logL = -HUGE(1._RP)
!        WRITE(77,*) obj%par
!        RETURN
!    END IF
!END DO

obj%DpredELL(1,1:NDAT_ELL) = DpredELL(NDAT_ELL:1:-1)
obj%DresELL(1,1:NDAT_ELL) = obj%DobsELL(1,1:NDAT_ELL)-obj%DpredELL(1,1:NDAT_ELL)
!IF (I_ABS_ELL == 0) THEN
!    obj%DresELL(1,1:NDAT_ELL) = obj%DobsELL(1,1:NDAT_ELL)-obj%DpredELL(1,1:NDAT_ELL)
!ELSEIF (I_LOG10_ELL == 0) THEN
!    obj%DresELL(1,1:NDAT_ELL) = ABS(obj%DobsELL(1,1:NDAT_ELL))-ABS(obj%DpredELL(1,1:NDAT_ELL))
!ELSE
!    obj%DresELL(1,1:NDAT_ELL) = LOG10(ABS(obj%DobsELL(1,1:NDAT_ELL)))-LOG10(ABS(obj%DpredELL(1,1:NDAT_ELL)))
!END IF

ibadlogL = 0
IF(IAR == 1)THEN
  obj%DarELL  = 0._RP
  !!
  !!  Compute autoregressive model
  !!
  CALL ARPRED_ELL(obj,1,1,NDAT_ELL)
  !! Recompute predicted data as ith autoregressive model
  obj%DresELL = obj%DresELL-obj%DarELL

  !! Check if predicted AR model data are outside max allowed bounds
  CALL CHECKBOUNDS_ARMXELL(obj,ibadlogL)
ENDIF

!!
!!  Compute log likelihood
!!
IF(ibadlogL == 0)THEN
  EtmpELL = 0._RP
  IF(ICOV_ELL == 0)THEN
  !!
  !! implicitly sample over sigma
  !!
  DO imod = 1,NMODE_ELL
    EtmpELL(imod) = -REAL(NDAT_ELL,RP)/2._RP * LOG( SUM(obj%DresELL(imod,:)**2._RP) / REAL(NDAT_ELL,RP) )
  ENDDO
  ELSEIF(ICOV_ELL == 1)THEN
    !!
    !! Sample over sigma (one per mode)
    !!
    !NDAT_ELL_tmp = 81
    NDAT_ELL_tmp = NDAT_ELL
    DO imod = 1,NMODE_ELL
      EtmpELL(imod) = LOG(1._RP/(2._RP*PI2)**(REAL(NDAT_ELL_tmp,RP)/2._RP)) &
                    -(SUM(obj%DresELL(imod,:)**2._RP)/(2._RP*obj%sdparELL(imod)**2._RP)&
                    +REAL(NDAT_ELL_tmp,RP)*LOG(obj%sdparELL(imod)))
    ENDDO
  ELSEIF(ICOV_ELL == 2)THEN
    !! Empirical cov mat estimate and magnitude scaling (scaling is done via std
    !! dev parameter... (see Dettmer et al. 2014 GJI)  scaling: Cd = sC'd not Cd = (s^2)C'd 
    DO imod = 1,NMODE_ELL
      EtmpELL(imod) = -DOT_PRODUCT(MATMUL(obj%DresELL(imod,:),CdiELL),obj%DresELL(imod,:)) !TRANSPOSE(obj%DresSWD(imod,:)))
      EtmpELL(imod) = EtmpELL(imod)/(2._RP*obj%sdparELL(imod)) &
           -REAL(NDAT_ELL,RP)/2._RP*LOG(obj%sdparELL(imod))
      !EtmpELL(imod) = LOG(1._RP/(2._RP*PI2)**(REAL(NDAT_ELL,RP)/2._RP)) &
      !              -(SUM(obj%DresELL(imod,:)**2._RP)/(2._RP*obj%sdparELL(imod)**2._RP)&
      !              +REAL(NDAT_ELL,RP)*LOG(obj%sdparELL(imod)))
    ENDDO
  ENDIF
  logL = SUM(EtmpELL)
  IF(ieee_is_nan(logL))THEN
    logL = -HUGE(1._RP)
  ENDIF
ELSE
   logL = -HUGE(1._RP)
   ibadlogL = 0
ENDIF
IF (IMAP == 1) WRITE(*,*) 'logL_ELL = ', logL
RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD_ELL

!!=======================================================================



!!=======================================================================



!!=======================================================================

SUBROUTINE MAKE_CURMOD(obj,curmod)
!!=======================================================================
!!
!! Build curmod from model
!!
USE RJMCMC_COM
USE ieee_arithmetic
IMPLICIT NONE
INCLUDE 'params.h'

INTEGER(KIND=IB)                 :: ipar,ilay,id,iparcur, NPL2, ijump
TYPE (objstruc)                  :: obj
REAL(KIND=SP),DIMENSION(maxlay,10):: curmod

NPL2 = NPL
ijump = 0

curmod = 0._RP
!! curmod = (/ thick, rho, alph, beta, %P, %S, tr, pl, st, di /)
id=0
DO ilay=1,obj%nunique+1
  !curmod(ilay,:) = curmod_glob
  DO ipar=1,NPL2
    id=id+1
    iparcur = idxpar(ipar)
    IF(ilay <= obj%nunique)THEN
      IF(iparcur == 1)THEN
        curmod(ilay,iparcur)    = REAL(obj%hiface(ilay),SP)  !! Model works with layer thickness!
      ELSE
        curmod(ilay,iparcur)    = REAL(obj%par(id),SP)
      ENDIF
    ELSE
      IF(iparcur == 1)THEN
        curmod(ilay,iparcur) = 0._SP
        id=id-1
      ELSE
        curmod(ilay,iparcur)    = REAL(obj%par(id),SP)
      ENDIF
    ENDIF
  ENDDO
ENDDO
IF(I_VPVS == 1)THEN
  !! for sampling Vs:
  curmod(1:obj%nunique+1,1) = curmod(1:obj%nunique+1,1)*1000._SP
  curmod(1:obj%nunique+1,4) = curmod(1:obj%nunique+1,4)*1000._SP
  curmod(1:obj%nunique+1,3) = curmod(1:obj%nunique+1,3)*curmod(1:obj%nunique+1,4)
ELSEIF(I_VPVS == -1)THEN
  !! for sampling Vs:
  curmod(1:obj%nunique+1,1) = curmod(1:obj%nunique+1,1)*1000._SP
  curmod(1:obj%nunique+1,4) = curmod(1:obj%nunique+1,4)*1000._SP
  curmod(1:obj%nunique+1,3) = curmod(1:obj%nunique+1,4)*1.75_SP
ELSE
  curmod(1:obj%nunique+1,1:4) = curmod(1:obj%nunique+1,1:4)*1000._SP
ENDIF

!! Birch Law for analytical density relationship
!curmod(1:obj%nunique+1,2)   = 1000.*0.77+0.32*curmod(:,3)
!! What Thomas uses:
curmod(1:obj%nunique+1,2) = (2.35+0.036*((curmod(1:obj%nunique+1,3)/1000.)-3.0)**2.)*1000.

!!
!! Convert angles to radians (trend, plunge, strike, dip)
curmod(1:obj%nunique+1,7)  = curmod(1:obj%nunique+1,7)/180._RP * PI
curmod(1:obj%nunique+1,8)  = curmod(1:obj%nunique+1,8)/180._RP * PI
curmod(1:obj%nunique+1,9)  = curmod(1:obj%nunique+1,9)/180._RP * PI
curmod(1:obj%nunique+1,10) = curmod(1:obj%nunique+1,10)/180._RP * PI

!IF(IMAP == 1)THEN
!  DO ilay=1,obj%nunique+1
!    WRITE(*,206)ilay,curmod(ilay,1:10)
!  ENDDO
!  206   FORMAT(I3,10F12.4)
!ENDIF
RETURN
END SUBROUTINE MAKE_CURMOD
!!=======================================================================


!!=======================================================================

SUBROUTINE INTERPLAYER_novar(obj)
!!=======================================================================
!!
!! This interpolates 1D layer nodes onto obj%par array for forward model
!! This does not allow for variable layer complexity.
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB) :: ivo,ivo2,ipar,ilay,iface,itmp,ntot,ilim
INTEGER(KIND=IB) :: NPL_tmp
TYPE(objstruc) :: obj
REAL(KIND=RP),DIMENSION(NPL*NLMX,NPL):: partmp
REAL(KIND=RP),DIMENSION(NLMX,2)  :: vorotmp
REAL(KIND=RP)                    :: vref,vpvsref
REAL(KIND=DRP),DIMENSION(NLMX,NPL):: tmpsort

!  obj%k       = 3
!  obj%NFP     = (obj%k * NPL) + (NPL-1)
!  obj%voro    = 0._RP
!  obj%voroidx = 0
!  !! True parameters for simulation:
!  obj%voro(1,:) = (/  0.0_RP, 6.0_RP, 1.8_RP,  0.0_RP/)
!  obj%voro(2,:) = (/ 10.0_RP, 7.0_RP, 0.0_RP,  5.0_RP/)
!  obj%voro(3,:) = (/ 30.0_RP, 8.0_RP, 0.0_RP, 30.0_RP/)
!  obj%voroidx(1,:) = (/ 1, 1, 1, 0/)
!  obj%voroidx(2,:) = (/ 1, 1, 0, 1/)
!  obj%voroidx(3,:) = (/ 1, 1, 0, 1/)

!PRINT*,''
!PRINT*,''
!PRINT*,'Before:',rank
!DO ilay = 1,obj%k
!  WRITE(*,203)ilay,obj%voro(ilay,:)
!ENDDO
!PRINT*,''
!PRINT*,''

!! Sort node according to increasing depth:
obj%nunique = obj%k-1

tmpsort = 0._DRP
tmpsort(1:obj%k,:) = REAL(obj%voro(1:obj%k,:),DRP)
CALL QSORTC2D(tmpsort(1:obj%k,:),obj%voroidx(1:obj%k,:))
obj%voro(1:obj%k,:) = REAL(tmpsort(1:obj%k,:),RP)

obj%ziface = 0._RP
obj%ziface(1:obj%k-1) = obj%voro(2:obj%k,1)

obj%hiface = 0._RP
obj%hiface(1) = obj%ziface(1)
DO ivo = 2,obj%k-1
  obj%hiface(ivo) = obj%ziface(ivo)-obj%ziface(ivo-1)
ENDDO

partmp = 0._RP
partmp(1:obj%k,1) = obj%ziface(1:obj%k)
partmp(1:obj%k,2:NPL) = obj%voro(1:obj%k,2:NPL)

IF(I_VREF == 1)THEN
  DO ivo = 1,obj%k
    !!WRITE(*,*) 'nlayer,z = ', ivo, obj%voro(ivo,1)
    CALL GETREF(vref,vpvsref,obj%voro(ivo,1))
    partmp(ivo,2) = vref + partmp(ivo,2)      
    IF(I_VPVS==1) partmp(ivo,3) = vpvsref + partmp(ivo,3)   !!   : no third column
    !PRINT*,ivo,obj%voro(ivo,1),obj%voro(ivo,2),vref
  ENDDO 
END IF
!!
!! Apply reference profile:
!!
!IF(I_VREF == 0)THEN
 ! partmp(1:obj%k,2:NPL) = obj%voro(1:obj%k,2:NPL)
!ELSE
 ! DO ivo = 1,obj%k
  !  CALL GETREF(vref,vpvsref,obj%voro(ivo,1))
   ! partmp(ivo,2) = vref + obj%voro(ivo,2)
   ! partmp(ivo,3) = vpvsref + obj%voro(ivo,3)
   ! !PRINT*,ivo,obj%voro(ivo,1),obj%voro(ivo,2),vref
 ! ENDDO 
!ENDIF

obj%par = 0._RP
DO ilay = 1,obj%k-1
  obj%par((ilay-1)*NPL+1:ilay*NPL) = partmp(ilay,:)
ENDDO
obj%par((obj%k-1)*NPL+1:(obj%k*NPL)-1) = partmp(ilay,2:)

203 FORMAT(I4,20F8.2)

END SUBROUTINE INTERPLAYER_novar
!!=======================================================================

SUBROUTINE INTERPLAYER(obj)
!!=======================================================================
!!
!! This interpolates 1D layer nodes onto obj%par array for forward model
!! I.e., some layers are duplicates when nodes are not populated:
!!  Parameter 1:  Parameter 2:     Parameter 3:
!!   --o------      ---o---       ------o--------
!!     |               |                |
!!   ------o--      -------       --o------------
!!         |           |            |
!!         |           |            |
!!         |           |            |
!!   ---------      -------       ------------o--
!!         |           |                      |
!!   ---o-----      -------       ---------------
!!      |              |                      |
!!      |              |                      |
!!   ---------      -------       ---------------
!!
!! The layer node always defines the volume partition below the 
!! node position. The node position defines the interface.
!! The first layer is always fixed at 0 and all nodes are populated.
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB) :: ivo,ivo2,ipar,ilay,iface,itmp,ntot
INTEGER(KIND=IB) :: NPL_tmp
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION(NPL-1)  :: niface
INTEGER(KIND=RP),DIMENSION(NPL*NLMX):: ifaceidx
REAL(KIND=RP),DIMENSION(NPL*NLMX,NPL):: partmp
REAL(KIND=RP),DIMENSION(NLMX,2)  :: vorotmp
REAL(KIND=RP),DIMENSION(NLMX,NPL-1)  :: voroh
REAL(KIND=RP),DIMENSION(NLMX+1,NPL-1,2):: voroz

REAL(KIND=RP),DIMENSION(NLMX,NPL)::  voro(NLMX,NPL),voroidx(NLMX,NPL)
REAL(KIND=RP)                    :: vref,vpvsref
REAL(KIND=DRP),DIMENSION(NLMX,NPL):: tmpsort
REAL(KIND=DRP),DIMENSION(NLMX*NPL):: tmpsort2

!  obj%k       = 3
!  obj%NFP     = (obj%k * NPL) + (NPL-1)
!  obj%voro    = 0._RP
!  obj%voroidx = 0
!  !! True parameters for simulation:
!  obj%voro(1,:) = (/  0.0_RP, 6.0_RP, 1.8_RP,  0.0_RP/)
!  obj%voro(2,:) = (/ 10.0_RP, 7.0_RP, 0.0_RP,  5.0_RP/)
!  obj%voro(3,:) = (/ 30.0_RP, 8.0_RP, 0.0_RP, 30.0_RP/)
!  obj%voroidx(1,:) = (/ 1, 1, 1, 0/)
!  obj%voroidx(2,:) = (/ 1, 1, 0, 1/)
!  obj%voroidx(3,:) = (/ 1, 1, 0, 1/)

!PRINT*,''
!PRINT*,''
!PRINT*,'Before:',rank
!DO ilay = 1,obj%k
!  WRITE(*,203)ilay,obj%voro(ilay,:)
!ENDDO
!PRINT*,''
!PRINT*,''

!obj%voroidx(1,:) = (/ 1, 1, 1/)
!obj%voroidx(2,:) = (/ 1, 1, 0/)
!obj%voroidx(3,:) = (/ 1, 1, 0/)

IF(IDIP == 1)THEN
  NPL_tmp = NPL-1
ELSE
  NPL_tmp = NPL
ENDIF

tmpsort = 0._DRP
tmpsort(1:obj%k,:) = REAL(obj%voro(1:obj%k,:),DRP)
CALL QSORTC2D(tmpsort(1:obj%k,:),obj%voroidx(1:obj%k,:))
obj%voro(1:obj%k,:) = REAL(tmpsort(1:obj%k,:),RP)

!!
!! Use local variable for voro and voroidx to allow easy change
!! from perturbation value to perturbation+vref
!!
voro    = 0._RP
voroidx = 0
voro    = obj%voro
voroidx = obj%voroidx

!!
!! Apply reference profile:
!!
IF(I_VREF == 1)THEN
  DO ivo = 1,obj%k
    CALL GETREF(vref,vpvsref,voro(ivo,1))
    !WRITE(*,201)ivo,voro(ivo,1),voro(ivo,2),voro(ivo,3)
    IF(voroidx(ivo,2) == 1)voro(ivo,2) = vref + voro(ivo,2)
    IF(I_VPVS==1) THEN 
      IF(voroidx(ivo,3) == 1)voro(ivo,3) = vpvsref + voro(ivo,3)
    END IF
    !WRITE(*,201)ivo,voro(ivo,1),voro(ivo,2),voro(ivo,3),vref,vpvsref
    !WRITE(*,*)''
  ENDDO
ENDIF
201 FORMAT(I,5F12.6)
!!
!!  Find interfaces for each parameter
!!
niface = 0._RP
DO ipar = 2,NPL
  niface(ipar-1) = SUM(voroidx(:,ipar))-1
ENDDO
obj%ziface = 0._RP
voroz  = 0._RP
voroh  = 0._RP
iface  = 0
DO ipar = 2,NPL
  vorotmp = 0._RP
  itmp = 0
  DO ivo = 1,obj%k
    IF(voroidx(ivo,ipar) == 1)THEN
      itmp = itmp + 1
      vorotmp(itmp,:) = (/voro(ivo,1),voro(ivo,ipar)/)
    ENDIF
  ENDDO
  IF(niface(ipar-1) > 0)THEN
    DO ivo = 1,niface(ipar-1)
      iface = iface + 1
      obj%ziface(iface) = vorotmp(ivo+1,1)
      voroz(ivo,ipar-1,1) = obj%ziface(iface)
      voroz(ivo,ipar-1,2) = vorotmp(ivo,2)
    ENDDO
    voroz(itmp,ipar-1,2) = vorotmp(itmp,2)
  ELSE
    voroz(1,ipar-1,2) = vorotmp(1,2)
  ENDIF
ENDDO

ntot = SUM(niface)
tmpsort2 = 0._DRP
tmpsort2(1:ntot) = REAL(obj%ziface(1:ntot),DRP)
CALL QSORTC1D(tmpsort2(1:ntot))
obj%ziface(1:ntot) = REAL(tmpsort2(1:ntot),RP)

!!
!!  Find unique interfaces
!!
obj%nunique = obj%k-1
IF(ntot > 1)THEN
  DO ivo = 1,ntot
    itmp = 0
    ifaceidx = 0
    DO ivo2 = ivo+1,ntot
      IF(obj%ziface(ivo) /= obj%ziface(ivo2))THEN
        itmp = itmp + 1
        ifaceidx(itmp) = ivo2
      ENDIF
    ENDDO
    IF(itmp > 0)THEN
      obj%ziface(ivo+1:ivo+itmp) = obj%ziface(ifaceidx(1:itmp))
      obj%ziface(ivo+itmp+1:) = 0._RP
    ENDIF
  ENDDO
ENDIF

obj%hiface = 0._RP
obj%hiface(1) = obj%ziface(1)
DO ivo = 2,obj%nunique
  obj%hiface(ivo) = obj%ziface(ivo)-obj%ziface(ivo-1)
ENDDO

partmp = 0._RP
partmp(1:obj%nunique,1) = obj%ziface(1:obj%nunique)
DO ipar = 2,NPL_tmp
  ivo = 1
  IF(niface(ipar-1) /= 0)THEN
    DO ilay = 1,obj%nunique
      IF(obj%ziface(ilay) <= voroz(ivo,ipar-1,1))THEN
        partmp(ilay,ipar) = voroz(ivo,ipar-1,2)
      ELSE
        ivo = ivo + 1
        partmp(ilay,ipar) = voroz(ivo,ipar-1,2)
      ENDIF
      IF(obj%ziface(ilay) >= voroz(niface(ipar-1),ipar-1,1))THEN
        partmp(ilay+1:obj%nunique+1,ipar) = voroz(ivo+1,ipar-1,2)
        EXIT
      ENDIF
    ENDDO
  ELSE
    partmp(1:obj%nunique+1,ipar) = voroz(1,ipar-1,2)
  ENDIF
ENDDO
partmp(1:obj%nunique,1) = obj%hiface(1:obj%nunique)
IF(IDIP == 1)THEN
  ipar = NPL
  ivo = 1
  partmp(:,ipar) = 0._RP
  IF(SUM(voroidx(:,ipar)) > 0)THEN
    DO ilay = 2,obj%nunique+1
      IF(voroidx(ilay,ipar) == 1)THEN
        partmp(ilay,ipar) = voro(ilay,ipar)
      ENDIF
    ENDDO
  ENDIF
ENDIF

obj%par = 0._RP
DO ilay = 1,obj%nunique
  obj%par((ilay-1)*NPL+1:ilay*NPL) = partmp(ilay,:)
ENDDO
obj%par(obj%nunique*NPL+1:((obj%nunique+1)*NPL)-1) = partmp(ilay,2:)

!IF(obj%voro(1,1)<0._RP)THEN
!PRINT*,''
!PRINT*,'After:',rank
!DO ilay = 1,obj%nunique+1
!  WRITE(*,203)ilay,partmp(ilay,:)
!ENDDO
!CALL PRINTPAR(obj)
!PRINT*,''
!PRINT*,obj%par
!PRINT*,''
!STOP
!PRINT*,'ntot',ntot
!PRINT*,'unique',obj%nunique
!PRINT*,'niface',niface
!PRINT*,'ziface',obj%ziface(1:ntot)
!DO ilay = 1,obj%nunique+1
!  WRITE(*,203)ilay,obj%par((ilay-1)*NPL+1:ilay*NPL)
!ENDDO
!PRINT*,''
!PRINT*,''
!  PRINT*,'INTERP'
!  STOP
!ENDIF
!
!STOP
!200 FORMAT(2I4,20F8.2)
!201 FORMAT(a,20F8.2)
!202 FORMAT(a,20I4)
203 FORMAT(I4,20F8.2)

END SUBROUTINE INTERPLAYER
!=======================================================================

SUBROUTINE GETREF(vref,vpvsref,z)
!!==============================================================================
!!
!! Reads observed data.
!!
USE RJMCMC_COM
IMPLICIT NONE
REAL(KIND=RP)   :: z,vref,vpvsref,grad,dz
INTEGER(KIND=IB):: ipar,iint
IF(z >= vel_ref(1,NVELREF))THEN
  vref = vel_ref(2,NVELREF)
  vpvsref = vel_ref(3,NVELREF)
ELSEIF(z == 0.)THEN
  vref = vel_ref(2,1)
  vpvsref = vel_ref(3,1)
ELSE
  iint = 0
  DO ipar=1,NVELREF
    IF((z-vel_ref(1,ipar)) <= 0.)THEN
      EXIT
    ENDIF
    iint = iint + 1
  ENDDO
  grad = (vel_ref(2,iint+1)-vel_ref(2,iint))/(vel_ref(1,iint+1)-vel_ref(1,iint))
  dz = (z-vel_ref(1,iint))
  vref = vel_ref(2,iint) + dz*grad
  grad = (vel_ref(3,iint+1)-vel_ref(3,iint))/(vel_ref(1,iint+1)-vel_ref(1,iint))
  vpvsref = vel_ref(3,iint) + dz*grad
ENDIF
RETURN
END SUBROUTINE GETREF
!!=======================================================================


!!=======================================================================

SUBROUTINE ARPRED_SWD(obj,ifr,idata,idatb)
!=======================================================================
!!
!! Autoregressive model to model data error correlations.
!! AR process is computed forward and backward and average is used.
!!
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj
INTEGER          :: i,j,k,ifr,idata,idatb
REAL(KIND=RP),DIMENSION(idatb-idata+1)::dres1,dar1
IF(obj%idxarSWD(ifr) == 1)THEN
   k = 1
   obj%DarSWD(ifr,idata)=0._RP          ! Matlab sets first point to zero...

   !!
   !! Real part:
   !!
   dres1 = 0._RP
   dres1 = obj%DresSWD(ifr,idata:idatb)

   dar1(1)=0._RP          ! Matlab sets first point to zero...
   DO i=2,idatb-idata+1
      dar1(i) = 0
      IF(k >= i)THEN
         DO j=1,i-1
            dar1(i) = dar1(i) + obj%arparSWD((ifr-1)+j) * dres1(i-j)
         ENDDO
      ELSE
         DO j=1,k
            dar1(i) = dar1(i) + obj%arparSWD((ifr-1)+j) * dres1(i-j)
         ENDDO
      ENDIF
   ENDDO
   obj%DarSWD(ifr,idata:idatb) = dar1
   obj%DarSWD(ifr,idata) = 0._RP
   obj%DarSWD(ifr,idatb) = 0._RP
ENDIF
END SUBROUTINE ARPRED_SWD
!!=======================================================================

SUBROUTINE ARPRED_ELL(obj,ifr,idata,idatb)
!=======================================================================
!!
!! Autoregressive model to model data error correlations.
!! AR process is computed forward and backward and average is used.
!!
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj
INTEGER          :: i,j,k,ifr,idata,idatb
REAL(KIND=RP),DIMENSION(idatb-idata+1)::dres1,dar1
IF(obj%idxarELL(ifr) == 1)THEN
   k = 1
   obj%DarELL(ifr,idata)=0._RP          ! Matlab sets first point to zero...

   !!
   !! Real part:
   !!
   dres1 = 0._RP
   dres1 = obj%DresELL(ifr,idata:idatb)

   dar1(1)=0._RP          ! Matlab sets first point to zero...
   DO i=2,idatb-idata+1
      dar1(i) = 0
      IF(k >= i)THEN
         DO j=1,i-1
            dar1(i) = dar1(i) + obj%arparELL((ifr-1)+j) * dres1(i-j)
         ENDDO
      ELSE
         DO j=1,k
            dar1(i) = dar1(i) + obj%arparELL((ifr-1)+j) * dres1(i-j)
         ENDDO
      ENDIF
   ENDDO
   obj%DarELL(ifr,idata:idatb) = dar1
   obj%DarELL(ifr,idata) = 0._RP
   obj%DarELL(ifr,idatb) = 0._RP
ENDIF
END SUBROUTINE ARPRED_ELL
!!=======================================================================


!!=======================================================================

SUBROUTINE CHECKBOUNDS_ARMXSWD(obj,ibadlogL)
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ifr,ibadlogL
TYPE(objstruc):: obj

DO ifr = 1,NMODE
   IF(MAXVAL(obj%DarSWD(ifr,:)) > armxSWD)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
   IF(MINVAL(obj%DarSWD(ifr,:)) < -armxSWD)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
ENDDO

RETURN
END SUBROUTINE CHECKBOUNDS_ARMXSWD
!!=======================================================================

SUBROUTINE CHECKBOUNDS_ARMXELL(obj,ibadlogL)
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ifr,ibadlogL
TYPE(objstruc):: obj

DO ifr = 1,NMODE_ELL
   IF(MAXVAL(obj%DarELL(ifr,:)) > armxELL)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
   IF(MINVAL(obj%DarELL(ifr,:)) < -armxELL)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
ENDDO

RETURN
END SUBROUTINE CHECKBOUNDS_ARMXELL
!=======================================================================

SUBROUTINE LOGLHOOD2(obj)
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj

!!
!!  Compute log likelihood
!!
obj%logL = 1._RP

RETURN
END SUBROUTINE LOGLHOOD2
!!=======================================================================
!!EOF
