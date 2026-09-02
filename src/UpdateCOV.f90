SUBROUTINE UpdateCOV(obj)

USE RJMCMC_COM
IMPLICIT NONE

TYPE(objstruc) :: obj   !!INTENT(INOUT)

REAL(KIND=RP) :: Ferro_norm
INTEGER(KIND=IB) :: is, ie, is2, ie2, ierror                      


IF ( cov_converged_datasets(1)==0 ) THEN
    
    IF (ICOV_SWD==2) THEN
        CdSWD_old = CdSWD
        CdiSWD_old = CdiSWD
    END IF
    is = 1
    ie = NMODE2*NDAT_SWD
    CALL COVmatEst(CdSWD, CdiSWD, sampleDres(:,is:ie), sampleDres(:,ncount3), NsampleDres, NDAT_SWD, nfrac_SWD, MAX_NAVE_SWD, damp_power_SWD, inonstat_SWD, iunbiased_SWD, imr_SWD, ierror)
    IF (ierror/=0) THEN
        WRITE(*,*) 'non positive definite SWD cov estimate at iter: ',icovIter
        WRITE(*,*) 'algorithm continues with previous cov matrix'
        IF (ICOV_SWD==2) THEN
            CdSWD = CdSWD_old       !!!The best approach is to use pointers instead of copying arrays
            CdiSWD = CdiSWD_old
        END IF
    ELSE IF (ICOV_SWD==2) THEN
        IF (iconverge_criterion<=2) THEN
            IF (iconverge_criterion_SWD==1) THEN
                covIter_errSWD = MAXVAL(ABS(CdSWD-CdSWD_old))  / MAXVAL(ABS(CdSWD_old))
            ELSEIF (iconverge_criterion_SWD==2) THEN            
                covIter_errSWD =  MAXVAL(ABS(CdSWD-CdSWD_old)) 
            ELSEIF (iconverge_criterion_SWD==3) THEN            
                covIter_errSWD =  Ferro_norm(CdSWD-CdSWD_old,NDAT_SWD) / Ferro_norm(CdSWD_old,NDAT_SWD)
            ELSEIF (iconverge_criterion_SWD==4) THEN            
                covIter_errSWD =  Ferro_norm(CdSWD-CdSWD_old,NDAT_SWD)
            END IF
            IF (rank==src) WRITE(*,*) 'SWD covariance convergence criterion value: ', covIter_errSWD
            IF (covIter_errSWD < converge_threshold_SWD) cov_converged_datasets(1)=1
        END IF
    ELSE !!ierror
        ICOV_SWD = 2
        ISD_SWD = ISD_SWD_covIter
        sdmn(1) = sdmn_covIter(1)
        sdmx(1) = sdmx_covIter(1)
        obj%sdparSWD = sdpar_covIter(1) 
        minlimsdSWD   = sdmn(1)
        maxlimsdSWD   = sdmx(1)
        pertsdsdscSWD = 10._RP
        maxpertsdSWD  = maxlimsdSWD-minlimsdSWD
        pertsdsdSWD   = maxpertsdSWD/pertsdsdscSWD
    END IF  !!ierror

END IF

IF (  cov_converged_datasets(2)==0 ) THEN
    
    IF (ICOV_ELL==2) THEN
        CdELL_old = CdELL
        CdiELL_old = CdiELL
    END IF
    is = NMODE2*NDAT_SWD + 1
    ie = NMODE2*NDAT_SWD + NMODE_ELL2*NDAT_ELL
    CALL COVmatEst(CdELL, CdiELL, sampleDres(:,is:ie), sampleDres(:,ncount3), NsampleDres, NDAT_ELL, nfrac_ELL, MAX_NAVE_ELL, damp_power_ELL, inonstat_ELL, iunbiased_ELL, imr_ELL, ierror)
    IF (ierror/=0) THEN
        WRITE(*,*) 'non positive definite ELL cov estimate at iter: ',icovIter
        WRITE(*,*) 'algorithm continues with previous cov matrix'
    IF (ICOV_ELL==2) THEN
        CdELL = CdELL_old
        CdiELL = CdiELL_old
    END IF
    ELSE IF (ICOV_ELL==2) THEN
        IF (iconverge_criterion<=2) THEN
            IF (iconverge_criterion_ELL==1) THEN
                covIter_errELL = MAXVAL(ABS(CdELL-CdELL_old))  / MAXVAL(ABS(CdELL_old))
            ELSEIF (iconverge_criterion_ELL==2) THEN            
                covIter_errELL =  MAXVAL(ABS(CdELL-CdELL_old)) 
            ELSEIF (iconverge_criterion_ELL==3) THEN            
                covIter_errELL =  Ferro_norm(CdELL-CdELL_old,NDAT_ELL) / Ferro_norm(CdELL_old,NDAT_ELL)
            ELSEIF (iconverge_criterion_ELL==4) THEN            
                covIter_errELL =  Ferro_norm(CdELL-CdELL_old,NDAT_ELL)
            END IF
            IF (rank==src) WRITE(*,*) 'ELL covariance convergence criterion value: ', covIter_errELL
            IF (covIter_errELL < converge_threshold_ELL) cov_converged_datasets(2)=1
        END IF
    ELSE !!ierror
        ICOV_ELL = 2
        ISD_ELL = ISD_ELL_covIter
        sdmn(2) = sdmn_covIter(2)
        sdmx(2) = sdmx_covIter(2)
        obj%sdparELL = sdpar_covIter(2)
        minlimsdELL   = sdmn(2)
        maxlimsdELL   = sdmx(2)
        pertsdsdscELL = 10._RP
        maxpertsdELL  = maxlimsdELL-minlimsdELL
        pertsdsdELL   = maxpertsdELL/pertsdsdscELL
    END IF !!ierror

END IF

!! check for cov convergence
IF (iconverge_criterion<=2) THEN
   
    IF ( ALL( cov_converged_datasets == 1) ) THEN
            cov_converged = .TRUE.
            IF(rank==src) CALL SAVECOVS()
    ELSE 
        IF ( iconverge_criterion==1 ) THEN
            WHERE ( ICOViter_datasets==1 ) cov_converged_datasets=0
        END IF
    END IF

!ELSEIF (iconverge_criterion==3) THEN
   
    
!    covIter_err = ( Ferro_norm(Cd-Cd_old) + Ferro_norm(CdSWD-CdSWD_old) + Ferro_norm(CdELL-CdELL_old) + ZFerro_norm(CdMT-CdMT_old) ) &
 !               / ( Ferro_norm(Cd_old) + Ferro_norm(CdSWD_old) + Ferro_norm(CdELL_old) + ZFerro_norm(CdMT_old ) )
    
END IF  

END SUBROUTINE UpdateCOV
!!----------------------------

FUNCTION Ferro_norm(A, N)

USE DATA_TYPE
IMPLICIT NONE

INTEGER(KIND=IB) :: N
REAL(KIND=RP) :: A(N,N)
REAL(KIND=RP) :: Ferro_norm

Ferro_norm = SQRT( SUM(A**2) / REAL(N,KIND=RP) )

RETURN
END FUNCTION Ferro_norm
!!-----------------------------


!!-----------------------------

SUBROUTINE SAVECOVS()

USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: irow                      


IF (ICOV_iterUpdate_SWD==1) THEN

    OPEN(UNIT=417, FILE=infileCdSWD, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NDAT_SWD
        WRITE(417,100) CdSWD(irow,:)
    END DO 
    CLOSE(417)
    OPEN(UNIT=418, FILE=infileCdiSWD, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NDAT_SWD
        WRITE(418,100) CdiSWD(irow,:)
    END DO 
    CLOSE(418)

END IF

IF (ICOV_iterUpdate_ELL==1) THEN

    OPEN(UNIT=517, FILE=infileCdELL, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NDAT_ELL
        WRITE(517,100) CdELL(irow,:)
    END DO 
    CLOSE(517)
    OPEN(UNIT=518, FILE=infileCdiELL, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NDAT_ELL
        WRITE(518,100) CdiELL(irow,:)
    END DO 
    CLOSE(518)

END IF

100 FORMAT(1000000ES20.10)

END SUBROUTINE SAVECOVS
