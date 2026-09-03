!==============================================================================

SUBROUTINE READPARFILE()
!==============================================================================
!USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=RP):: ntmp,iaz,ik,ilim
REAL(KIND=RP)   :: vref, dVs, VpVsmin, VpVsmax, dVpVs
CHARACTER(LEN=256) :: kwline
CHARACTER(LEN=32)  :: kw
INTEGER(KIND=IB)   :: io_kw,ipos,ic_kw,ntmp2
INTERFACE
   FUNCTION LOGFACTORIAL(n)
     USE DATA_TYPE
     REAL(KIND=RP) :: LOGFACTORIAL
     REAL(KIND=RP),INTENT(IN):: n
   END FUNCTION LOGFACTORIAL
END INTERFACE

parfile        = filebase(1:filebaselen) // '_parameter.dat'!! Inversion parameter file
!!
!! Read parameter file
!!
OPEN(UNIT=20,FILE=parfile,FORM='formatted',STATUS='OLD',ACTION='READ')
READ(20,*) IMAP       !!  1 predict data for the map_voro model and exit
READ(20,*) IMAGSCALE  !!  2 magnitude-scaled error model
READ(20,*) ENOS       !!  3 even-numbered order statistics prior on k
READ(20,*) IPOIPR     !!  4 Poisson prior on k
READ(20,*) IAR        !!  5 autoregressive error model
READ(20,*) I_VARPAR   !!  6 variable layer complexity (trans-D)
READ(20,*) IBD_SINGLE !!  7 birth/death for single parameters onto nodes
READ(20,*) I_SWD      !!  8 invert SWD data
READ(20,*) I_ELL      !!  9 invert ELL data
READ(20,*) I_VREF     !! 10 perturbations around <base>_vel_ref.txt
READ(20,*) I_VPVS     !! 11 sample Vp/Vs ratio
READ(20,*) ISMPPRIOR  !! 12 sample the prior
READ(20,*) ISETSEED   !! 13 fixed random seed
READ(20,*) IEXCHANGE  !! 14 parallel-tempering exchange moves
READ(20,*) NDAT_SWD   !! 15 max no. SWD data per curve
READ(20,*) NMODE      !! 16 no. SWD curves
READ(20,*) NDAT_ELL   !! 17 no. ELL data
READ(20,*) NMODE_ELL  !! 18 no. ELL modes
READ(20,*) NLMN       !! 19 min number of nodes
READ(20,*) NLMX       !! 20 max number of nodes
READ(20,*) ICHAINTHIN !! 21 chain thinning interval
READ(20,*) NKEEP      !! 22 no. models to keep before writing
READ(20,*) NPTCHAINS1 !! 23 no. chains at T=1
READ(20,*) dTlog      !! 24 temperature increment
READ(20,*) lambda     !! 25 lambda for Poisson prior on k
READ(20,*) hmx        !! 26 max. partition depth [km]
READ(20,*) hmin       !! 27 min. layer thickness [km]
READ(20,*) armxSWD    !! 28 max. AR prediction size SWD
READ(20,*) armxELL    !! 29 max. AR prediction size ELL
READ(20,*) TCHCKPT    !! 30 checkpointing interval in s
READ(20,*) dVs        !! 31 Vs one sided prior width (relative to background model)
READ(20,*) dVpVs      !! 32 VpVs ratio one sided prior width
READ(20,*) sdmn       !! 33 sigma prior lower bounds (SWD, ELL)
READ(20,*) sdmx       !! 34 sigma prior upper bounds (SWD, ELL)
READ(20,*) ISD_SWD    !! 35 1 = hierarchical sigma sampled (SWD)
READ(20,*) ISD_ELL    !! 36 1 = hierarchical sigma sampled (ELL)
READ(20,*) ICOV_SWD   !! 37 SWD likelihood (0 implicit, 1 hierarchical, 2 Cdi file, 3 sd file)
READ(20,*) ICOV_ELL   !! 38 ELL likelihood
READ(20,*) ELL_verbose!! 39 
READ(20,*) ELL_prec   !! 40
READ(20,*) I_ABS_ELL  !! 41 
READ(20,*) I_LOG10_ELL!! 42 
READ(20,*) I_SAMPLING_TYPE_ELL !! 43
READ(20,*) I_SET_STEP_ELL !! 44
READ(20,*) STEP_SIZE_ELL  !! 45
READ(20,*) I_SET_COUNT_ELL !! 46
READ(20,*) COUNT_ELL  !! 47 
READ(20,*) I_SET_RANGE_ELL !! 48
!READ(20,*) VpVsmin   !! minimum VpVs ratio
!READ(20,*) VpVsmax   !! maximum VpVs ratio
!!
!! ---- OPTIONAL KEYWORD LINES (multimode-raydsp branch) ----
!! Any line after the 48 positional ones of the form  KEYWORD value(s)  is
!! parsed; everything else (blank, comments, leftover legacy lines) is ignored,
!! so every pre-existing parameter file keeps working unchanged.
!!   DVSCON  x            max |adjacent-layer dVs| in km/s, indicator prior
!!                        checked before the forward (Kennett 2023/2026;
!!                        BayHunter lvz/hvz parity). Absent or < 0 = off.
!!   MODE_OF m1 m2 ...    Rayleigh mode number of each SWD curve slot (NMODE
!!                        ascending integers >= 0). Slot files are named by
!!                        mode: <base>_SWD.dat (mode 0), <base>_SWD_M<m>.dat.
!!                        Absent = 0 1 ... NMODE-1.  E.g. "MODE_OF 0 2" fits
!!                        the fundamental + second higher mode with no R1.
!!   IGRP    0|1          0 = phase velocity (default), 1 = group velocity.
!!   SWD_WARM 0|1|-1      warm-started root scan: march each mode's root from
!!                        the previous period instead of rescanning from cmin
!!                        (3-8x fewer propagator calls). -1 (default) enables
!!                        it only when DVSCON > 0, the regime where it is
!!                        exact; 1 forces it on, 0 off.
!!   SWD_SCAN cmin cmax dc [dc_over]  DISPER80 root-scan window and step in
!!                        km/s. Default 2.0 6.5 0.05 (crustal); near-surface
!!                        work needs e.g. 0.08 1.6 0.005. Overtones scan at
!!                        dc_over (default dc/5; give it explicitly, e.g.
!!                        0.001, for bit-reproducible runs across builds).
!!
ALLOCATE( MODE_OF(NMODE) )
DO ntmp2 = 1,NMODE
  MODE_OF(ntmp2) = ntmp2-1
ENDDO
DO
  READ(20,'(A)',IOSTAT=io_kw) kwline
  IF(io_kw /= 0) EXIT
  ipos = INDEX(kwline,'!')
  IF(ipos > 0) kwline = kwline(1:ipos-1)
  kwline = ADJUSTL(kwline)
  IF(LEN_TRIM(kwline) == 0) CYCLE
  ipos = INDEX(kwline,' ')
  IF(ipos <= 1) CYCLE
  kw = kwline(1:ipos-1)
  DO ic_kw = 1,LEN_TRIM(kw)
    IF(kw(ic_kw:ic_kw) >= 'a' .AND. kw(ic_kw:ic_kw) <= 'z') kw(ic_kw:ic_kw) = ACHAR(IACHAR(kw(ic_kw:ic_kw))-32)
  ENDDO
  SELECT CASE (TRIM(kw))
  CASE ('DVSCON')
    READ(kwline(ipos:),*,IOSTAT=io_kw) DVSCON
  CASE ('MODE_OF')
    READ(kwline(ipos:),*,IOSTAT=io_kw) MODE_OF
    IF(io_kw /= 0)THEN
      WRITE(6,*) 'ERROR: MODE_OF needs NMODE =',NMODE,' integers: ',TRIM(kwline)
      STOP
    ENDIF
  CASE ('SWD_WARM')
    READ(kwline(ipos:),*,IOSTAT=io_kw) SWD_WARM
  CASE ('IGRP')
    READ(kwline(ipos:),*,IOSTAT=io_kw) IGRP
  CASE ('SWD_SCAN')
    READ(kwline(ipos:),*,IOSTAT=io_kw) SWD_CMIN,SWD_CMAX,SWD_DC,SWD_DC_OVER
    IF(io_kw > 0)THEN
      WRITE(6,*) 'ERROR: SWD_SCAN needs cmin cmax dc [dc_overtones]: ',TRIM(kwline)
      STOP
    ENDIF
    IF(io_kw < 0)THEN   !! only three values given: overtone step defaults to dc/5
      READ(kwline(ipos:),*,IOSTAT=io_kw) SWD_CMIN,SWD_CMAX,SWD_DC
      SWD_DC_OVER = -1._RP
    ENDIF
  CASE DEFAULT
    !! not a keyword (e.g. legacy trailing lines) - ignored
  END SELECT
  io_kw = 0
ENDDO
DO ntmp2 = 2,NMODE
  IF(MODE_OF(ntmp2) <= MODE_OF(ntmp2-1) .OR. MODE_OF(1) < 0)THEN
    WRITE(6,*) 'ERROR: MODE_OF must be ascending and >= 0:',MODE_OF
    STOP
  ENDIF
ENDDO
CLOSE(20)

IF (I_SWD==1) THEN
    NMODE2 = NMODE
ELSE
    NMODE2 = 0_IB
ENDIF

IF (I_ELL==1) THEN
    NMODE_ELL2 = NMODE_ELL
ELSE
    NMODE_ELL2 = 0_IB
ENDIF



kmin = NLMN
kmax = NLMX
!kmin = 1
!kmax = 8

!! Poisson Prior on k:
ALLOCATE(pk(kmax))
DO ik = kmin,kmax
  pk(ik)  = EXP(-lambda)*lambda**REAL(ik,RP)/EXP(LOGFACTORIAL(REAL(ik,RP)))
ENDDO

infileSWD      = filebase(1:filebaselen) // '_SWD.dat'
infile_sdSWD   = filebase(1:filebaselen) // '_sdSWD.dat'
infileELL      = filebase(1:filebaselen) // '_ELL.dat'
infileCdiSWD   = filebase(1:filebaselen) // '_CdiSWD.dat'
infileCdiELL   = filebase(1:filebaselen) // '_CdiELL.dat'
infileref      = filebase(1:filebaselen) // '_vel_ref.txt'
logfile        = filebase(1:filebaselen) // '_RJMH.log'
seedfile       = filebase(1:filebaselen) // '_seeds.log'
mapfile        = filebase(1:filebaselen) // '_map_voro.dat'
arfileSWD      = filebase(1:filebaselen) // '_maparSWD.dat'
predfileSWD    = filebase(1:filebaselen) // '_mappredSWD.dat'
obsfileSWD     = filebase(1:filebaselen) // '_obsSWD.dat'
arfileELL      = filebase(1:filebaselen) // '_maparELL.dat'
predfileELL    = filebase(1:filebaselen) // '_mappredELL.dat'
obsfileELL     = filebase(1:filebaselen) // '_obsELL.dat'
!covfile        = filebase(1:filebaselen) // '_cov.txt'
sdfile         = filebase(1:filebaselen) // '_sigma.txt'
samplefile     = filebase(1:filebaselen) // '_voro_sample.txt'
stepsizefile     = filebase(1:filebaselen) // '_stepsize.txt'
!! Parameters per Voronoi node: depth, Vs, [Vp/Vs]
IF(I_VPVS == 1) THEN
  NPL = 3
ELSE
  NPL = 2
END IF
NFPMX  = NLMX * NPL
NFPMX2 = NLMX * NPL * NPL

ioutside = 0;ireject  = 0; iaccept = 0; iaccept_delay = 0; ireject_delay = 0
i_sdpert = 0;ishearfail = 0 ;i_ref_nlay = 0

!!
!!
!! Read velocity reference file
!!
201 FORMAT(a64)
202 FORMAT(a28,I4,I4)
203 FORMAT(4F12.3)
IF(I_VREF == 1)THEN
  OPEN(UNIT=20,FILE=infileref,FORM='formatted',STATUS='OLD',ACTION='READ')
  READ(20,*) NVELREF, ntmp
  NPREM = ntmp - NVELREF
  ALLOCATE(vel_ref(4,NVELREF),vel_prem(4,NPREM))
  IF(rank == src)WRITE(6,201) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
  IF(rank == src)WRITE(6,202) ' Velocity reference model:  ',NVELREF,NPREM
  IF(rank == src)WRITE(6,201) ' Reference:                                                     '
  IF(rank == src)WRITE(6,201) '   Depth(km)   Vs (km/s)      VpVs        Density               '
  !!  Read reference model to max sampling depth
  DO iaz = 1,NVELREF
    READ(20,*) vel_ref(:,iaz)
    IF(rank == src)WRITE(6,203) vel_ref(:,iaz)
  ENDDO
  !!  Read PREM beyond that
  IF(rank == src)WRITE(6,201) ' Prem (deep reference):                                         '
  DO iaz = 1,NPREM
    READ(20,*) vel_prem(:,iaz)
    IF(rank == src)WRITE(6,203) vel_prem(:,iaz)
  ENDDO
  CLOSE(20)
  IF(rank == src)WRITE(6,201) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
ENDIF

ALLOCATE(idxpar(NPL),sdevm((NLMX*NPL)+NPL-1,NLMX))
sdevm = 0._RP

ALLOCATE(minlim(NPL),maxlim(NPL),maxpert(NPL),pertsd(NPL),pertsdsc(NPL))
minlim   = 0._RP;maxlim   = 0._RP;maxpert  = 0._RP
pertsd   = 0._RP;pertsdsc = 30._RP

ALLOCATE(minlimarSWD(NMODE),maxlimarSWD(NMODE),maxpertarSWD(NMODE))
ALLOCATE(pertarsdSWD(NMODE),pertarsdscSWD(NMODE))
minlimarSWD  = 0._RP;maxlimarSWD  = 0._RP;maxpertarSWD = 0._RP
pertarsdSWD  = 0._RP;pertarsdscSWD= 18._RP    !! better to pertarSWDsd and pertarSWDsdsc, change them in the future

ALLOCATE(minlimarELL(NMODE_ELL),maxlimarELL(NMODE_ELL),maxpertarELL(NMODE_ELL))
ALLOCATE(pertarsdELL(NMODE_ELL),pertarsdscELL(NMODE_ELL))
minlimarELL  = 0._RP;maxlimarELL  = 0._RP;maxpertarELL = 0._RP
pertarsdELL  = 0._RP;pertarsdscELL= 18._RP

ALLOCATE(minlimsdSWD(NMODE),maxlimsdSWD(NMODE),maxpertsdSWD(NMODE))
ALLOCATE(pertsdsdSWD(NMODE),pertsdsdscSWD(NMODE))
minlimsdSWD  = 0._RP;maxlimsdSWD  = 0._RP;maxpertsdSWD = 0._RP
pertsdsdSWD  = 0._RP;pertsdsdscSWD= 18._RP    !! better to pertsdSWDsd and pertsdSWDsdsc. change them in the future

ALLOCATE(minlimsdELL(NMODE_ELL),maxlimsdELL(NMODE_ELL),maxpertsdELL(NMODE_ELL))
ALLOCATE(pertsdsdELL(NMODE_ELL),pertsdsdscELL(NMODE_ELL))
minlimsdELL  = 0._RP;maxlimsdELL  = 0._RP;maxpertsdELL = 0._RP
pertsdsdELL  = 0._RP;pertsdsdscELL= 18._RP

!! curmod columns: 1 = thickness, 3 = Vp (or Vp/Vs before scaling), 4 = Vs
IF(I_VPVS == 1)THEN
  idxpar   = (/ 1, 4, 3 /)        !! Sample Vs and Vp/Vs
ELSEIF(I_VPVS == -1)THEN
  idxpar(1:2)   = (/ 1, 4 /)      !! Sample Vs (Vp = 1.75 Vs)
ELSE
  idxpar   = (/ 1, 3, 4 /)        !! Sample Vp
ENDIF
!!
!!  Prior bounds
!! (Note: Density is empirical through Birch's Law)
!!
!! Without dip:
!!           h     vs     vp/vs
ilim = 1
minlim(ilim) = hmin
maxlim(ilim) = hmx
ilim = ilim + 1
minlim(ilim) = -dVs
maxlim(ilim) = dVs
ilim = ilim + 1
IF(I_VPVS == 1)THEN
  !! Sample Vs and VpVs ratio
  minlim(ilim) =  -dVpVs
  maxlim(ilim) =  dVpVs
  ilim = ilim + 1
ENDIF

maxpert = maxlim-minlim
pertsd = maxpert/pertsdsc

IF(IAR == 1)THEN
  !! Set prior and proposal scaling for AR model:
  minlimarSWD   = -0.5000_RP
  maxlimarSWD   =  0.90_RP
  pertarsdscSWD =  10._RP
  maxpertarSWD  = maxlimarSWD-minlimarSWD
  pertarsdSWD   = maxpertarSWD/pertarsdscSWD
  minlimarELL   = -0.5000_RP
  maxlimarELL   =  0.90_RP
  pertarsdscELL =  10._RP
  maxpertarELL  = maxlimarELL-minlimarELL
  pertarsdELL   = maxpertarELL/pertarsdscELL
ENDIF

IF(ICOV_SWD >= 1)THEN
  !! Set prior and proposal scaling for data error standard deviations:
  minlimsdSWD   = sdmn(1)
  maxlimsdSWD   = sdmx(1)
  pertsdsdscSWD = 10._RP  !! also it is set in UpdateCOV()
  maxpertsdSWD  = maxlimsdSWD-minlimsdSWD
  pertsdsdSWD   = maxpertsdSWD/pertsdsdscSWD
ENDIF
IF(ICOV_ELL >= 1)THEN
  minlimsdELL   = sdmn(2)
  maxlimsdELL   = sdmx(2)
  pertsdsdscELL = 10._RP  !! also it is set in UpdateCOV()
  maxpertsdELL  = maxlimsdELL-minlimsdELL
  pertsdsdELL   = maxpertsdELL/pertsdsdscELL
ENDIF

END SUBROUTINE READPARFILE
!==============================================================================
!!
!! Write some info:
!!
!!=============================================================================

SUBROUTINE PRINTPAR2()

USE RJMCMC_COM
IMPLICIT NONE

!IF(rank == src)THEN
  WRITE(6,*) 'IMAP      = ', IMAP
  WRITE(6,*) 'ICOV_SWD  = ', ICOV_SWD
  WRITE(6,*) 'ICOV_ELL  = ', ICOV_ELL
  WRITE(6,*) 'I_SWD     = ', I_SWD
  WRITE(6,*) 'I_ELL     = ', I_ELL
  WRITE(6,*) 'IAR       = ', IAR
  WRITE(6,*) 'I_VPVS    = ', I_VPVS
  WRITE(6,*) 'ISMPPRIOR = ', ISMPPRIOR
  WRITE(6,*) 'ISETSEED  = ', ISETSEED
  WRITE(6,*) 'IEXCHANGE = ', IEXCHANGE
  WRITE(6,*) 'NLMN      = ', NLMN       !! Max number of layers
  WRITE(6,*) 'NLMX      = ', NLMX       !! Max number of layers
  WRITE(6,*) 'ICHAINTHIN= ', ICHAINTHIN !! Chain thinning interval
  WRITE(6,*) 'NKEEP     = ', NKEEP      !! No. models to keep before writing
  IF (icovIter>0_IB) WRITE(6,*) 'NKEEP2    = ', NKEEP2      !! No. models to keep before writing
  WRITE(6,*) 'NPTCHAINS1= ', NPTCHAINS1 !! Chain thinning interval
  WRITE(6,*) 'dTlog     = ', dTlog      !! Temperature increment
  WRITE(6,*) 'hmx       = ', hmx        !! Max. partition depth
  WRITE(6,*) 'hmin      = ', hmin       !! Min. layer thickness (must be small enough to not violate detailed balance)
  WRITE(6,*) 'TCHCKPT   = ', TCHCKPT    !! Checkpointing interval in s
  WRITE(6,*) 'ISD_SWD    = ', ISD_SWD                      
  WRITE(6,*) 'ISD_ELL    = ', ISD_ELL                      
  WRITE(6,*) 'Sample file: ',samplefile
  IF(ICOV_iterUpdate==1) WRITE(6,*) 'Sample file covIter   : ', samplefile_covIter
  WRITE(6,*) ''
  WRITE(6,*) 'minlim:  '
  WRITE(6,203) minlim
  WRITE(6,*) 'maxlim:  '
  WRITE(6,203) maxlim
  WRITE(6,*) 'pertsdsc:'
  WRITE(6,203) pertsdsc
  WRITE(6,*) 'minlim sigma(std) SWD:  '
  WRITE(6,203) minlimsdSWD
  WRITE(6,*) 'maxlim sigma(std) SWD:  '
  WRITE(6,203) maxlimsdSWD
  WRITE(6,*) 'pertsdsdscSWD:'
  WRITE(6,203) pertsdsdscSWD
  WRITE(6,*) 'minlim sigma(std) ELL:  '
  WRITE(6,203) minlimsdELL
  WRITE(6,*) 'maxlim sigma(std) ELL:  '
  WRITE(6,203) maxlimsdELL
  WRITE(6,*) 'pertsdsdscELL:'
  WRITE(6,203) pertsdsdscELL
  !WRITE(6,*) 'Done reading data.'
  WRITE(6,*) '--- multimode SWD keywords ---'
  WRITE(6,*) 'DVSCON     = ', DVSCON
  WRITE(6,*) 'MODE_OF    = ', MODE_OF
  WRITE(6,*) 'IGRP       = ', IGRP
  WRITE(6,*) 'SWD_SCAN   = ', SWD_CMIN, SWD_CMAX, SWD_DC, SWD_DC_OVER
  WRITE(6,*) 'SWD_WARM   = ', SWD_WARM
  IF (icovIter==0_IB) WRITE(6,*) 'Done reading parameter file.'
  WRITE(6,*) ''
  WRITE(6,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
  WRITE(6,*) ''
!ENDIF
CALL FLUSH(6)

203 FORMAT(100ES12.3)

RETURN
END SUBROUTINE PRINTPAR2   
!==============================================================================

SUBROUTINE READDATA(obj)
!==============================================================================
!USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=RP):: iaz,idat,io
INTEGER(KIND=IB):: imode,nswd_m
CHARACTER(LEN=100) :: swdfile
CHARACTER(LEN=8)   :: modestr
TYPE(objstruc)  :: obj

IF(I_SWD == 1)THEN
  !!
  !! Surface wave dispersion data:
  !!
  !!
  !! One file per curve slot, named by its Rayleigh MODE number (MODE_OF):
  !!   mode 0 -> <base>_SWD.dat, mode m -> <base>_SWD_M<m>.dat
  !! Each is read to EOF; the count goes into NDAT_MODE, so curves may have
  !! different numbers of points on different frequency grids (NDAT_SWD is
  !! only the array width).
  !!
  IF(.NOT.ALLOCATED(NDAT_MODE)) ALLOCATE( NDAT_MODE(NMODE) )
  NDAT_MODE = 0
  obj%periods = 0._RP
  obj%DobsSWD = 0._RP
  DO imode = 1,NMODE
    IF(MODE_OF(imode) == 0)THEN
      swdfile = infileSWD
    ELSE
      WRITE(modestr,'(I0)') MODE_OF(imode)
      swdfile = filebase(1:filebaselen) // '_SWD_M' // TRIM(modestr) // '.dat'
    ENDIF
    OPEN(20,FILE=swdfile,FORM='formatted',STATUS='OLD',ACTION='READ',IOSTAT=io)
    IF(io /= 0)THEN
      WRITE(6,*) 'ERROR: cannot open SWD data file for mode ',MODE_OF(imode),': ',TRIM(swdfile)
      STOP
    ENDIF
    nswd_m = 0
    DO idat=1,NDAT_SWD
       READ(20,*,IOSTAT=io) obj%periods(imode,idat),obj%DobsSWD(imode,idat)
       IF (io > 0) THEN
         STOP "Check input.  Something was wrong"
       ELSEIF (io < 0) THEN
         EXIT
       ELSE
         nswd_m = nswd_m + 1
       ENDIF
    ENDDO
    CLOSE(20)
    IF(nswd_m == 0)THEN
      WRITE(6,*) 'ERROR: no data read for SWD mode ',MODE_OF(imode),' from ',TRIM(swdfile)
      STOP
    ENDIF
    NDAT_MODE(imode) = nswd_m
    WRITE(6,*) 'SWD mode ',MODE_OF(imode),':',nswd_m,' points from ',TRIM(swdfile)
    IF(nswd_m == NDAT_SWD)THEN
      WRITE(6,*) '  NOTE: hit NDAT_SWD; increase it if this file has more rows.'
    ENDIF
  ENDDO
  IF(ICOV_SWD == 2)THEN
    !! SWD case, read CdSWD^-1 from file:
    OPEN(UNIT=30,FILE=infileCdiSWD,FORM='formatted',STATUS='OLD',ACTION='READ')
    !!
    !! Cov is NDAT_SWD by NDAT_SWD
    !!
    !!ALLOCATE(CdiSWD(NDAT_SWD,NDAT_SWD))
    CdiSWD = 0._RP
    DO iaz = 1,NDAT_SWD
      READ(30,*) CdiSWD(iaz,1:NDAT_SWD)
    ENDDO
    CLOSE(30)
  ELSEIF (ICOV_SWD == 3) THEN
    ALLOCATE(sdSWD(NMODE,NDAT_SWD))
    sdSWD = 0._RP
    OPEN(20,FILE=infile_sdSWD,FORM='formatted',STATUS='OLD',ACTION='READ')
    DO idat=1,NDAT_SWD
       READ(20,*,IOSTAT=io) sdSWD(1,idat)
       IF (io > 0) THEN
         STOP "Check input.  Something was wrong"
       ELSEIF (io < 0) THEN
         EXIT
       ELSE
       !  ndatad=ndatad+1
       ENDIF
       !if (i==ndatadmax) stop "number of Dispersion data >= ndatadmax"
    ENDDO
    CLOSE(20)
  ENDIF!!ICOV_SWD
ENDIF
IF(I_ELL == 1)THEN
  !!
  !! Surface wave dispersion data:
  !!
  OPEN(20,FILE=infileELL,FORM='formatted',STATUS='OLD',ACTION='READ')
  !ndat = 0
  DO idat=1,NDAT_ELL
     READ(20,*,IOSTAT=io) obj%periods_ELL(1,idat),obj%DobsELL(1,idat)
     IF (io > 0) THEN
       STOP "Check input.  Something was wrong"
     ELSEIF (io < 0) THEN
       EXIT
     ELSE
     !  ndatad=ndatad+1
     ENDIF
     !if (i==ndatadmax) stop "number of Dispersion data >= ndatadmax"
  ENDDO
  CLOSE(20)! close the file
  IF(ICOV_ELL == 2)THEN
    !! ELL case, read CdELL^-1 from file:
    OPEN(UNIT=30,FILE=infileCdiELL,FORM='formatted',STATUS='OLD',ACTION='READ')
    !!
    !! Cov is NDAT_ELL by NDAT_ELL
    !!
    !!ALLOCATE(CdiELL(NDAT_ELL,NDAT_ELL))
    CdiELL = 0._RP
    DO iaz = 1,NDAT_ELL
      READ(30,*) CdiELL(iaz,1:NDAT_ELL)
    ENDDO
    CLOSE(30)
  ENDIF
ENDIF
RETURN
END SUBROUTINE READDATA
!==============================================================================

SUBROUTINE PRINTPAR(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER :: ivo
TYPE(objstruc) :: obj

WRITE(6,*) 'Voronoi nodes:',obj%k
WRITE(6,*) 'NFP:',obj%NFP

DO ivo=1,obj%k
   WRITE(6,201) obj%voro(ivo,1:NPL)
   WRITE(6,205) obj%voroidx(ivo,1:NPL)
ENDDO
WRITE(6,*) 'Layer parameter vector:',obj%nunique,'layers'
DO ivo=1,obj%nunique
   WRITE(6,201) obj%par((ivo-1)*NPL+1:ivo*NPL)
ENDDO
WRITE(6,202) '            ',obj%par(obj%nunique*NPL+1:obj%nunique*NPL+(NPL-1))
WRITE(6,*) 'Partition:'
WRITE(6,203) obj%ziface(1:obj%nunique)
WRITE(6,*) 'Layers:'
WRITE(6,203) obj%hiface(1:obj%nunique)
IF(ICOV_SWD >= 1) THEN
  IF(I_SWD == 1) WRITE(6,206) 'sigma SWD = ',obj%sdparSWD
END IF
IF(ICOV_ELL >= 1) THEN
  IF(I_ELL == 1) WRITE(6,206) 'sigma ELL = ',obj%sdparELL
END IF
IF(IAR == 1)THEN
   WRITE(6,*) 'AR parameters:'
   WRITE(6,206) 'SWD:',obj%arparSWD
   WRITE(6,206) 'ELL:',obj%arparELL
ENDIF

201 FORMAT(6F12.4)
205 FORMAT(6I12)
202 FORMAT(A12,8F12.4)
203 FORMAT(11F12.4)
204 FORMAT(4F16.4)
!206 FORMAT(a,128F12.4)
206 FORMAT(a,128ES12.3)
END SUBROUTINE PRINTPAR
!!=======================================================================
RECURSIVE FUNCTION LOGFACTORIAL(n)  RESULT(fact)
!-----Factorial------------------------------------------------------
!!=======================================================================

USE DATA_TYPE
IMPLICIT NONE
REAL(KIND=RP) :: fact
REAL(KIND=RP), INTENT(IN) :: n

IF (n == 0) THEN
   fact = 0
ELSE
   fact = LOG(n) + LOGFACTORIAL(n-1)
END IF

END FUNCTION LOGFACTORIAL
!==============================================================================

SUBROUTINE READCOVPARFILE()
!!=============================================================================

USE RJMCMC_COM
IMPLICIT NONE


covparfile = filebase(1:filebaselen) // '_covparameter.dat'!! Inversion iterative cov parameter file
!!
!! Read parameter file
!!
OPEN(UNIT=20,FILE=covparfile,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
READ(20,*) Icov_iterUpdate_SWD
READ(20,*) Icov_iterUpdate_ELL
READ(20,*) covIter_zero_nsamples
READ(20,*) covIter_period 
READ(20,*) MAXcovIter
READ(20,*) ICOVest    
READ(20,*) CHAINTHIN_COVest_period_zeroIter
READ(20,*) CHAINTHIN_COVest_period_nonzeroIter
READ(20,*) ISD_SWD_covIter
READ(20,*) ISD_ELL_covIter
READ(20,*) sdmn_covIter 
READ(20,*) sdmx_covIter 
READ(20,*) sdpar_covIter 
READ(20,*) NKEEP_covIter 
READ(20,*) NKEEP_covIter_res 
READ(20,*) iSAVEsample_covIter 
READ(20,*) iSAVEsample_only_zeroIter 
READ(20,*) iMAP_calc 
READ(20,*) iconverge_criterion 
READ(20,*) iconverge_criterion_SWD 
READ(20,*) iconverge_criterion_ELL 
READ(20,*) converge_threshold_SWD 
READ(20,*) converge_threshold_ELL 
READ(20,*) nfrac_SWD
READ(20,*) MAX_NAVE_SWD
READ(20,*) inonstat_SWD
READ(20,*) iunbiased_SWD
READ(20,*) imr_SWD
READ(20,*) damp_power_SWD
READ(20,*) nfrac_ELL
READ(20,*) MAX_NAVE_ELL
READ(20,*) inonstat_ELL
READ(20,*) iunbiased_ELL
READ(20,*) imr_ELL
READ(20,*) damp_power_ELL
CLOSE(20) 

samplefile_covIter         = filebase(1:filebaselen) // '_voro_sample_covIter.txt'
samplefile_res_covIter     = filebase(1:filebaselen) // '_voro_sample_res_covIter.txt'
infileCdSWD                = filebase(1:filebaselen) // '_CdSWD.dat'
infileCdELL                = filebase(1:filebaselen) // '_CdELL.dat'

!!
!! Write some info:
!!
IF(rank == src)THEN
  WRITE(6,*) 'Icov_iterUpdate_SWD           = ', Icov_iterUpdate_SWD
  WRITE(6,*) 'Icov_iterUpdate_ELL           = ', Icov_iterUpdate_ELL
  WRITE(6,*) 'covIter_zero_nsamples         = ', covIter_zero_nsamples
  WRITE(6,*) 'covIter_period                = ', covIter_period
  WRITE(6,*) 'MAXcovIter                    = ', MAXcovIter 
  WRITE(6,*) 'ICOVest                       = ', ICOVest
  WRITE(6,*) 'CHAINTHIN_COVest_period_zeroIter     = ', CHAINTHIN_COVest_period_zeroIter
  WRITE(6,*) 'CHAINTHIN_COVest_period_nonzeroIter  = ', CHAINTHIN_COVest_period_nonzeroIter
  WRITE(6,*) 'NKEEP_covIter                 = ', NKEEP_covIter
  WRITE(6,*) 'NKEEP_covIter_res             = ', NKEEP_covIter_res
  WRITE(6,*) 'iSAVEsample_covIter           = ', iSAVEsample_covIter
  WRITE(6,*) 'iSAVEsample_only_zeroIter     = ', iSAVEsample_only_zeroIter
  WRITE(6,*) 'iMAP_calc                     = ', iMAP_calc 
  WRITE(6,*) 'iconverge_criterion           = ', iconverge_criterion
  WRITE(6,*) 'iconverge_criterion_SWD       = ', iconverge_criterion_SWD
  WRITE(6,*) 'iconverge_criterion_ELL       = ', iconverge_criterion_ELL
  WRITE(6,202) 'converge_threshold_SWD         = ', converge_threshold_SWD
  WRITE(6,202) 'converge_threshold_ELL         = ', converge_threshold_ELL
  WRITE(6,*) 'nfrac_SWD                     = ', nfrac_SWD
  WRITE(6,*) 'MAX_NAVE_SWD                  = ', MAX_NAVE_SWD
  WRITE(6,*) 'inonstat_SWD                  = ', inonstat_SWD
  WRITE(6,*) 'iunbiased_SWD                 = ', iunbiased_SWD
  WRITE(6,*) 'imr_SWD                       = ', imr_SWD  
  WRITE(6,202) 'damp_power_SWD                 = ', damp_power_SWD   
  WRITE(6,*) 'nfrac_ELL                     = ', nfrac_ELL
  WRITE(6,*) 'MAX_NAVE_ELL                  = ', MAX_NAVE_ELL
  WRITE(6,*) 'inonstat_ELL                  = ', inonstat_ELL
  WRITE(6,*) 'iunbiased_ELL                 = ', iunbiased_ELL
  WRITE(6,*) 'imr_ELL                       = ', imr_ELL  
  WRITE(6,202) 'damp_power_ELL                 = ', damp_power_ELL   
  WRITE(6,*) 'ISD_SWD_covIter    = ', ISD_SWD_covIter
  WRITE(6,*) 'ISD_ELL_covIter    = ', ISD_ELL_covIter
  WRITE(6,*) ''
  WRITE(6,*) 'min scaling factor:  '
  WRITE(6,203) sdmn_covIter
  WRITE(6,*) 'max scaling factor:  '
  WRITE(6,203) sdmx_covIter
  WRITE(6,*) 'starting scaling factor:  '
  WRITE(6,203) sdpar_covIter

  WRITE(6,*) 'Done reading covariance parameter file.'
  WRITE(6,*) ''
  WRITE(6,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
  WRITE(6,*) ''

END IF

202 FORMAT(a,100ES12.3)
203 FORMAT(100ES12.3)

END SUBROUTINE READCOVPARFILE
!!============================================================================

!EOF
