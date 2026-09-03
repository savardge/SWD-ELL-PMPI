!==============================================================================
MODULE RJMCMC_COM
!! Global module of the SWD-ELL trans-dimensional (rjMcMC) inversion.
!! RF and MT machinery removed (SWD-ELL-PMPI); see git history for the
!! original RF-SWD-ELL-MT code.
   USE MPI
   USE DATA_TYPE
   IMPLICIT NONE
!!
!! General switches
!!
   INTEGER(KIND=IB) :: IMAP       !! 1 = predict data for the map_voro model and exit
   INTEGER(KIND=IB) :: ICOV_SWD   !! SWD likelihood: 0 = implicit sigma, 1 = hierarchical sigma per curve, 2 = Cdi file, 3 = sd file
   INTEGER(KIND=IB) :: ICOV_ELL   !! ELL likelihood: same coding
   INTEGER(KIND=IB) :: IMAGSCALE  !! 1 = magnitude-scaled error model (SWD/ELL)
   INTEGER(KIND=IB) :: ENOS       !! 1 = even-numbered order statistics prior on node depths
   INTEGER(KIND=IB) :: IPOIPR     !! 1 = Poisson prior on k
   INTEGER(KIND=IB) :: IAR        !! 1 = autoregressive error model
   INTEGER(KIND=IB) :: I_VARPAR   !! 1 = variable layer complexity (trans-D)
   INTEGER(KIND=IB) :: IBD_SINGLE !! 1 = birth/death for single parameters onto nodes
   INTEGER(KIND=IB) :: I_SWD      !! 1 = invert SWD data
   INTEGER(KIND=IB) :: I_ELL      !! 1 = invert ELL data
   INTEGER(KIND=IB) :: I_VREF     !! 1 = node velocities are perturbations around <base>_vel_ref.txt
   INTEGER(KIND=IB) :: I_VPVS     !! 1 = sample Vp/Vs ratio (instead of Vp); -1 = fixed Vp/Vs = 1.75
   INTEGER(KIND=IB) :: ISMPPRIOR  !! 1 = sample the prior (logL = const)
   INTEGER(KIND=IB) :: ISETSEED   !! 1 = fixed random seed table
   INTEGER(KIND=IB) :: IEXCHANGE  !! 1 = parallel-tempering exchange moves
   INTEGER(KIND=IB),PARAMETER :: IDIP = 0     !! dipping layers were an RF-only feature (removed)
   INTEGER(KIND=IB) :: ISD_SWD    !! 1 = sample hierarchical sigma of SWD curves
   INTEGER(KIND=IB) :: ISD_ELL    !! 1 = sample hierarchical sigma of ELL curves
   INTEGER(KIND=IB) :: I_ABS_ELL  !! 1 = inverts abs of ellipticities
   INTEGER(KIND=IB) :: I_LOG10_ELL!! 1 = inverts log10 of absolute ellipticities
!!
!! Model and data dimensions
!!
   INTEGER(KIND=IB)            :: NDAT_SWD     ! Max number of SWD data per curve (array width)
   INTEGER(KIND=IB)            :: NMODE        ! Number of SWD curves inverted jointly
   INTEGER(KIND=IB)            :: NDAT_ELL     ! Number ellipticity data
   INTEGER(KIND=IB)            :: NMODE_ELL    ! Number ellipticity modes
   !! ---- multimode SWD (multimode-raydsp branch) ----
   INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:) :: NDAT_MODE  !! actual no. of data per SWD curve slot (NDAT_SWD = max)
   INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:) :: MODE_OF    !! Rayleigh mode number of each curve slot (keyword MODE_OF; default 0..NMODE-1)
   REAL(KIND=RP)    :: DVSCON   = -1._RP   !! max |adjacent-layer dVs| [km/s] indicator prior (keyword DVSCON; < 0 = off)
   INTEGER(KIND=IB) :: IGRP     = 0        !! 1 = group velocity, 0 = phase velocity (keyword IGRP)
   REAL(KIND=RP)    :: SWD_CMIN = 2.0_RP   !! DISPER80 phase-speed scan window [km/s] (keyword SWD_SCAN cmin cmax dc [dc_over])
   REAL(KIND=RP)    :: SWD_CMAX = 6.5_RP   !!   defaults = the original crustal values of dispersion.f90
   REAL(KIND=RP)    :: SWD_DC   = 0.05_RP  !!   scan step [km/s] for the fundamental
   REAL(KIND=RP)    :: SWD_DC_OVER = -1._RP !!   scan step for overtones (< 0 = dc/5)
   INTEGER(KIND=IB) :: SWD_WARM = -1     !! warm-started root scan (keyword SWD_WARM):
                                          !! 1 = on, 0 = off, -1 = auto (on iff DVSCON > 0)
   INTEGER(KIND=IB)            :: NLMN         ! Min number of layers
   INTEGER(KIND=IB)            :: NLMX         ! Max number of layers
   INTEGER(KIND=IB)            :: NPL          ! No. parameters per layer
   CHARACTER(len=64) :: filebasefile      = 'filebase.txt'
!!
!! Forward model parameter indexing (curmod columns: thick rho alpha beta ...)
!!
  INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxpar
  REAL(KIND=RP)             :: hmx                !! Max partition depth in km
  REAL(KIND=RP),DIMENSION(2):: sdmn               !! Min standard deviation (SWD, ELL)
  REAL(KIND=RP),DIMENSION(2):: sdmx               !! Max standard deviation (SWD, ELL)
!! 
!! Forward gpell specific
!!
  INTEGER(KIND=IB)      :: ELL_verbose
  REAL(KIND=RP)         :: ELL_prec
  INTEGER(KIND=IB)      :: I_SAMPLING_TYPE_ELL
  INTEGER(KIND=IB)      :: I_SET_STEP_ELL
  REAL(KIND=RP)         :: STEP_SIZE_ELL
  INTEGER(KIND=IB)      :: I_SET_COUNT_ELL
  INTEGER(KIND=IB)      :: COUNT_ELL  
  INTEGER(KIND=IB)      :: I_SET_RANGE_ELL
!!
!!  Prior variables and good seeding model
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlim
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxlim
   INTEGER(KIND=IB)            :: kmin     = 0       ! Min number of layers
   INTEGER(KIND=IB)            :: kmax     = 0       ! Max number of layers
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pk       ! Poisson prior on k
   REAL(KIND=RP)               :: lambda             ! Lambda parameter for Poisson prior on k
   REAL(KIND=RP)               :: hmin               ! Min allowed layer thickness
   REAL(KIND=RP),PARAMETER     :: fact     = 1.00_RP ! factor for rotated space perturbation
   REAL(KIND=RP),PARAMETER     :: factdelay= 1.50_RP ! shrinking factor for delayed rejection (>1.)
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxpert
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsd
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsc
!!
!!  Autoregressive model prior variables:
!!
   REAL(KIND=RP)              :: armxSWD    = 0.5_RP           ! Max AR model range (SWD)
   REAL(KIND=RP)              :: armxELL    = 0.5_RP           ! Max AR model range (ELL)
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimarSWD, maxlimarSWD, maxpertarSWD
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertarsdSWD, pertarsdscSWD
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimarELL, maxlimarELL, maxpertarELL
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertarsdELL, pertarsdscELL
!!
!!  Standard deviation prior variables:
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimsdSWD, maxlimsdSWD, maxpertsdSWD
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsdSWD, pertsdsdscSWD
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimsdELL, maxlimsdELL, maxpertsdELL
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsdELL, pertsdsdscELL
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: sdSWD 
   CHARACTER(len=64) :: filebase
   INTEGER(KIND=IB)  :: filebaselen
   CHARACTER(LEN=100) :: infileSWD
   CHARACTER(LEN=100) :: infile_sdSWD
   CHARACTER(LEN=100) :: infileELL
   CHARACTER(LEN=100) :: infileref
   CHARACTER(LEN=100) :: infileCdiSWD
   CHARACTER(LEN=100) :: infileCdiELL
   CHARACTER(LEN=100) :: parfile
   CHARACTER(LEN=64) :: logfile
   CHARACTER(LEN=64) :: seedfile
   CHARACTER(len=64) :: mapfile
   CHARACTER(LEN=64)  :: obsfileSWD
   CHARACTER(LEN=64)  :: arfileSWD
   CHARACTER(LEN=64)  :: predfileSWD
   CHARACTER(LEN=64)  :: obsfileELL
   CHARACTER(LEN=64)  :: arfileELL
   CHARACTER(LEN=64)  :: predfileELL
   CHARACTER(len=64) :: sdfile
   CHARACTER(len=64) :: samplefile
   CHARACTER(len=64) :: stepsizefile
!!
!! Velocity reference model
!!
  INTEGER(KIND=IB)                        :: NVELREF, NPREM
  REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: vel_ref, vel_prem
!!
!! Iterative covariance estimate parameters (datasets: 1 = SWD, 2 = ELL)
!!
   CHARACTER(len=64)                            :: samplefile_covIter
   CHARACTER(len=64)                            :: samplefile_res_covIter
   CHARACTER(LEN=100)                           :: covparfile
   CHARACTER(LEN=100)                           :: infileCdSWD
   CHARACTER(LEN=100)                           :: infileCdELL
   INTEGER(KIND=IB)                             :: usample_covIter, usample_res_covIter, units, unit2
   INTEGER(KIND=IB)                             :: ICOVest                          
   INTEGER(KIND=IB)                             :: icovIter
   INTEGER(KIND=IB)                             :: covIter_zero_nsamples, covITER_period, MAXcovIter
   INTEGER(KIND=IB)                             :: ICOV_iterUpdate, ICOV_iterUpdate_SWD, ICOV_iterUpdate_ELL
   LOGICAL                                      :: cov_converged
   INTEGER(KIND=IB)                             :: ISD_SWD_covIter, ISD_ELL_covIter
   INTEGER(KIND=IB)                             :: CHAINTHIN_COVest_period_zeroIter, CHAINTHIN_COVest_period_nonzeroIter
   INTEGER(KIND=IB)                             :: NMODE2, NMODE_ELL2, ncount3                                  
   REAL(KIND=RP),DIMENSION(2)                   :: sdmn_covIter               
   REAL(KIND=RP),DIMENSION(2)                   :: sdmx_covIter              
   REAL(KIND=RP),DIMENSION(2)                   :: sdpar_covIter      
   INTEGER(KIND=IB)                             :: NKEEP2, NKEEP_covIter, NKEEP3, NKEEP_covIter_res        
   INTEGER(KIND=IB)                             :: iSAVEsample_covIter, iSAVEsample_only_zeroIter
   INTEGER(KIND=IB)                             :: iMAP_calc           
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE     :: CdSWD, CdSWD_old, CdiSWD_old      
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE     :: CdELL, CdELL_old, CdiELL_old      
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE     :: sampleDres, sample2      
   INTEGER(KIND=IB)                             :: NsampleDres
   INTEGER(KIND=IB)                             :: iconverge_criterion, iconverge_criterion_SWD, iconverge_criterion_ELL
   REAL(KIND=RP)                                :: converge_threshold_SWD, converge_threshold_ELL
   REAL(KIND=RP)                                :: covIter_errSWD, covIter_errELL
   INTEGER(KIND=IB), DIMENSION(2)               :: cov_converged_datasets
   INTEGER(KIND=IB), DIMENSION(2)               :: ICOViter_datasets     
   INTEGER(KIND=IB)                             :: nfrac_SWD, MAX_NAVE_SWD, inonstat_SWD, iunbiased_SWD, imr_SWD
   REAL(KIND=RP)                                :: damp_power_SWD
   INTEGER(KIND=IB)                             :: nfrac_ELL, MAX_NAVE_ELL, inonstat_ELL, iunbiased_ELL, imr_ELL
   REAL(KIND=RP)                                :: damp_power_ELL
!!
!! Parallel Tempering parameters
!!
  INTEGER(KIND=IB)                         :: NPTCHAINS1            !! # chains T=1
  REAL(KIND=RP)                            :: dTlog                 ! Temperature increment
  INTEGER(KIND=IB)                         :: NT                    !! # tempering levels (temperatures)
  INTEGER(KIND=IB)                         :: NPTCHAINS             !! # parallel tempering chains
  INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: NCHAINT
  INTEGER(KIND=IB)                         :: ncswap     = 0_IB     !! Temp swap accepted counter
  INTEGER(KIND=IB)                         :: ncswapprop = 0_IB     !! Temp swap proposed counter
  REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)   :: beta_pt               !! Temperature array parallel tempering
  INTEGER(KIND=IB)                         :: ibirth  = 0,ideath  = 0
  INTEGER(KIND=IB)                         :: ibirths = 0,ideaths = 0
!!
!!  Sampling specific parameters
!!
   INTEGER(KIND=IB)           :: NFPMX
   INTEGER(KIND=IB)           :: NFPMX2
   INTEGER(KIND=IB)           :: ioutside   = 0
   INTEGER(KIND=IB)           :: ireject    = 0, iaccept = 0, iaccept_delay = 0, ireject_delay = 0
   INTEGER(KIND=IB)           :: i_bd,i_bds     ! Birth-Death track (0=MCMC, 1=birth, 2=death)
   INTEGER(KIND=IB)           :: i_sdpert = 0   ! if sigma is perturbed, don't compute forward model
   INTEGER(KIND=IB)           :: ishearfail = 0
   INTEGER(KIND=IB)           :: i_ref_nlay = 0
!!
!!  Convergence parameters
!!
   INTEGER(KIND=IB)       :: iconv    = 0       ! Convergence switch slaves
   INTEGER(KIND=IB)       :: iconv2   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iconv3   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iarfail  = 0       ! Tracks failure of AR model when predicted AR series too large
!!
!! RJMCMC parameters
!!
   INTEGER(KIND=IB),PARAMETER    :: NCHAIN     = 1E9_IB  ! # iterations (max # MCMC steps)
   INTEGER(KIND=IB)              :: ICHAINTHIN = 1E0_IB  ! Chain thinning interval
   INTEGER(KIND=IB)              :: NKEEP      = 1E1_IB  ! Number models to keep before writing
   INTEGER(KIND=IB),PARAMETER    :: NAP        = 10      ! Misc parameters in sample (for bookeeping)
   INTEGER(KIND=IB),PARAMETER    :: NDM        = 100     ! No. steps in lin rot est
   INTEGER(KIND=IB)          :: TCHCKPT              !! No. seconds (integer value) between checkpoints
   INTEGER(KIND=IB)          :: icheckpoint          !! No. of checkpoints to data (read from checkpoint/status.txt)
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)   :: sdevm  ! Std dev for perturbations
!!
!!  Structures for objects and data 
!!
  INTEGER :: imcmc1 = 1   !! Counter for models at T=1 (needs to survive checkpointing!)
  INTEGER :: imcmc2 = 1   !! Counter for mcmc steps to scale diminishing adaptation (needs to survive checkpointing!)
  INTEGER :: NFIELD = 35  !! The number of fields in objstruc (SWD+ELL)
  INTEGER :: objtype1     !! Name of objtype for MPI sending
  INTEGER :: objtype2     !! Name of objtype for MPI sending
  INTEGER :: objtype3     !! Name of objtype for MPI sending
   TYPE :: objstruc
      INTEGER(KIND=IB)                        :: k          ! No. nodes
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: voro       ! 1D voronoi nodes
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:,:):: voroidx    ! 1D voronoi nodes
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: par     ! Forward parameters
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: hiface
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: ziface
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparSWD      !! Std dev SWD data (one per curve)
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparELL      !! Std dev ELL data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveSWD      !! Std dev SWD running average (AR discrimination)
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveELL      !! Std dev ELL running average
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arparSWD      !! AR model parameters SWD
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arparELL      !! AR model parameters ELL
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxarSWD   !! AR on/off index (1=on)
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxarELL   !! AR on/off index (1=on)
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: gvoroidx   !! Index of live parameters on birth/death node
      INTEGER(KIND=IB)                        :: nunique     !! Number of unique interfaces
      INTEGER(KIND=IB)                        :: NFP         !! Number forward parameters
      REAL(KIND=RP)                           :: beta
      REAL(KIND=RP)                           :: logL        !! log likelihood
      REAL(KIND=RP)                           :: logPr       !! log Prior probability ratio
      REAL(KIND=RP)                           :: tcmp
      INTEGER(KIND=IB)                        :: ireject_bd = 0
      INTEGER(KIND=IB)                        :: iaccept_bd = 0
      INTEGER(KIND=IB)                        :: ireject_bds = 0
      INTEGER(KIND=IB)                        :: iaccept_bds = 0
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsSWD     !! Observed data SWD
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredSWD    !! Predicted SWD data for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresSWD     !! SWD data residuals for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarSWD      !! SWD autoregressive model predicted data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: periods     !! SWD periods (s) per curve
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsELL     !! Observed data ELL
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredELL    !! Predicted ELL data for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresELL     !! ELL data residuals for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarELL      !! ELL autoregressive model predicted data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: periods_ELL !! ELL periods (s) per curve
   END TYPE objstruc


!!
!! Structure for covariance matrices (only applies for ICOV >= 2)
!!
!  TYPE :: covstruc
!    REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: Cdi   ! Inverse covariance matrix
!  END TYPE covstruc
!!
!! File units for IO
!!

!!
!! Covariance matrices (only applies for ICOV >= 2)
!!
REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: CdiSWD  ! Inverse covariance matrix SWD data
REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: CdiELL  ! Inverse covariance matrix ELL data
!!

  INTEGER(KIND=IB) :: usample, ustep, ulog
  INTEGER(KIND=IB) :: reclen

   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE      :: icount

!!
!!  Global variables
!!
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE   :: sample                    ! Posterior sample
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE     :: tmpmap                    ! temporary for reading map
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE     :: buf_save_snd,buf_save_rcv ! Buffers for MPI sending
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE  :: buffer1                   !
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE     :: buffer2,buffer3           !

!!
!!  MPI global variables
!!
   INTEGER(KIND=IB)            :: rank,NTHREAD,ierr
   INTEGER(KIND=IB)            :: ncount1,ncount2
   INTEGER(KIND=IB), PARAMETER :: src = 0_IB
   INTEGER                     :: to,from,COMM
   INTEGER                     :: status(MPI_STATUS_SIZE)
   INTEGER(KIND=IB)            :: isize1,isize2,isize3
   INTERFACE
      FUNCTION RANDPERM(num)
         USE data_type, ONLY : IB
         IMPLICIT NONE
         INTEGER(KIND=IB), INTENT(IN) :: num
         INTEGER(KIND=IB), DIMENSION(num) :: RANDPERM
      END FUNCTION RANDPERM
   END INTERFACE
   REAL(KIND=RP) :: tsave1, tsave2              ! Overall time 


  CONTAINS
  !==============================================================================
  integer function newunit(unit)
  !==============================================================================
    integer, intent(out), optional :: unit
    integer, parameter :: LUN_MIN=10, LUN_MAX=1000
    logical :: opened
    integer :: lun
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
      inquire(unit=lun,opened=opened)
      if (.not. opened) then
        newunit=lun
        exit
      end if
    end do
    if (present(unit)) unit=newunit
  end function newunit

END MODULE RJMCMC_COM
!=======================================================================
