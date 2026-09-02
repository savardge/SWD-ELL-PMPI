!=======================================================================

SUBROUTINE MAKE_MPI_STRUC_SP(obj,objtype)
!!=======================================================================
!!
!! MPI structured type mirroring objstruc (35 fields, SWD+ELL).
!! Offsets come from mpi_get_address; note that objstruc components are
!! therefore addressed ABSOLUTELY -- never re-allocate them (see COPY_OBJ).
!!
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc) :: obj
INTEGER(KIND=IB) :: ifield,objtype
INTEGER(KIND=IB) :: oldtypes(NFIELD), blockcounts(NFIELD)
INTEGER(KIND=MPI_ADDRESS_KIND) :: offsets(NFIELD)

!!                k   voro     voroidx           par          hiface     ziface
blockcounts = (/ 1,  NLMX*NPL, NLMX*NPL , (NLMX+1)*NPL*NPL,  NPL*NLMX,   NPL*NLMX,&
!!               sdparSWD  sdparELL   sdaveSWD  sdaveELL   arparSWD arparELL  idxarSWD  idxarELL  gvoroidx nunique NFP
                 NMODE,   NMODE_ELL,   NMODE,   NMODE_ELL,   NMODE,  NMODE_ELL,  NMODE,  NMODE_ELL,   NPL-1,    1,    1, &
!!               beta    logL     logPr tcmp  ireject_bd  iaccept_bd ireject_bds  iaccept_bds
                  1,       1,       1,    1,      1,          1,         1,          1,&
!!                DobsSWD           DpredSWD        DresSWD       DarSWD           periods
                 NMODE*NDAT_SWD, NMODE*NDAT_SWD, NMODE*NDAT_SWD, NMODE*NDAT_SWD, NMODE*NDAT_SWD,&
!!                      DobsELL           DpredELL            DresELL             DarELL           periods_ELL
                 NMODE_ELL*NDAT_ELL, NMODE_ELL*NDAT_ELL, NMODE_ELL*NDAT_ELL, NMODE_ELL*NDAT_ELL, NMODE_ELL*NDAT_ELL /)

oldtypes(1)     = MPI_INTEGER
oldtypes(2)     = MPI_DOUBLE_PRECISION
oldtypes(3)     = MPI_INTEGER
oldtypes(4:12)  = MPI_DOUBLE_PRECISION
oldtypes(13:17) = MPI_INTEGER
oldtypes(18:21) = MPI_DOUBLE_PRECISION
oldtypes(22:25) = MPI_INTEGER
oldtypes(26:35) = MPI_DOUBLE_PRECISION

call mpi_get_address(obj%k,offsets(1),ierr)
call mpi_get_address(obj%voro,offsets(2),ierr)
call mpi_get_address(obj%voroidx,offsets(3),ierr)
call mpi_get_address(obj%par,offsets(4),ierr)
call mpi_get_address(obj%hiface,offsets(5),ierr)
call mpi_get_address(obj%ziface,offsets(6),ierr)
call mpi_get_address(obj%sdparSWD,offsets(7),ierr)
call mpi_get_address(obj%sdparELL,offsets(8),ierr)
call mpi_get_address(obj%sdaveSWD,offsets(9),ierr)
call mpi_get_address(obj%sdaveELL,offsets(10),ierr)
call mpi_get_address(obj%arparSWD,offsets(11),ierr)
call mpi_get_address(obj%arparELL,offsets(12),ierr)
call mpi_get_address(obj%idxarSWD,offsets(13),ierr)
call mpi_get_address(obj%idxarELL,offsets(14),ierr)
call mpi_get_address(obj%gvoroidx,offsets(15),ierr)
call mpi_get_address(obj%nunique,offsets(16),ierr)
call mpi_get_address(obj%NFP,offsets(17),ierr)
call mpi_get_address(obj%beta,offsets(18),ierr)
call mpi_get_address(obj%logL,offsets(19),ierr)
call mpi_get_address(obj%logPr,offsets(20),ierr)
call mpi_get_address(obj%tcmp,offsets(21),ierr)
call mpi_get_address(obj%ireject_bd,offsets(22),ierr)
call mpi_get_address(obj%iaccept_bd,offsets(23),ierr)
call mpi_get_address(obj%ireject_bds,offsets(24),ierr)
call mpi_get_address(obj%iaccept_bds,offsets(25),ierr)
call mpi_get_address(obj%DobsSWD,offsets(26),ierr)
call mpi_get_address(obj%DpredSWD,offsets(27),ierr)
call mpi_get_address(obj%DresSWD,offsets(28),ierr)
call mpi_get_address(obj%DarSWD,offsets(29),ierr)
call mpi_get_address(obj%periods,offsets(30),ierr)
call mpi_get_address(obj%DobsELL,offsets(31),ierr)
call mpi_get_address(obj%DpredELL,offsets(32),ierr)
call mpi_get_address(obj%DresELL,offsets(33),ierr)
call mpi_get_address(obj%DarELL,offsets(34),ierr)
call mpi_get_address(obj%periods_ELL,offsets(35),ierr)

DO ifield=2,SIZE(offsets)
  offsets(ifield) = offsets(ifield) - offsets(1)
ENDDO
offsets(1) = 0

call MPI_TYPE_CREATE_STRUCT( NFIELD, blockcounts, offsets, oldtypes, objtype,ierr)
call MPI_TYPE_COMMIT(objtype, ierr)

RETURN
END SUBROUTINE MAKE_MPI_STRUC_SP
!!==============================================================================

SUBROUTINE ALLOC_OBJ(obj)
!!==============================================================================
!!
!! Allocates memory (every component, unconditionally, so that COPY_OBJ can
!! always do plain value assignment).
!!
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc) :: obj

ALLOCATE( obj%voro(NLMX,NPL),obj%voroidx(NLMX,NPL),obj%par((NLMX+1)*NPL*NPL) )
ALLOCATE( obj%hiface(NLMX*NPL),obj%ziface(NLMX*NPL) )
ALLOCATE( obj%sdaveSWD(NMODE),obj%sdaveELL(NMODE_ELL) )
ALLOCATE( obj%sdparSWD(NMODE),obj%arparSWD(NMODE),obj%idxarSWD(NMODE) )
ALLOCATE( obj%sdparELL(NMODE_ELL),obj%arparELL(NMODE_ELL),obj%idxarELL(NMODE_ELL) )
ALLOCATE( obj%gvoroidx(NPL-1) )
!! SWD data:
ALLOCATE( obj%DobsSWD(NMODE,NDAT_SWD),obj%DpredSWD(NMODE,NDAT_SWD),obj%DresSWD(NMODE,NDAT_SWD) )
ALLOCATE( obj%DarSWD(NMODE,NDAT_SWD),obj%periods(NMODE,NDAT_SWD) )
!! ELL data:
ALLOCATE( obj%DobsELL(NMODE_ELL,NDAT_ELL),obj%DpredELL(NMODE_ELL,NDAT_ELL),obj%DresELL(NMODE_ELL,NDAT_ELL) )
ALLOCATE( obj%DarELL(NMODE_ELL,NDAT_ELL),obj%periods_ELL(NMODE_ELL,NDAT_ELL) )

obj%voro     = 0._RP
obj%voroidx  = 0
obj%par      = 0._RP
obj%hiface   = 0._RP
obj%ziface   = 0._RP
obj%sdparSWD = 0._RP
obj%sdparELL = 0._RP
obj%sdaveSWD = 0._RP
obj%sdaveELL = 0._RP
obj%arparSWD = 0._RP
obj%arparELL = 0._RP
obj%idxarSWD = 0
obj%idxarELL = 0
obj%gvoroidx = 0
obj%DobsSWD  = 0._RP
obj%DpredSWD = 0._RP
obj%DresSWD  = 0._RP
obj%DarSWD   = 0._RP
obj%periods  = 0._RP
obj%DobsELL  = 0._RP
obj%DpredELL = 0._RP
obj%DresELL  = 0._RP
obj%DarELL   = 0._RP
obj%periods_ELL  = 0._RP
END SUBROUTINE ALLOC_OBJ
!=======================================================================

SUBROUTINE ALLOC_COVmat()
!!==============================================================================

USE RJMCMC_COM
IMPLICIT NONE

IF (I_SWD == 1) THEN
    IF (ICOV_SWD==2) ALLOCATE( CdiSWD(NDAT_SWD,NDAT_SWD) )
    IF (ICOV_iterUpdate_SWD==1) THEN
        ALLOCATE( CdSWD(NDAT_SWD,NDAT_SWD) )
        ALLOCATE( CdSWD_old(NDAT_SWD,NDAT_SWD) )
        IF (.NOT.ALLOCATED(CdiSWD)) ALLOCATE( CdiSWD(NDAT_SWD,NDAT_SWD) )
        ALLOCATE( CdiSWD_old(NDAT_SWD,NDAT_SWD) )
    END IF
END IF !!I_SWD

IF (I_ELL == 1) THEN
    IF (ICOV_ELL==2) ALLOCATE( CdiELL(NDAT_ELL,NDAT_ELL) )
    IF (ICOV_iterUpdate_ELL==1) THEN
        ALLOCATE( CdELL(NDAT_ELL,NDAT_ELL) )
        ALLOCATE( CdELL_old(NDAT_ELL,NDAT_ELL) )
        IF (.NOT.ALLOCATED(CdiELL)) ALLOCATE( CdiELL(NDAT_ELL,NDAT_ELL) )
        ALLOCATE( CdiELL_old(NDAT_ELL,NDAT_ELL) )
    END IF
END IF !!I_ELL

END SUBROUTINE ALLOC_COVmat
!!=============================================================================
