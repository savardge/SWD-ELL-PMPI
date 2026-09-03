!! Standalone driver around the production SWD forward model, for validation
!! against an independent implementation (disba).
!!
!! Calls dispersion() exactly as LOGLHOOD_SWD does -- same warm start, same
!! per-period validity flags -- once with IGRP=0 (phase) and once with IGRP=1
!! (group), for every requested mode.
!!
!! Usage:
!!   ./disp_driver models.txt periods.txt maxmode cmin cmax dc dc_over iwarm
!!
!! models.txt : repeated blocks of
!!                nl
!!                thick rho vp vs      (nl lines; km, g/cm3, km/s, km/s)
!! periods.txt: one period [s] per line, ASCENDING (the warm start requires it)
!!
!! stdout (CSV, one line per model/mode/period):
!!   imodel,mode,iperiod,period,ivalid,c,u
!! with c,u = 0 and ivalid = 0 where the mode has no trapped root.
PROGRAM disp_driver
  IMPLICIT NONE
  INTEGER, PARAMETER :: NTMAX = 500, NLMAX = 500, NMODMAX = 100000
  REAL    :: peri(NTMAX), cph(NTMAX), cgr(NTMAX)
  INTEGER :: ivalid(NTMAX), ivalidg(NTMAX)
  REAL    :: thick(NLMAX), rho(NLMAX), vp(NLMAX), vs(NLMAX)
  REAL    :: cmin, cmax, dc, dcov
  INTEGER :: nt, nl, imod, iper, ier, mode, maxmode, iwarm, i
  CHARACTER(LEN=256) :: fmod, fper, arg

  CALL GET_COMMAND_ARGUMENT(1, fmod)
  CALL GET_COMMAND_ARGUMENT(2, fper)
  CALL GET_COMMAND_ARGUMENT(3, arg); READ(arg,*) maxmode
  CALL GET_COMMAND_ARGUMENT(4, arg); READ(arg,*) cmin
  CALL GET_COMMAND_ARGUMENT(5, arg); READ(arg,*) cmax
  CALL GET_COMMAND_ARGUMENT(6, arg); READ(arg,*) dc
  CALL GET_COMMAND_ARGUMENT(7, arg); READ(arg,*) dcov
  CALL GET_COMMAND_ARGUMENT(8, arg); READ(arg,*) iwarm

  OPEN(10, FILE=fper, STATUS='OLD', ACTION='READ')
  nt = 0
  DO
    READ(10,*,IOSTAT=ier) peri(nt+1)
    IF (ier /= 0) EXIT
    nt = nt + 1
  ENDDO
  CLOSE(10)

  WRITE(*,'(A)') 'imodel,mode,iperiod,period,ivalid,c,u'

  OPEN(11, FILE=fmod, STATUS='OLD', ACTION='READ')
  DO imod = 1, NMODMAX
    READ(11,*,IOSTAT=ier) nl
    IF (ier /= 0) EXIT
    DO i = 1, nl
      READ(11,*) thick(i), rho(i), vp(i), vs(i)
    ENDDO

    DO mode = 0, maxmode
      !! Phase velocity (IGRP = 0)
      CALL dispersion(nl, rho, vp, vs, thick, cph, peri, nt, 0, ier, mode, &
                      ivalid, cmin, cmax, dc, dcov, iwarm)
      !! Group velocity (IGRP = 1) -- same root, same scan, different output
      CALL dispersion(nl, rho, vp, vs, thick, cgr, peri, nt, 1, ier, mode, &
                      ivalidg, cmin, cmax, dc, dcov, iwarm)

      DO iper = 1, nt
        !! The two passes must agree on which periods carry a root; if they
        !! ever disagree, mark the period invalid so the comparison sees it.
        IF (ivalid(iper) /= ivalidg(iper)) ivalid(iper) = 0
        IF (ivalid(iper) == 0) THEN
          cph(iper) = 0.
          cgr(iper) = 0.
        ENDIF
        WRITE(*,'(I0,",",I0,",",I0,",",ES16.9,",",I0,",",ES16.9,",",ES16.9)') &
              imod, mode, iper, peri(iper), ivalid(iper), cph(iper), cgr(iper)
      ENDDO
    ENDDO
  ENDDO
  CLOSE(11)
END PROGRAM disp_driver
