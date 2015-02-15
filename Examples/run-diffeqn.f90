MODULE run_diffeqn_a
  USE nag_f77_d_chapter, ONLY: D02PVF, D02PCF

  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! Nag variables
  REAL (KIND=DP), ALLOCATABLE :: work(:), statem(:)
  INTEGER :: lenwrk, ifail

CONTAINS


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_diffeqn
  ! VARIABLES:
  !           state
  !           tmax
  ! SYNOPSIS:
  ! Set up storage and run D02PVF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_diffeqn(state,tmax)
    REAL (KIND=DP), INTENT(IN) ::state(:), tmax
    REAL (KIND=DP) :: tol, thresh(SIZE(state))
    INTEGER :: neq, lenwrk, method, IFAIL


    neq=SIZE(state);  lenwrk=40*neq
    IF (ALLOCATED(work)) DEALLOCATE(work); ALLOCATE(work(lenwrk))
    IF (ALLOCATED(statem)) DEALLOCATE(statem); ALLOCATE(statem(neq))


    method=2
    tol=1.0D-3;  thresh=1.0D-3
    IFAIL=0

    CALL D02PVF(neq, 0.0D0, state, tmax, tol, thresh, method, 'U', .FALSE., &
         & 0.0D0, work, lenwrk, ifail)
    
  END SUBROUTINE setup_diffeqn
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   step_diffeqn
  ! VARIABLES:
  !           state - Must be unchanged outside
  !           time  - Time wanted for next state
  ! SYNOPSIS:
  ! Repeatedly call D02PCF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE step_diffeqn(state,time, diffeqn)
    REAL (KIND=DP), INTENT(INOUT) ::state(:)
    REAL (KIND=DP), INTENT(IN) :: time

    INTERFACE
       SUBROUTINE diffeqn(t,s,sp)
         REAL (KIND(1.0D0)), INTENT(IN) ::t,s(SIZE(state))
         REAL (KIND(1.0D0)), INTENT(OUT) ::sp(SIZE(STATE))
       end SUBROUTINE diffeqn
    end INTERFACE


    REAL (KIND=DP) :: tgot, dstate(SIZE(state))
    INTEGER :: IFAIL

    ! The following just ensures that we don't give up merely because
    ! it took more than 200 steps, or whathever unavoidable limit
    ! nag imposes.  The two error codes are "A considerable amount of
    ! work has been expended since the last call" and "the problem
    ! appears to be stiff"
    IFAIL=1 
   CALL D02PCF(diffeqn, time, tgot, state, dstate, &
         & statem, work, ifail)

    DO WHILE ((IFAIL .EQ. 3) .OR. (IFAIL .EQ. 4)) 
       IFAIL=1
       CALL D02PCF(diffeqn, time, tgot, state, dstate, &
            & statem, work, ifail)
    end DO

    IF (IFAIL .NE. 0) THEN
       WRITE(*,'("**** D02PCF failed with error code",I5,X,F8.5" ****")') &
            & IFAIL, time
       STOP
    end IF


  END SUBROUTINE step_diffeqn

END MODULE run_diffeqn_a

MODULE run_diffeqn_b
  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  INTEGER :: iwork, n
  REAL (KIND=DP), ALLOCATABLE :: work(:)
  REAL (KIND=DP) :: tlast

  ! Alternate NAG drivers

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_diffeqn
  ! VARIABLES:
  !           state
  !           tmax
  ! SYNOPSIS:
  ! Set up storage and run D02PVF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_diffeqn(state,tmax)
    REAL (KIND=DP), INTENT(IN) ::state(:), tmax

    n=SIZE(state)

    iwork=28+21*n

    ALLOCATE(work(iwork))

    ! REcord initial time (And avoid error on unused tmax)
    tlast=tmax; tlast=0.0D0

  end SUBROUTINE setup_diffeqn
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   step_diffeqn
  ! VARIABLES:
  !           state - Must be unchanged outside
  !           time  - Time wanted for next state
  ! SYNOPSIS:
  ! Repeatedly call D02PCF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE step_diffeqn(state,time, diffeqn)
    USE nag_f77_d_chapter, ONLY: d02cjf
    REAL (KIND=DP), INTENT(INOUT) ::state(:)
    REAL (KIND=DP), INTENT(IN) :: time

    INTERFACE
       SUBROUTINE diffeqn(t,s,sp)
         REAL (KIND(1.0D0)), INTENT(IN) ::t,s(SIZE(state))
         REAL (KIND(1.0D0)), INTENT(OUT) ::sp(SIZE(STATE))
       end SUBROUTINE diffeqn
    end INTERFACE

    EXTERNAL :: d02cjx
    REAL (KIND=DP), EXTERNAL :: d02cjw

    INTEGER :: ifail
    REAL (KIND=DP) :: tol

    tol=1.0D-3
    ifail=0

    CALL d02CJF(tlast, time, n, state, diffeqn, &
         & tol, 'M', D02CJX, D02CJW, work, ifail)

  end SUBROUTINE step_diffeqn
  
  
END MODULE run_diffeqn_b
