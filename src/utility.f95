MODULE utility
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=KIND(1.0D0)

  ! Utility functions, no specific connection to current problem, or
  ! general mathematics.
  REAL (KIND=DP), PARAMETER :: pi=3.14159265358979323844D0
  COMPLEX (KIND=DP), PARAMETER :: ii=(0.0D0, 1.0D0), one=(1.0D0, 0.0D0)
  
CONTAINS



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   evaluate_integral
  ! VARIABLES:
  !           type - Vector of Matrix
  ! SYNOPSIS:
  ! Evaluate an element of the overlap matrix
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION evaluate_integral(integrand)
    REAL (KIND=DP) ::evaluate_integral
    INTERFACE
       FUNCTION integrand(x)
         REAL (KIND=KIND(1.0D0)), INTENT(IN) ::x
         REAL (KIND=KIND(1.0D0)) ::integrand
       END FUNCTION integrand
    end INTERFACE
    
    REAL (KIND=DP) :: epsabs, epsrel, abserr
    INTEGER :: ifail
    INTEGER, PARAMETER :: lw=2000, liw=500
    INTEGER :: iwork(liw)
    REAL (KIND=DP) :: work(lw)


    epsabs=1.0D-5; epsrel=1.0D-4;

    CALL d01AMF(integrand, 0.0D0, 2, epsabs, epsrel, &
         & evaluate_integral, abserr, work, lw, iwork, liw, ifail)

  END FUNCTION evaluate_integral


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   die
  ! VARIABLES:
  !           msg
  ! SYNOPSIS:
  ! Die with given message
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE die(msg)
    CHARACTER (LEN=*), INTENT(IN) ::msg

    WRITE(*,*) msg
    STOP
    
  END SUBROUTINE die

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   assert
  ! SYNOPSIS:
  ! Assert that a condition is true, or die with given error message
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE assert(condition,error)
    LOGICAL, INTENT(IN) ::condition
    CHARACTER (LEN=*), INTENT(IN) ::error

    IF (.NOT. condition) THEN
       WRITE(*,*) error
       STOP
    end IF

  END SUBROUTINE assert


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   get_range
  ! VARIABLES:
  !           min
  !           max
  !           steps
  !           prefix
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_range(prefix, min,max,steps)
    REAL (KIND=DP), INTENT(OUT) ::min,max
    INTEGER, INTENT(OUT) ::steps
    CHARACTER, INTENT(IN) ::prefix

    WRITE(*,*) prefix, " (Min, Max, Steps)"
    READ(*,*) min, max, steps

  END SUBROUTINE get_range

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   interpolate
  ! VARIABLES:
  !           min, max - Range of variable
  !           i        - Index
  !           steps    - Number of points
  !           log      - Whether to use a log scale
  ! SYNOPSIS:
  !   Return a number interpolated between min and max, according to
  !   i/steps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION interpolate(min,max,i,steps,logarithmic)
    REAL (KIND=DP), INTENT(IN) ::min,max
    INTEGER, INTENT(IN) ::i,steps
    REAL (KIND=DP) ::interpolate
    LOGICAL :: logarithmic

    IF (steps .EQ. 1) THEN
       interpolate=min
    ELSE
       IF (logarithmic) THEN
          interpolate = min * exp( log(max/min)*(i-1)/(steps-1))
       ELSE
          interpolate=min + (max-min)*(i-1)/(steps-1)
       end IF
    end IF

    IF (i .EQ. 1) THEN
       interpolate = min
    ELSE IF (i .EQ. steps) THEN
       interpolate = max
    end IF
    

  END FUNCTION interpolate


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   read_logical
  ! SYNOPSIS:
  ! Read a string and test for y/N
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION read_logical()
    LOGICAL (KIND=DP) ::read_logical
    CHARACTER :: str*4
    
    READ(*,*) str

    read_logical=parse_logical(str)

  END FUNCTION read_logical
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   parse_logical
  ! VARIABLES:
  !           str
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION parse_logical(str)
    CHARACTER, INTENT(IN) ::str
    LOGICAL (KIND=DP) ::parse_logical

    parse_logical=(SCAN(str,'yYtT') .GT. 0) 
    
  END FUNCTION parse_logical



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   dinvnorm
  ! VARIABLES: (z)
  ! SYNOPSIS:
  ! ren-raw chen, rutgers business school
  ! normal inverse
  ! translate from
  ! http://home.online.no/~pjacklam/notes/invnorm
  ! a routine written by john herrero
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION dinvnorm(p)
    REAL (KIND=DP), INTENT(IN) :: p

    REAL (KIND=DP) ::dinvnorm, q, r

    REAL (KIND=DP), PARAMETER ::&
         & a1=-39.6968302866538, &
         & a2=220.946098424521, &
         & a3=-275.928510446969, &
         & a4=138.357751867269, &
         & a5=-30.6647980661472, &
         & a6=2.50662827745924, &
         & b1=-54.4760987982241, &
         & b2=161.585836858041, &
         & b3=-155.698979859887, &
         & b4=66.8013118877197, &
         & b5=-13.2806815528857, &
         & c1=-0.00778489400243029, &
         & c2=-0.322396458041136, &
         & c3=-2.40075827716184, &
         & c4=-2.54973253934373, &
         & c5=4.37466414146497, &
         & c6=2.93816398269878, &
         & d1=0.00778469570904146, &
         & d2=0.32246712907004, &
         & d3=2.445134137143, &
         & d4=3.75440866190742, &
         & p_low=0.02425
    REAL (KIND=DP) ::  p_high
    p_high=1-p_low

    CALL ASSERT((p .GE. 0.0D0) .AND. (p .LT. 1.0D0),&
         &  "Invalid range for argument to ")

    IF(p.LT.p_low) THEN
       q=sqrt(-2*log(p))
       dinvnorm=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/&
            & ((((d1*q+d2)*q+d3)*q+d4)*q+1)
    ELSE IF(p.le.p_high) THEN
       q=p-0.5; r=q*q
       dinvnorm=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/&
            &(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
    ELSE 
       q=sqrt(-2*log(1-p))
       dinvnorm=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/&
            &((((d1*q+d2)*q+d3)*q+d4)*q+1)
    end IF

  END FUNCTION dinvnorm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   DMATMUL
  ! VARIABLES:
  !           mat   - Matrix
  !           diag  - Diagonal part
  !           type  - 'L' for pre, 'R' for post
  ! SYNOPSIS:
  ! PRe- or post-multiply matrix by diagonal matrix.  Note that
  ! within this function we do not need to know what the indexing of
  ! mat and diag is, so the following is written with standard fortran
  ! indexing, i.e. start at 1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION DMATMUL(mat,diag,type) RESULT(OUT)
    COMPLEX (KIND=DP), INTENT(IN) ::mat(:,:)
    REAL (KIND=DP) :: diag(:)
    COMPLEX (KIND=DP) :: out(SIZE(mat,1),SIZE(mat,2))
    CHARACTER :: type

    INTEGER :: n

    out=0.0D0

    SELECT CASE(type)
    case('L')
       DO n=1, SIZE(diag)
          out(n,:) = diag(n) * mat(n,:)
       end DO
    case('R')
       DO n=1, SIZE(diag)
          out(:,n) = mat(:,n) * diag(n)
       end DO
    CASE DEFAULT
       WRITE(*,*) "Invalid type in multiplication"
       STOP
    end SELECT


  end FUNCTION DMATMUL


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! NAME:   HERMITIAN
    ! VARIABLES:
    !           matrix
    ! SYNOPSIS:
    ! Hermitize:  (a+a^dagger)/2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION HERMITIAN(matrix)
      COMPLEX (KIND=DP), INTENT(IN) ::matrix(:,:)
      COMPLEX (KIND=DP) ::HERMITIAN(SIZE(matrix,1),SIZE(matrix,2))

      hermitian=0.5D0*(matrix + CONJG(TRANSPOSE(matrix)))

    END FUNCTION HERMITIAN

END MODULE utility
