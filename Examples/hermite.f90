MODULE hermite
  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  REAL (KIND=DP), ALLOCATABLE, TARGET :: Coeff(:,:), norm(:)
  

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   SetupHermiteCoeff
  ! VARIABLES:
  !           nmax - Maximum order
  ! SYNOPSIS:
  ! Set up coefficients
  ! http://mathworld.wolfram.com/HermitePolyCoeffnomial.html
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SetupHermiteCoeff(nmax)
    USE utility, ONLY: pi

    INTEGER, INTENT(IN) ::nmax

    INTEGER :: n, i
    REAL (KIND=DP), ALLOCATABLE ::hp1(:), hp2(:)
    REAL (KIND=DP), POINTER :: hp(:)
    REAL (KIND=DP) :: nfac
    
    ALLOCATE(Coeff(0:nmax, nmax+1), norm(0:nmax))
    
    DO n=0, nmax

       ALLOCATE(hp1(n+1), hp2(n+1))

       hp=>COEFF(n,1:n+1)
       
       IF(n .EQ. 0) THEN
          hp(1)=1.0
       ELSE IF(n .EQ. 1) THEN
          hp(1)=2.0
          hp(2)=0.0
       ELSE
          hp1(1:n+1)=0.0
          hp1(1:n)=2.0*Coeff(n-1, 1:n)

          hp2(1:n+1)=0.0
          hp2(3:)=2.0*(n-1)*Coeff(n-2, 1:n-1)

          hp=hp1-hp2

       END IF

       DEALLOCATE(hp1,hp2)

       nfac=1.0D0
       DO i=1, n
          nfac=nfac*i
       end DO

       norm(n) = sqrt(pi) *nfac * (2.0D0**n)

    end DO

  end SUBROUTINE SetupHermiteCoeff
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   wf
  ! VARIABLES:
  !           x - Coordinate
  !           n  - WF index
  ! SYNOPSIS:
  ! Evaluate HO wf at position x.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION wf(x,n)

    REAL (KIND=DP), INTENT(IN) ::x
    INTEGER, INTENT(IN) ::n

    REAL (KIND=DP) :: wf, parts(0:n)
    REAL (KIND=DP), POINTER :: hp(:)
    INTEGER :: i

    hp=>Coeff(n, 1:n+1)

    DO i=0, n
       parts(i) = hp(n+1-i) * (x**i)
    end DO


    wf=SUM(parts)*exp(-0.5*x*x)/sqrt(norm(n))

  END FUNCTION wf


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   matrix_density, vector_density
  ! VARIABLES:
  !           nmat  - Photon density matrix
  !           x     -  Position
  ! Synopsis:
  ! Calculate the density from the matrix of photons, or vector
  ! for gain profile.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION matrix_density(nmat,x)
    COMPLEX (KIND=DP), INTENT(IN) :: nmat(0:, 0:)
    REAL (KIND=DP), INTENT(IN) :: x
    REAL (KIND=DP) :: wfvec(0:SIZE(nmat,1)-1)
    REAL (KIND=DP) ::matrix_density

    INTEGER :: n

    wfvec = (/ (wf(x,n), n=0, SIZE(nmat,1)-1) /)

    ! Note lack of complex wavefunctions.
    matrix_density = REAL(DOT_PRODUCT(wfvec, MATMUL(nmat, wfvec)))
    
  END FUNCTION matrix_density

  FUNCTION vector_density(uvec,x)
    REAL (KIND=DP), INTENT(IN) :: uvec(0:)
    REAL (KIND=DP), INTENT(IN) :: x
    REAL (KIND=DP) :: wfvec(0:SIZE(uvec)-1)
    REAL (KIND=DP) ::vector_density

    INTEGER :: g

    wfvec = (/ (wf(x,g), g=0, SIZE(uvec)-1) /)

    ! Note lack of complex wavefunctions.
    vector_density = SUM(wfvec*uvec)
    
  END FUNCTION vector_density

END MODULE hermite
