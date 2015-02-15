PROGRAM ft
  USE nag_library, ONLY: c06ecf

  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  REAL (KIND=DP), PARAMETER :: pi=3.14159265358979323844D0


  REAL (KIND=DP) :: dt, t, dummy(3), w, dw
  REAL (KIND=DP), ALLOCATABLE :: rl(:), il(:)
  COMPLEX (KIND=DP) :: lam

  INTEGER :: n,l, ifail
  INTEGER :: IN=27, OUT=28
  CHARACTER :: input_filename*128, output_filename*128
  !
  WRITE(*,*) "Filename, number of points"
  READ(*,*) input_filename, n
  WRITE(*,*) "Output filenamne"
  READ(*,*) output_filename

  ALLOCATE(rl(n), il(n))

  OPEN(IN,file=input_filename)
  READ(IN,*)

  DO l=1, n
     READ(IN,*) t, rl(l), il(l), dummy(1:3)
  end DO
  CLOSE(IN)

  ! Use last time to give total time
  dw = 2*pi/t
  dt = t/(N-1)

  ifail=0
  CALL c06ecf(rl,il,n,ifail)

  OPEN(OUT, file=output_filename)
  ! Write out in order of increaseing frequency, with zero
  ! near the middle
  DO l=1+n/2, n-1
     lam = dt*sqrt(1.0D0*N)*CMPLX(rl(l),il(l),KIND(1.0D0))
     
     w=dw*(l-n)
     WRITE(OUT,'(4(D12.5,X))') w, ABS(lam)**2, lam
  end DO
  
  DO l=0, n/2
     lam = dt*sqrt(1.0D0*N)*CMPLX(rl(l),il(l),KIND(1.0D0))

     w=dw*l
     WRITE(OUT,'(4(D12.5,X))') w, ABS(lam)**2, lam
  end DO
  CLOSE(OUT)
  
END PROGRAM ft
