MODULE rs_basis
  USE params, ONLY: nmodes

  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! Real space for excitation profile: pre-calculates the matrix elements.
  
  ! Number of excitation modes (points)
  INTEGER :: gmodes

  ! Consequent requires states
  INTEGER :: statesize

  ! Storage of values
  REAL (KIND=DP), ALLOCATABLE :: avec(:,:,:)
  REAL (KIND=DP) :: dx, x0

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_basis_overlaps
  ! SYNOPSIS:
  ! Set up, or read in, as appropriate, overlaps between basis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_basis_overlaps()
    USE hermite, ONLY: SetupHermiteCoeff
    USE utility, ONLY: parse_logical
    

    CHARACTER :: fname*128, gmodestr*3, nmodestr*3,dxstr*8, readstr

    LOGICAL :: do_read, exists

    WRITE(*,*) "Number, step: dx.  Read from file if possible?"
    READ(*,*) gmodes, dx, readstr
    do_read=parse_logical(readstr)

    ! Calculate the offset, so that x0 + (gmodes-1)*dx =  - x0
    ! Note that gmode sis really the number of points, i.e we index
    ! from 1.
    x0 = - 0.5*(gmodes-1)*dx

    ! Consequent size of state (noting that nmatrix is complex).
    statesize = gmodes + 2*(nmodes+1)*(nmodes+1)


    ! Create the filename for current pattern of modes.  If desired,
    ! we read from this.  Always we write to this.
    WRITE(gmodestr,'(I3)') gmodes
    WRITE(dxstr, '(F8.5)') dx
    WRITE(nmodestr,'(I3)') nmodes

    fname="rs_basis"//TRIM(ADJUSTL(nmodestr))//"_"//&
         &            TRIM(ADJUSTL(gmodestr))//"_dx="// &
         &            TRIM(ADJUSTL(dxstr))//".dat"

    ALLOCATE(avec(0:nmodes,0:nmodes,gmodes))
    
    ! First make sure the hermitie coefficients are donw -- we need
    ! this even if reusing overlaps as we use it in calculating real
    ! space representations.
    CALL SetupHermiteCoeff(MAX(nmodes,gmodes))


    ! See whether file eixsts to read things from.
    INQUIRE(FILE=fname, EXIST=exists)
    IF (do_read .AND. exists) THEN
       WRITE(*,*) "Reading overlaps from file"
       CALL read_basis_overlaps(fname)
    ELSE
       WRITE(*,*) "Calculating overlaps"
       CALL calculate_basis_overlaps()
    end IF
    
    ! ALways write out results
    CALL write_basis_overlaps(fname)
    
  END SUBROUTINE setup_basis_overlaps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   calculate_basis_overlaps
  ! SYNOPSIS:
  ! Calculate all basis overlaps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calculate_basis_overlaps()
    USE hermite, ONLY: wf

    INTEGER :: n1, n2, g
    REAL (KIND=DP) :: x

    DO n1=0, nmodes
       DO n2=0, nmodes
          
          DO g=1, gmodes
             x = x0 + (g-1)*dx

             avec(n1,n2,g) = wf(x, n1) * wf(x,n2)

          end DO
          
       end DO
    end DO
    

  end SUBROUTINE calculate_basis_overlaps
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   write_basis_overlaps
  ! VARIABLES:
  !           fname - Filename to write to
  ! SYNOPSIS:
  ! Write out values
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE write_basis_overlaps(fname)
    CHARACTER, INTENT(IN) ::fname*128

    INTEGER :: OUT=27

    OPEN(OUT,FILE=fname)
    WRITE(OUT,*) nmodes, gmodes
    WRITE(OUT,*) dx
    WRITE(OUT,*) avec
    CLOSE(OUT)
    
  END SUBROUTINE write_basis_overlaps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   read_basis_overlaps
  ! VARIABLES:
  !           fname - Filename to read to
  ! SYNOPSIS:
  ! Read in values
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_basis_overlaps(fname)
    USE utility, ONLY: assert
    CHARACTER, INTENT(IN) ::fname*128

    INTEGER :: IN=27, d_nmodes, d_gmodes

    OPEN(IN,FILE=fname)
    READ(IN,*) d_nmodes, d_gmodes
    CALL assert(nmodes.EQ.d_nmodes, "Mismatched nmodes")
    CALL assert(gmodes.EQ.d_gmodes, "Mismatched nmodes")
    READ(IN,*) dx
    READ(IN,*) avec
    CLOSE(IN)
    
  END SUBROUTINE read_basis_overlaps


END MODULE rs_basis
