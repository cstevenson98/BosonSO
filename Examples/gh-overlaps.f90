MODULE gh_overlaps
  USE params, ONLY: nmodes
  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! Deals with the overlaps between wavefunctions
  
  ! Number of  excitation modes.  
  INTEGER :: gmodes

  ! Consequent requires states
  INTEGER :: statesize

  ! Storage of overlaps
  REAL (KIND=DP), ALLOCATABLE :: a0vec(:), avec(:,:,:), amat(:,:,:,:)

  ! Private storage of which integral we are evaluating
  INTEGER, PRIVATE :: n1, n2, g1, g2
  CHARACTER, PRIVATE :: int_type

  LOGICAL :: DEBUG=.TRUE.

CONTAINS
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   overlap_integrand
  ! VARIABLES:
  !           x
  ! SYNOPSIS:
  ! Basic integrand function for overlaps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION overlap_integrand(x)
    USE hermite, ONLY: wf

    REAL (KIND=DP), INTENT(IN) ::x
    REAL (KIND=DP) ::overlap_integrand
    
    SELECT CASE(int_type) 
    CASE('0')
       ! Flat profile overlap 
       overlap_integrand=wf(x,g1)
    CASE('T')
       ! Test 
       overlap_integrand=wf(x,n1)*wf(x,n2)
    CASE('V')
       ! Vector of matrices
       overlap_integrand=wf(x,n1)*wf(x,n2)*wf(x,g1)
    CASE('M')
       ! Matrix of matrices
       overlap_integrand=wf(x,n1)*wf(x,n2)*wf(x,g1)*wf(x,g2)
    end SELECT
    
  END FUNCTION overlap_integrand

  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   evaluate_overlap
  ! VARIABLES:
  !           type
  ! SYNOPSIS:
  ! Evaluate elements of the overlap matrix
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION evaluate_overlap(type)
    USE utility, ONLY: die, evaluate_integral

    CHARACTER, INTENT(IN) ::type
    REAL (KIND=DP) ::evaluate_overlap

    INTEGER :: pcheck
    
    ! Set module level type
    int_type=type

    SELECT CASE(int_type)
    CASE('0')
       pcheck=g1
    CASE('T')
       pcheck=n1+n2
    CASE('V')
       pcheck=n1+n2+g1
    CASE('M')
       pcheck=n1+n2+g1+g2
    CASE DEFAULT
       CALL die("Invalid type of integral requested")
    end SELECT
    
    ! NB: Could check parity, but not really worth it for factor
    ! of two in setup.
    evaluate_overlap = evaluate_integral(overlap_integrand)
    

  END FUNCTION evaluate_overlap



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_basis_overlaps
  ! SYNOPSIS:
  ! Set up, or read in, as appropriate, overlaps between basis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_basis_overlaps()
    USE hermite, ONLY: SetupHermiteCoeff
    USE utility, ONLY: parse_logical

    CHARACTER :: fname*128, gmodestr*3, nmodestr*3, readstr

    LOGICAL :: do_read, exists

    WRITE(*,*) "Maximum excitation basis states.  Read from file if possible?"
    READ(*,*) gmodes, readstr
    do_read=parse_logical(readstr)
    
    ! Consequent size of state (noting that nmatrix is complex).
    statesize = (gmodes+1) + 2*(nmodes+1)*(nmodes+1)

    ! Create the filename for current pattern of modes.  If desired,
    ! we read from this.  Always we write to this.
    WRITE(gmodestr,'(I3)') gmodes
    WRITE(nmodestr,'(I3)') nmodes

    fname="basis_overlaps_"//TRIM(ADJUSTL(nmodestr))//"_"//&
         &                   TRIM(ADJUSTL(gmodestr))//".dat"


    ALLOCATE(a0vec(0:gmodes),&
         &   avec(0:nmodes,0:nmodes,0:gmodes), &
         &   amat(0:nmodes,0:nmodes,0:gmodes,0:gmodes))


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
    REAL (KIND=DP) :: test

    DO g1=0, gmodes
       a0vec(g1) = evaluate_overlap('0')
    end DO
   

    DO n1=0, nmodes
       DO n2=0, n1

          IF (DEBUG) THEN
             test=evaluate_overlap('T')
             WRITE(*,'("Test overlaps:",2(I3,X),F8.5)')  n1,n2,test
          end IF
          
          DO g1=0, gmodes
             avec(n1,n2,g1) = evaluate_overlap('V')
             avec(n2,n1,g1) = avec(n1,n2,g1)
          end DO

          DO g1=0, gmodes
             DO g2=0, n2
                amat(n1,n2,g1,g2) = evaluate_overlap('M')
                amat(n2,n1,g1,g2) = amat(n1,n2,g1,g2)
                amat(n1,n2,g2,g1) = amat(n1,n2,g1,g2)
                amat(n2,n1,g2,g1) = amat(n1,n2,g1,g2)
             end DO
             
          end DO

       end DO
    end DO
    
    
  END SUBROUTINE calculate_basis_overlaps

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
    WRITE(OUT,*) a0vec
    WRITE(OUT,*) avec
    WRITE(OUT,*) amat
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

    READ(IN,*) a0vec
    READ(IN,*) avec
    READ(IN,*) amat
    CLOSE(IN)
    
  END SUBROUTINE read_basis_overlaps


END MODULE gh_overlaps
