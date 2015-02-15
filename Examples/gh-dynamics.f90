MODULE setup_excitation_gh
  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! This module encodes the initial conditions
  
  ! Parameters of excitation profile:
  REAL (KIND=DP) :: x0, width, f0
  
  ! Module level label of which excitation mode we are resolving onto.
  INTEGER, PRIVATE :: g

CONTAINS
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   excitation
  ! VARIABLES:
  !           x - Position
  ! SYNOPSIS:
  ! Function to integrate to give fractional excitation at each
  ! point in space, so as to get projection onto basis states.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION excitation(x)
    REAL (KIND=DP), INTENT(IN) ::x
    REAL (KIND=DP) ::excitation

    excitation = f0 * exp( - 0.5D0*( (x-x0)/width)**2)
    
  END FUNCTION excitation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   integrand
  ! VARIABLES:
  !           x - Position
  ! SYNOPSIS:
  ! Overlap of excitation profile with basis states
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION integrand(x)
    USE hermite, ONLY: wf

    REAL (KIND=DP), INTENT(IN) ::x
    REAL (KIND=DP) ::integrand

    integrand = excitation(x) * wf(x,g)
    
  END FUNCTION integrand

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_initial_excitation
  ! VARIABLES:
  !           uvec - Excitation, resolved onto normal modes, to write
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_initial_excitation(uvec)
    USE gh_overlaps, ONLY: gmodes
    USE utility, ONLY: evaluate_integral

    REAL (KIND=DP), INTENT(OUT) ::uvec(0:gmodes)
    
    ! Read in parameters for excitation profile
    WRITE(*,*)  "Initial excitaiton: Center location; width, amplitude"
    READ(*,*) x0, width, f0

    ! Calculate resolution onto basis states.  n is module level
    DO g=0, gmodes
       uvec(g) = evaluate_integral(integrand)
    end DO
    
    
  END SUBROUTINE setup_initial_excitation



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   compare_excitation
  ! VARIABLES:
  !           rmax - Max radius
  !           nr   - Number of points
  !           uvec - Current exctiation
  !           fname - To write to
  ! SYNOPSIS:
  ! Output initial excitation, and resolved version
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE compare_excitation(rmax,nr,uvec,fname)
    USE hermite, ONLY: vector_density
    USE utility, ONLY: interpolate

    REAL (KIND=DP), INTENT(IN) ::rmax, uvec(:)
    INTEGER, INTENT(IN) ::nr
    CHARACTER (LEN=*) :: fname

    INTEGER :: ir, OUT=40
    REAL (KIND=DP) :: r

    OPEN(OUT,FILE=fname)

    DO ir=1, nr
       r = INTERPOLATE(-rmax, rmax, ir, nr, .FALSE.)
       WRITE(OUT,*)  r, vector_density(uvec,r), excitation(r)
    end DO
    CLOSE(OUT)

  END SUBROUTINE compare_excitation

END MODULE setup_excitation_gh


PROGRAM gh_dynamics
  USE utility, ONLY: interpolate
  USE params, ONLY: setup_params, nmodes, ngl0
  USE gh_overlaps, ONLY: setup_basis_overlaps, statesize, gmodes
  USE hermite, ONLY: vector_density, matrix_density
  USE setup_excitation_gh, ONLY: setup_initial_excitation, compare_excitation
  USE gh_eom, ONLY: setup_eom, diffeqn, encode_state, decode_state, dump_state
  USE run_diffeqn_a, ONLY: setup_diffeqn, step_diffeqn


  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=KIND(1.0D0)

  REAL (KIND=DP) :: t, tmax, r, rmax
  REAL (KIND=DP), ALLOCATABLE ::  uvec(:), s(:), ds(:)
  COMPLEX (KIND=DP), ALLOCATABLE ::  nmat(:,:)
  
  INTEGER :: it, nt, ir, nr
  INTEGER :: RHO_FILE=37, EXC_FILE=38

  LOGICAL :: PROGRESS=.TRUE.
  CHARACTER :: PMARKER*5=' ->  '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Programme to time evolve the photon density matrix and
  ! excition fraction from some initial conditions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set up, or read in basic properties and matrices
  CALL setup_params()


  CALL setup_basis_overlaps()
  IF (PROGRESS) WRITE(*,*) PMARKER,"Calculated overlaps"
  CALL setup_eom()
  IF (PROGRESS) WRITE(*,*) PMARKER,"Done EOM initialisation"

  WRITE(*,*) "Time range to run simulation: tmax, nsteps"
  READ(*,*) tmax, nt
  WRITE(*,*) "Range of positions to print; number of points"
  READ(*,*) rmax, nr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set up the initial conditions
  ALLOCATE(uvec(0:gmodes), nmat(0:nmodes,0:nmodes), &
       & s(statesize), ds(statesize))

  ! No photons
  nmat=0.0D0
  ! Excitation profile according to parameters
  CALL setup_initial_excitation(uvec)

  CALL compare_excitation(rmax, nr, uvec, "initial_excitation.dat")

  ! Convert state into encoded form
  CALL encode_state(nmat, uvec, s)

  ! Reset files
  OPEN(RHO_FILE, FILE="dump-state.dat"); CLOSE(RHO_FILE)
  OPEN(RHO_FILE, FILE="dump-diff.dat"); CLOSE(RHO_FILE)

  CALL dump_state(s, "dump-state.dat")

  ! Set up the differential equation solver
  CALL setup_diffeqn(s, tmax)
  IF (PROGRESS) WRITE(*,*) PMARKER,"Done EOM setup"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Iterate over timesteps, and time evolve equations
  OPEN(RHO_FILE,FILE="density.dat")
  OPEN(EXC_FILE,FILE="excitation.dat")
  DO it=2, nt

     ! Find the next state of the equation
     t = INTERPOLATE(0.0D0, tmax, it, nt, .FALSE.)
     CALL step_diffeqn(s,t,diffeqn)

     ! Debugging: dump the state and the differential of the state:
     CALL dump_state(s, "dump-state.dat")
     CALL diffeqn(0.0D0, s, ds)
     CALL dump_state(ds, "dump-diff.dat")


     IF (PROGRESS) WRITE(*,*) PMARKER, "At step", it

     ! Convert state into matrix and vector
     CALL decode_state(s, nmat, uvec)

     ! Write out the corresponding photon density profile.
     DO ir=1, nr
        r = INTERPOLATE(-rmax, rmax, ir, nr, .FALSE.)
        WRITE(RHO_FILE,*) t, r, matrix_density(nmat,r)*ngl0
        WRITE(EXC_FILE,*) t, r, vector_density(uvec,r)
     end DO
     WRITE(RHO_FILE,*)
     WRITE(EXC_FILE,*)

  end DO
  CLOSE(RHO_FILE)
  CLOSE(EXC_FILE)

END PROGRAM gh_dynamics

