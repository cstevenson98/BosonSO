MODULE setup_excitation_rs
  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! Encodes initial conditions
  
CONTAINS


  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_initial_excitation
  ! VARIABLES:
  !           fvec - Excitation, resolved onto normal modes, to write
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_initial_excitation(fvec)
    USE rs_basis, ONLY: gmodes, dx, x0
    USE utility, ONLY: evaluate_integral


    REAL (KIND=DP), INTENT(OUT) ::fvec(gmodes)
    REAL (KIND=DP) :: x1, width, f0
    INTEGER :: g

    ! Read in parameters for excitation profile
    WRITE(*,*)  "Initial excitaiton: Center location; width, amplitude"
    READ(*,*) x1, width, f0

    DO g=1, gmodes
       fvec(g) = excitation(x0 + (g-1)*dx)
    end DO



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

    excitation = f0 * exp( - 0.5D0*( (x-x1)/width)**2)
    
  END FUNCTION excitation


  end SUBROUTINE setup_initial_excitation
  
END MODULE setup_excitation_rs

PROGRAM rs_dynamics
  USE utility, ONLY: interpolate, assert
  USE params, ONLY: setup_params, nmodes, ngl0
  USE rs_basis, ONLY: setup_basis_overlaps, statesize, gmodes, x0, dx
  USE hermite, ONLY: matrix_density
  USE setup_excitation_rs, ONLY: setup_initial_excitation
  USE rs_eom, ONLY: setup_eom, diffeqn, encode_state, decode_state, dump_state
  USE run_diffeqn_b, ONLY: setup_diffeqn, step_diffeqn

    USE params, ONLY: nmodes
    USE hermite, ONLY: wf


  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=KIND(1.0D0)

  REAL (KIND=DP) :: t, tmax, r, rmax
  REAL (KIND=DP), ALLOCATABLE ::  fvec(:), s(:), ds(:)
  COMPLEX (KIND=DP), ALLOCATABLE ::  nmat(:,:)
  
  INTEGER :: it, nt, ir, nr, n
  INTEGER :: RHO_FILE=37, EXC_FILE=38

  LOGICAL :: PROGRESS=.TRUE.
  CHARACTER :: PMARKER*5=' ->  '

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Programme to time evolve the photon density matrix and
  ! excition fraction from some initial conditions.  Working in
  ! real space.
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


  
  OPEN(RHO_FILE, FILE="wfs.dat")
  DO n=0, nmodes
     DO ir=1, gmodes
        r = x0 + dx*(ir-1)
        WRITE(RHO_FILE,*)  r, wf(r,n), n
     end DO
     WRITE(RHO_FILE,*)
  end DO
  CLOSE(RHO_FILE)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set up the initial conditions
  ALLOCATE(fvec(gmodes), nmat(0:nmodes,0:nmodes), s(statesize), ds(statesize))

  
  ! No photons
  nmat=0.0D0
  ! Excitation profile according to parameters
  CALL setup_initial_excitation(fvec)


  ! Convert state into encoded form
  CALL encode_state(nmat, fvec, s)

  ! Set up the differential equation solver
  CALL setup_diffeqn(s, tmax)
  IF (PROGRESS) WRITE(*,*) PMARKER,"Done EOM setup"

  OPEN(RHO_FILE, STATUS="REPLACE", FILE="dump-state-rs.dat"); CLOSE(RHO_FILE)
  OPEN(RHO_FILE, STATUS="REPLACE", FILE="dump-diff-rs.dat"); CLOSE(RHO_FILE)

  CALL dump_state(s, "dump-state-rs.dat")
  CALL diffeqn(0.0D0, s, ds)
  CALL dump_state(ds, "dump-diff-rs.dat")

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Iterate over timesteps, and time evolve equations
  OPEN(RHO_FILE,FILE="density.dat")
  OPEN(EXC_FILE,FILE="excitation.dat")
  DO it=1, nt

     ! Find the next state of the equation
     t = INTERPOLATE(0.0D0, tmax, it, nt, .FALSE.)
     IF (t .GT. 0.0D0) CALL step_diffeqn(s,t,diffeqn)

     IF (PROGRESS) WRITE(*,*) PMARKER, "At step", it


     CALL diffeqn(0.0D0, s, ds)
     CALL dump_state(s, "dump-state-rs.dat")
     CALL dump_state(ds, "dump-diff-rs.dat")

     ! Convert state into matrix and vector parts:
     CALL decode_state(s, nmat, fvec)

     ! Write out the corresponding photon density profile.
     DO ir=1, nr
        r = INTERPOLATE(-rmax, rmax, ir, nr, .FALSE.)
        WRITE(RHO_FILE,*) t, r, matrix_density(nmat,r)*ngl0
     end DO
     WRITE(RHO_FILE,*)

     DO ir=1, gmodes
        r = x0 + dx*(ir-1)
        WRITE(EXC_FILE,*) t, r, fvec(ir)
     end DO
     WRITE(EXC_FILE,*)

     CALL ASSERT(COUNT(fvec.LT.0.0D0) .EQ. 0, "Negative excitation exists")

  end DO
  CLOSE(RHO_FILE)
  CLOSE(EXC_FILE)


  
END PROGRAM rs_dynamics
