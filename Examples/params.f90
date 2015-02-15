MODULE params
  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! Keeps track of the basic parameters in the model
  ! TLS loss and pumping; photon loss
  REAL (KIND=DP) :: GammaDown, GammaUp, kappa
  ! Number of gain molecules in one harmonic oscillator length
  REAL (KIND=DP) :: ngl0
  ! Harmonic oscillator frequency
  REAL (KIND=DP) :: omega
  
  ! Number of modes to consider
  INTEGER :: nmodes
  ! TLS gain profile, with prefactor of gain molecule density, i.e.
  ! corresponding to old N*Gamma(\pm delta-m).
  REAL (KIND=DP), ALLOCATABLE :: Gp(:), Gm(:)

CONTAINS


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_params
  ! VARIABLES:
  ! SYNOPSIS:
  !  Read in parameters from file.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_params()
    CHARACTER :: calcGstr

    WRITE(*,*) "Loss parametesr: Gamma_down, Gamma_up, kappa"
    READ(*,*) GammaDown, GammaUp, kappa
    WRITE(*,*) "Gain density (HO units), HO frequency"
    READ(*,*) ngl0, omega

    WRITE(*,*) "Max mode number; how to calculate G [D/F]"
    READ(*,*) nmodes, calcGstr


    ALLOCATE(GP(0:nmodes), GM(0:nmodes))

    SELECT CASE(calcGSTR)
    CASE('F')
       WRITE(*, *) "Reading gain profie from file"
       CALL read_gain_profile()

    CASE DEFAULT
       WRITE(*,*) "Using dummy model to allocate distribution"
       CALL setup_dummy_gain_profile()
       CALL write_gain_profile("gain-profiles.dat")
    end SELECT

    
  END SUBROUTINE setup_params
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   read_gain_profile
  ! SYNOPSIS:
  ! Read in gain profile from file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_gain_profile()
    USE utility, ONLY: assert

    INTEGER :: IN=27, file_nmodes, n
    CHARACTER :: fname*128
    
    WRITE(*,*) "Filename to read from"
    READ(*,*) fname
    
    OPEN(IN, FILE=fname)
    READ(IN,*) file_nmodes
    CALL assert(file_nmodes.EQ.nmodes, "Incorrect number of modes in file")
    DO n=0, nmodes
       READ(IN,*) GP(n), GM(n)
    end DO
    
    CLOSE(IN)
    
  END SUBROUTINE read_gain_profile

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   write_gain_profile
  ! SYNOPSIS:
  ! Write gain profile to file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE write_gain_profile(fname)
    USE utility, ONLY: assert
    CHARACTER(LEN=*) :: fname

    INTEGER :: OUT=27, file_nmodes, n
    
    OPEN(OUT, FILE=fname)
    WRITE(OUT,*) nmodes
    DO n=0, nmodes
       WRITE(OUT,*) GP(n), GM(n)
    end DO
    
    CLOSE(OUT)
    
  END SUBROUTINE write_gain_profile


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_dummy_gain_profile
  ! SYNOPSIS:
  !  Set up gain profiel from a dumb weighted Gaussian profile.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_dummy_gain_profile()


    REAL (KIND=DP) :: T, d0, width, en, NGamma0
    INTEGER :: n

    WRITE(*,*) "Temperature [energy], delta_0, linewidth, N*Gamma0"
    READ(*,*) T, width, d0, NGamma0

    DO n=0, nmodes
       en = d0 + n*omega

       GP(n) = NGamma0*exp( - (en/width)**2  + 0.5*en/T)
       GM(n) = NGamma0*exp( - (en/width)**2  - 0.5*en/T)
       
    end DO
    

  END SUBROUTINE setup_dummy_gain_profile

END MODULE params
