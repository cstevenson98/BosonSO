MODULE gh_eom
  USE gh_overlaps, ONLY: statesize, gmodes
  USE params, ONLY: nmodes

  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! This module contains actual differential equation to be evolved,
  ! working in the basis where gain is projected onto Gauss Hermite
  ! functions

  ! Storage for helper matrices, and temporaray results.  aup, adn complex
  ! to allow using matrix multiplication in line with that for nmat.
  COMPLEX (KIND=DP), ALLOCATABLE :: aup(:,:), adn(:,:)
  REAL (KIND=DP), ALLOCATABLE :: h0(:), gtmp(:,:)
  REAL (KIND=DP), ALLOCATABLE :: amatpm(:,:,:,:)
  REAL (KIND=DP), ALLOCATABLE :: aspnt(:,:), avecp(:,:,:)


  
CONTAINS
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_eom
  ! VARIABLES:
  ! SYNOPSIS:
  ! Set up temp storage etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_eom()
    USE params, ONLY: omega, gm, gp
    USE gh_overlaps, ONLY: amat, avec

    INTEGER :: n, g1, g2

    ! Set up temporary storage
    ALLOCATE(aup(0:nmodes, 0:nmodes), adn(0:nmodes, 0:nmodes))
    ALLOCATE(gtmp(0:gmodes, 0:gmodes))

    ! Set up helper quantities, such as product of A matrix
    ! with G matrices, and spontaneous emission matrix (which
    ! is independent of photon state)
    ALLOCATE(amatpm(0:nmodes, 0:nmodes, 0:gmodes, 0:gmodes), &
         &    avecp(0:nmodes, 0:nmodes, 0:gmodes), &
         &    aspnt(0:gmodes, 0:gmodes), h0(0:nmodes))

    h0 = (/ (n*omega, n=0, nmodes) /)

    DO n=0, nmodes
       ! Amatpm appears in the evoluation of u; it is the thing
       ! which has its trace taken with n, and then gives a matrix
       ! into u space
       amatpm(:, n, :, :) = amat(:, n, :, :) * (gp(n) + gm(n))

       ! avecp also appears in u, but gives a vector in u space after
       ! the trace with n.  This corresponds to the rate of absorption,
       ! i.e. the part that doesn't depend on existing excited state
       ! population.
       avecp(:, n, :) = avec(:,n,:) * gp(n)

    end DO

    ! Matrix of spontaneous emission in excitation basis.
    DO g1=0, gmodes
       DO g2=0, gmodes
           aspnt(g1,g2) = SUM((/ ( amat(n,n,g1,g2)*gm(n), n=0, nmodes) /))
       end DO
    end DO
    
          
  END SUBROUTINE setup_eom


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   encode_state
  ! VARIABLES:
  !           nmat - Number matrix
  !           uvec - Gain vector
  !           s    - Flattened structure
  ! SYNOPSIS:
  ! Put things in flattened structure
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE encode_state(nmat,uvec,s)
    COMPLEX (KIND=DP), INTENT(IN) :: nmat(0:nmodes,0:nmodes)
    REAL (KIND=DP), INTENT(IN) :: uvec(0:gmodes)
    REAL (KIND=DP), INTENT(OUT) ::s(statesize)
    
    INTEGER :: sz

    sz=(nmodes+1)*(nmodes+1)

    s(1   :  sz) = RESHAPE(REAL(nmat),  (/ sz /))
    s(1+sz:2*sz) = RESHAPE(AIMAG(nmat), (/ sz /))
    s(2*sz+1:statesize) = uvec
    

  END SUBROUTINE encode_state


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   decode_state
  ! VARIABLES:
  !           s    - Flattened structure
  !           nmat - Number matrix
  !           uvec - Gain vector
  ! SYNOPSIS:
  ! Take things out of flattened structure.  FIXME: Should note hermiticity
  ! to halve the size of state.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE decode_state(s,nmat,uvec)
    USE utility, ONLY: ii
    REAL (KIND=DP), INTENT(IN) ::s(statesize)
    COMPLEX (KIND=DP), INTENT(OUT) :: nmat(0:nmodes,0:nmodes)
    REAL (KIND=DP), INTENT(OUT) :: uvec(0:gmodes)

    INTEGER :: sz

    sz=(nmodes+1)*(nmodes+1)

    nmat =    RESHAPE(s(1   :  sz), (/ nmodes+1, nmodes+1 /) ) + &
         & ii*RESHAPE(s(1+sz:2*sz), (/ nmodes+1, nmodes+1 /) )

    uvec = s(2*sz+1:statesize) 

  END SUBROUTINE decode_state



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   dump_state
  ! VARIABLES:
  !           s     - State
  !           fname - Filename to write
  ! SYNOPSIS:
  ! Dump state in human readable format to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE dump_state(s,fname)
    REAL (KIND=DP), INTENT(IN) ::s(:)
    CHARACTER(LEN=*), INTENT(IN) ::fname

    INTEGER :: OUT=47, n
    REAL    (KIND=DP), DIMENSION(0:gmodes)          :: uvec
    COMPLEX (KIND=DP), DIMENSION(0:nmodes,0:nmodes) :: nmat

    CALL decode_state(s,nmat,uvec)

    OPEN(OUT,FILE=fname,ACCESS='APPEND')
    WRITE(OUT,'("# Modes: Photon:",I3," Excitation: ",I3)') nmodes, gmodes
    WRITE(OUT,'("# Excitation vector")') 
    DO n=0, gmodes
       WRITE(OUT,'(D12.5)')  uvec(n)
    end DO
    

    WRITE(OUT,*) 
    WRITE(OUT,'("# Photon matrix")') 
    DO n=0, nmodes
       WRITE(OUT,'(50(2(F9.3,X),4X))')  nmat(n,:)
    end DO
    
    WRITE(OUT,*)
    WRITE(OUT,'("# HC check: ",D12.5)') SUM(ABS(nmat-CONJG(TRANSPOSE(nmat))))

    CLOSE(OUT)

  END SUBROUTINE dump_state


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   diffeqn
  ! VARIABLES:
  !           t   - Current time
  !           s   - Current state
  !           ds  - change of state
  ! SYNOPSIS:
  ! Differential equation to solve
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE diffeqn(t,s,ds)
    USE utility, ONLY: ii, hermitian, dmatmul
    USE params
    USE gh_overlaps, ONLY: statesize, a0vec, avec, amat

    REAL (KIND=DP), INTENT(IN) ::t,s(statesize)
    REAL (KIND=DP), INTENT(OUT) :: ds(statesize)
    REAL (KIND=DP) :: dummy


    REAL    (KIND=DP), DIMENSION(0:gmodes)          :: uvec, duvec
    COMPLEX (KIND=DP), DIMENSION(0:nmodes,0:nmodes) :: nmat, dnmat

    INTEGER :: n, m,  g1, g2

    ! To avoid unused variable error message
    dummy=t

    CALL decode_state(s,nmat,uvec)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Create a up, and adn = 1 - aup
    aup=0.0D0
    DO n=0, nmodes
       DO m=0, nmodes
          aup(n,m) = SUM(avec(n,m,:)*uvec)
       end DO
    end DO

    adn=0.0D0
    DO n=0, nmodes
       adn(n,n) = 1.0D0
    end DO
    adn=adn-aup
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Time evolution of nmat

    ! Commutator with H
    dnmat = 2*ii * DMATMUL(nmat, h0, 'R') - kappa * nmat 

    ! Anticommutator with G-, aup -- spontaneous emission.
    dnmat = dnmat + DMATMUL(aup, gm, 'R')/ngl0

    ! Other gain parts
    dnmat = dnmat + MATMUL( nmat, DMATMUL(aup, gm, 'R'))
    
    ! Loss parts
    dnmat = dnmat - MATMUL(DMATMUL(nmat, gp, 'R'), adn)

    dnmat=HERMITIAN(dnmat)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Time evolution of uvec
    duvec = - GammaDown*uvec + GammaUp*(a0vec - uvec)
    
    ! Spontaneous emission part 
    duvec = duvec - MATMUL(aspnt, uvec) / ngl0

    ! Stimulated emission and suppression of absorption.  Note that
    ! rates are all real, due to sum of trace and trace of transpose
    ! (which is complex conjugate by hermiticity)
    ! Note that this has both the gp and gm parts, viva amatpm
    DO g1=0, gmodes
       DO g2=0, gmodes
          ! Taking Trace (A.B) = sum_n ( sum_m A_{mn} B_{nm})
          gtmp(g1,g2) = REAL( &
               & SUM((/ ( SUM(nmat(n,:)*amatpm(:,n,g1,g2)), n=0, nmodes) /)) )
       end DO
    end DO

    duvec = duvec - MATMUL(gtmp, uvec) 
    

    ! Absorption parts; thes corresponding saturation of this 
    ! is in the definition of amatpm.
    ! Vector in u space first.  Note the abuse of gtmp.

    DO g1=0, gmodes
       gtmp(g1,0) = REAL(&
            & SUM((/ ( SUM(nmat(n,:)*avecp(:,n,g1)), n=0, nmodes) /)) )
    end DO
    duvec = duvec + gtmp(:,0) 


    CALL encode_state(dnmat, duvec, ds)

  END SUBROUTINE diffeqn

END MODULE gh_eom
