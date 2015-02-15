MODULE rs_eom
  USE rs_basis, ONLY: statesize, gmodes
  USE params, ONLY: nmodes

  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! This module contains the differential equation where gain
  ! is written in real space.

  ! Storage for helper matrices, and temporaray results.  aup is complex
  ! to allow using matrix multiplication in line with that for nmat.
  COMPLEX (KIND=DP), ALLOCATABLE :: aup(:,:), id(:,:)
  REAL (KIND=DP), ALLOCATABLE :: h0(:), avecm(:,:,:), avecp(:,:,:)
  REAL (KIND=DP), ALLOCATABLE :: gtup(:), gtdn(:), aspnt(:)
  
  
CONTAINS
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_eom
  ! VARIABLES:
  ! SYNOPSIS:
  ! Set up temp storage etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_eom()
    USE params, ONLY: omega, gm, gp
    USE rs_basis, ONLY:  avec

    INTEGER :: n, g

    ! Set up temporary storage
    ALLOCATE(aup(0:nmodes, 0:nmodes),Gtup(gmodes), Gtdn(gmodes))

    ! Set up helper functions
    ALLOCATE(id(0:nmodes,0:nmodes), h0(0:nmodes), aspnt(gmodes), &
         &   avecm(0:nmodes, 0:nmodes, gmodes), &
         &   avecp(0:nmodes, 0:nmodes, gmodes))


    h0 = (/ (n*omega, n=0, nmodes) /)

    id=0.0D0
    DO n=0, nmodes
       ! Id matrix
       id(n,n) = 1.0D0

       ! Post multiplication of the a vector by G plus and minus
       avecm(:, n, :) = avec(:, n, :) *  gm(n)
       avecp(:, n, :) = avec(:, n, :) *  gp(n)

    end DO
    
    ! Spontaneous emission rate:  trace of avecm
    DO g=1, gmodes
       aspnt(g) =  REAL( SUM((/(avecm(n,n,g),n=0, nmodes) /)) ) 
    end DO

  end SUBROUTINE setup_eom
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   encode_state
  ! VARIABLES:
  !           nmat - Number matrix
  !           fvec - Gain vector
  !           s    - Flattened structure
  ! SYNOPSIS:
  ! Put things in flattened structure
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE encode_state(nmat,fvec,s)
    COMPLEX (KIND=DP), INTENT(IN) :: nmat(0:nmodes,0:nmodes)
    REAL (KIND=DP), INTENT(IN) :: fvec(gmodes)
    REAL (KIND=DP), INTENT(OUT) ::s(statesize)
    
    INTEGER :: sz
    
    sz=(nmodes+1)*(nmodes+1)

    s(1   :  sz) = RESHAPE(REAL(nmat),  (/ sz /))
    s(1+sz:2*sz) = RESHAPE(AIMAG(nmat), (/ sz /))
    s(2*sz+1:statesize) = fvec
    

  END SUBROUTINE encode_state


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   decode_state
  ! VARIABLES:
  !           s    - Flattened structure
  !           nmat - Number matrix
  !           fvec - Gain vector
  ! SYNOPSIS:
  ! Take things out of flattened structure.  FIXME: Should note hermiticity
  ! to halve the size of state.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE decode_state(s,nmat,fvec)
    USE utility, ONLY: ii
    REAL (KIND=DP), INTENT(IN) ::s(statesize)
    COMPLEX (KIND=DP), INTENT(OUT) :: nmat(0:nmodes,0:nmodes)
    REAL (KIND=DP), INTENT(OUT) :: fvec(gmodes)

    INTEGER :: sz

    sz=(nmodes+1)*(nmodes+1)

    nmat =    RESHAPE(s(1   :  sz), (/ nmodes+1, nmodes+1 /) ) + &
         & ii*RESHAPE(s(1+sz:2*sz), (/ nmodes+1, nmodes+1 /) )

    fvec = s(2*sz+1:statesize) 

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
    REAL    (KIND=DP), DIMENSION(gmodes)            :: fvec
    COMPLEX (KIND=DP), DIMENSION(0:nmodes,0:nmodes) :: nmat

    CALL decode_state(s,nmat,fvec)

    OPEN(OUT,FILE=fname,ACCESS='APPEND')
    WRITE(OUT,'("# Modes: Photon:",I3," Excitation: ",I3)') nmodes, gmodes
    WRITE(OUT,'("# Excitation vector")') 
    DO n=1, gmodes
       WRITE(OUT,'(D12.5)')  fvec(n)
    end DO
    

    WRITE(OUT,*) 
    WRITE(OUT,'("# Photon matrix")') 
    DO n=0, nmodes
       WRITE(OUT,'(50(2(D9.2,X),4X))')  nmat(n,:)
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
  ! Note that we are using the version where n has been rescale,d
  ! so ngl0 appears with the spontaneous rates
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE diffeqn(t,s,ds)
    USE utility, ONLY: ii, hermitian, dmatmul
    USE params
    USE rs_basis, ONLY: avec, dx

    REAL (KIND=DP), INTENT(IN) ::t,s(statesize)
    REAL (KIND=DP), INTENT(OUT) :: ds(statesize)
    REAL (KIND=DP) :: dummy

    REAL    (KIND=DP), DIMENSION(gmodes)            :: fvec, dfvec
    COMPLEX (KIND=DP), DIMENSION(0:nmodes,0:nmodes) :: nmat, dnmat
    
    INTEGER :: n, m, g

    ! To avoid unused variable error message
    dummy=t

    CALL decode_state(s,nmat,fvec)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local helper quantities:
    ! a up and down.

    DO n=0, nmodes
       DO m=0, nmodes
          aup(n,m) = dx * SUM(avec(n,m,:) * fvec)
       end DO
    end DO
    
    ! Tilde Gamma -- local rates.
    DO g=1, gmodes
       gtdn(g) = GammaDown + aspnt(g)/ngl0 &
            &  + REAL(SUM((/(SUM(avecm(n,:,g)*nmat(:,n)), n=0, nmodes ) /)))

       gtup(g) = GammaUp &
            &  + REAL(SUM((/(SUM(avecp(n,:,g)*nmat(:,n)), n=0, nmodes ) /)))
    end DO
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Time evolution of photon number.  Note that we build it in parts
    ! first, and then symmetrise at the end.

    dnmat = 2*ii*DMATMUL(nmat, h0, 'R') - kappa * nmat 

    ! Emission: Spontaneous and stimulated
    dnmat = dnmat + MATMUL(nmat+id/ngl0, DMATMUL(aup, gm, 'R')) 

    ! Absorption
    dnmat = dnmat - MATMUL(DMATMUL(nmat, gp, 'R'), id - aup) 
    
    dnmat=HERMITIAN(dnmat)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Time evolution of state: Easy noe we have Gamma tilde
    
    dfvec = gtup - (gtdn + gtup) * fvec 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL encode_state(dnmat, dfvec, ds)


  end SUBROUTINE diffeqn
  


END MODULE rs_eom
