module state_control
!This module defines routines which handle the encoding
!and decoding of the matrix and scalar pair of the problem
!so that the NAG routine can deal with them. It also holds
!the definition of the differential equation describing
!the atom-light system within the cavity.
    use params
    use utility, only: ii
    implicit none

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Flattens matrix into a state vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine encode_state(mat, lam, encvec)
    use params
    implicit none
    
    complex(kind=kind(1.0D0)), intent(in) :: mat(-nmax:nmax,-mmax:mmax), lam
    real(kind=kind(1.0D0)), intent(out)   :: encvec(2*varsize)
    integer                               :: siz
    
    siz = (2*nmax+1)*(2*mmax+1)
    
    encvec(1    :  siz) = reshape(real(mat), (/ siz /))
    encvec(siz+1:2*siz) = reshape(aimag(mat),(/ siz /))
    encvec(2*siz+1:2*varsize) = (/ real(lam), aimag(lam) /)
    
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Extracts scalar and matrix from state vector
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine decode_state(encvec, mat, lam)
    use params
    implicit none
    
    real(kind=kind(1.0D0)), intent(in)     :: encvec(2*varsize)
    complex(kind=kind(1.0D0)), intent(out) :: mat(-nmax:nmax, -mmax:mmax), lam
    integer                                :: siz
    
    siz = (2*nmax+1)*(2*mmax+1)
    
    mat =    reshape(encvec(1    :  siz), (/ 2*nmax+1, 2*mmax+1 /)) + &
          ii*reshape(encvec(siz+1:2*siz), (/ 2*nmax+1, 2*mmax+1 /))
          
    lam = encvec(2*siz+1) + ii*encvec(2*varsize)
    
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Writes matrix to file for easy reading.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_state(state)
    use params
    implicit none
    
    real(kind=kind(1.0D0))    :: state(2*varsize)
    complex(kind=kind(1.0D0)) :: matrix(-nmax:nmax, -mmax:mmax),&
                                 lam
    integer                   :: i, j
    logical                   :: extant
    
    call decode_state(state, matrix, lam)
    
    inquire(file='write_mat.txt', exist=extant)

    if(extant) then
        open(1, file='write_mat.txt', status='old', position='append', action='write')
    else
        open(1,file='write_mat.txt',status='new', action='write')
    end if
    
    do i = -nmax, nmax
        write(1,'(100(2(F6.4,X),3X))') (matrix(i,j), j = -mmax, mmax)
        write(1,'(100(a))') ('-----', j = -2*mmax, 2*mmax)
    end do
    
    write(1,'(2(F5.2,X))') lam
    write(1,*) '-----'
    write(1,*) '#####'
        
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This routine takes the time and the current
!encoded state and returns the encoded state time
!derivative.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diffeqn(t, s, ds_dt)
    use params
    implicit none
    
    real(kind=kind(1.0D0)), intent(in)          :: t 
    real(kind=kind(1.0D0)), intent(in)          :: s(*) 
    real(kind=kind(1.0D0)), intent(out)         :: ds_dt(*)
    real(kind=kind(1.0D0))                      :: dummy
    
    complex(kind=kind(1.0D0)) :: mat(-nmax:nmax, -mmax:mmax)
    complex(kind=kind(1.0D0)) :: lam
    complex(kind=kind(1.0D0)) :: dmat(-nmax:nmax, -mmax:mmax)
    complex(kind=kind(1.0D0)) :: dlam
    
    integer                   :: i,j,k,l

    dummy = t

    !First, get back matrix and scalar.
    call decode_state(s, mat, lam)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Firstly, calculate the rate of change of lam:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    dlam = (OmegaTil - ii*(Kappa/UN)) * lam
    
    do i = -nmax, nmax  !Sum over
    do j = -mmax, mmax  !entire mat.
        do k = -1, 1, 2
            if(abs(i+2*k).gt.nmax) then
                !do nothing, ie. add zero to dlam
            else
                dlam = dlam + lam*conjg(mat(i+2*k,j))*mat(i,j) 
            end if
        end do
    end do
    end do

    do i = -nmax, nmax !Sum over
    do j = -mmax, mmax !entire mat.
        do k = -1, 1, 2
        do l = -1, 1, 2
            if(abs(i+k).gt.nmax .or. abs(j+l).gt.mmax) then
                !do nothing, ie. add zero to dlam
            else
                dlam = dlam + Pump*conjg(mat(i+k,j+l))*mat(i,j)
            end if
        end do
        end do
    end do
    end do
    
    dlam = -ii *(UN/Omega_r)*dlam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Now calculate the elements of the derivative of
!the momentum matrix, phi(i,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    do i = -nmax, nmax !Do this for
    do j = -mmax, mmax !entire mat.
            
        dmat(i,j) = (i**2+j**2)*mat(i,j)
                
        do k = -1, 1, 2
        if(abs(j+2*k).gt.mmax) then
            !do nothing
        else
            dmat(i,j) = dmat(i,j) - Pump*mat(i,j+2*k)
        end if
        end do
                
        do k = -1, 1, 2
        if(abs(i+2*k).gt.nmax) then
            !do nothing
        else
            dmat(i,j) = dmat(i,j) - (abs(lam)**2)*mat(i+2*k,j)
        end if
        end do
                
        do k = -1, 1, 2
        do l = -1, 1, 2
        if(abs(i+k).gt.nmax .or. abs(j+l).gt.mmax) then
            !do nothing
        else
            dmat(i,j) = dmat(i,j) - Pump*(lam+conjg(lam))*mat(i+k,j+l)
        end if
        end do
        end do
                
        dmat(i,j) = -ii*dmat(i,j)
        
    end do
    end do
    
    !Encode this into the output state.
    call encode_state(dmat, dlam, ds_dt)

end subroutine

end module
