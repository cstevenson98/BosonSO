!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Summer 2014 - C.Stevenson, St. Andrews 110012872
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This module defines routines which handle the 
!encoding and decoding of the matrix and scalar 
!pair of the problem so that the NAG routine can 
!deal with them. It also holds the definition of 
!the differential equation describing the 
!atom-light system within the cavity.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module state_control
    use params
    use utility, only: ii
    implicit none

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This routine takes the time and the current
!encoded state and returns the encoded state time
!derivative.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diffeqn(t, s, ds_dt)
    implicit none
    
    real(kind=kind(1.0D0)), intent(in)          :: t, s(*) 
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
    
    dlam = (OmegaTil - ii*(Kappa/E_0)) * lam
    
    do i = -nmax, nmax  !Sum over
    do j = -mmax, mmax  !entire mat.
        do k = -1, 1, 2
            if(abs(i+2*k).gt.nmax) then
                !do nothing, ie. add zero to dlam
            else
                dlam = dlam - lam*conjg(mat(i+2*k,j))*mat(i,j) 
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
                dlam = dlam - Pump*conjg(mat(i+k,j+l))*mat(i,j)
            end if
        end do
        end do
    end do
    end do
    
    dlam = -ii*E_0*dlam

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
            dmat(i,j) = dmat(i,j) - (Pump**2)*mat(i,j+2*k)
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
                
        dmat(i,j) = -ii*Omega_r*dmat(i,j)
        
    end do
    end do
    
    !Encode this into the output state.
    call encode_state(dmat, dlam, ds_dt)

end subroutine

end module
