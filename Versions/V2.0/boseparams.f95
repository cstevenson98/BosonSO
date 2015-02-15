!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Summer 2014 - C.Stevenson, St. Andrews 110012872
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This module contains all of the physical and 
!mathematical constants related to the problem.
!These are then able to be used globally.
!
!   nmax, mmax - x, y extent of momentum-space
!                matrix.        
!
!   varsize    - the number of complex variables
!                in problem.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
module params
    use utility
    implicit none
    
    !Flags for desired program operation.
    integer                :: Multiple, StateKeep, &
                              PumpEvolve, Up

    !Size of matrix, and total no. of complex variables
    integer                :: nmax, mmax, varsize
    
    !Constants of the problem.
    real(kind=kind(1.0D0)) :: E_0, Omega_0, Kappa, N_A

    !Variables of the problem, in the form of arrays
    !containing start and end values as well as step.
    real(kind=kind(1.0D0)) :: A_Omega(3), & !(/OmInit, OmEnd, OmStep/)
                              A_gsqrtN(3)!(/gNInit, gNEnd, gNStep/)

    !The variables which hold the current value within the
    !range of the above arrays.
    real(kind=kind(1.0D0)) :: Omega, gsqrtN, & 
                              Omega_r, OmegaTil, Pump

    !The state of interest
    real(kind=kind(1.0D0)), allocatable :: CurrentState(:)

contains

subroutine setup()
    implicit none

    logical :: extant
  
    read(*,*) Multiple, StateKeep, PumpEvolve, Up
    read(*,*) nmax, E_0, Omega_0, Kappa, N_A
    
    mmax = nmax
    varsize = (2*nmax+1)*(2*mmax+1)+1
    allocate(CurrentState(2*varsize))

    read(*,*) A_Omega, A_gsqrtN

    !A_Omega_r = 0.5*Omega_0 + (A_gsqrtN**2)/(8*E_0)
    !A_Pump = 0.5*A_gsqrtN/sqrt((A_Omega_r*abs(E_0)))
    !A_OmegaTil = A_Omega/E_0 + 0.5

    call read_continue()
   
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine checks if there is a continue 
!file/flag. If there is, this routine uses it to
!instantiate all of the consants else it set them
!equal to the first element of the arrays, uses
!initial condition for the 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_continue()

    real(kind=kind(1.0D0)) :: state(2*varsize)
    logical :: extant

    inquire(file='continue.dat', exist=extant)
    if(extant) then 
        open(23, file='continue.dat', status='old', action='read')
        read(23,*) gsqrtN, Omega
        if(StateKeep .eq. 1) then
        read(23,*) CurrentState
        else 
        call initial(CurrentState)
        end if
        close(23)
    else
        gsqrtN = A_gsqrtN(1)
        Omega = A_Omega(1)
        call initial(CurrentState)
    end if

    Omega_r =  0.5*Omega_0 + (gsqrtN**2)/(8*E_0)
    Pump = 0.5*gsqrtN/sqrt((Omega_r*abs(E_0)))
    OmegaTil = Omega/E_0 + 0.5
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Makes continue file/flag by writing current state
!and the pump and frequency.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_continue() !###########################################
    implicit none
    
    real(kind=kind(1.0D0)) :: state(2*varsize)
    logical :: extant

    if(Multiple .eq. 0) then
        return
    end if

    inquire(file='continue.dat', exist=extant)

    if(extant) then
        open(20, file='continue.dat', status='old', action='write')
    else
        open(20,file='continue.dat',status='new', action='write')
    end if
    
    gsqrtN = gsqrtN + (A_gsqrtN(2) - A_gsqrtN(1))/A_gsqrtN(3)

    if(gsqrtN .ge. A_gsqrtN(2)) then
        write(*,*) '-----------------------------------'
        Omega = Omega + (A_Omega(2)-A_Omega(1))/A_Omega(3)

        if(Omega .ge. A_Omega(2)) then
            !close(20, status='delete')
            return

        else
            gsqrtN = A_gsqrtN(1)
            write(20,*) gsqrtN, Omega
            write(20,*) CurrentState
            close(20)
            return

        end if
    else 
        write(20,*) gsqrtN, Omega
        write(20,*) CurrentState
        close(20)
        return
    end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Flattens matrix into a state vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine encode_state(mat, lam, encvec)
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
    implicit none
    
    real(kind=kind(1.0D0)), intent(in)     :: encvec(2*varsize)
    complex(kind=kind(1.0D0)), intent(out) :: mat(-nmax:nmax, -mmax:mmax), lam
    integer                                :: siz
    
    siz = (2*nmax+1)*(2*mmax+1)
    
    mat =    reshape(encvec(1    :  siz), (/ 2*nmax+1, 2*mmax+1 /)) + &
          ii*reshape(encvec(siz+1:2*siz), (/ 2*nmax+1, 2*mmax+1 /))
          
    lam = encvec(2*siz+1) + ii*encvec(2*varsize)
    
end subroutine

subroutine initial(state)
    implicit none

    real(kind=kind(1.0D0)) :: state(2*varsize)
    complex(kind=kind(1.0D0)) :: Phi(-nmax:nmax, -mmax:mmax), Lam

    Phi = (1-ii)/(2*Sqrt(2*N_A)); Phi(0,0) = 1.0
    Phi(1,0)=0.0; Phi(0,1)=0.0; Phi(-1,0)=0.0; Phi(0,-1)=0.0
    Lam = 0.0

    call encode_state(Phi, Lam, state)
end subroutine

end module
