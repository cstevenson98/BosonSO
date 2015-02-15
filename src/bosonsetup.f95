!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Summer 2014 - C.Stevenson, St. Andrews 110012872
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This module contains all of the physical and 
!mathematical constants related to the problem.
!These are then able to be used globally.
!
!   nmax, mmax - x, y extent of momentum-space
!                matrix.        
!
!   varsize    - the number of complex variables
!                in problem.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module params
    use utility
    implicit none
    
    !Flags for desired program operation.
    integer                :: StateKeep, TestFlag
    character(len=15)      :: RunType

    !Initial time of each run and the maximum possible time we
    !wish to run code for. The time step.
    real(kind=kind(1.0D0))    :: tinit, tend, tstep

    ! Convergence threshold
    real(kind=kind(1.0D0))    :: convergence_threshold

    !Size of matrix, and total no. of complex variables
    integer                :: nmax, mmax, varsize
    
    !Constants of the problem.
    real(kind=kind(1.0D0)) :: E_0, Kappa, N_A, Omega_r

    !Arrays to hold range on parameters pump and omegatil
    !to use, and the number of steps in the format:
    !{init, end, step count}
    real(kind=kind(1.0D0)) :: A_Pump(3), A_OmegaTil(3)

    !The variables which hold the current value within the
    !range of the above arrays.
    real(kind=kind(1.0D0)) :: Pump, OmegaTil

    !The state of interest
    real(kind=kind(1.0D0)), allocatable :: CurrentState(:)
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    
    !The follwing define the universal parameters which
    !contain the names for the output files.
    character(len=100) :: TimeOutput, SweepOutput
    
    !And the local min and max for the convergence criterion.
    real(kind=kind(1.0D0)) :: RunMax, RunMin
    

contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!This is the main routine which sets up the program
!for one run. Reads from file as well as some hard
!coding.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setup()
    implicit none

    logical :: extant
    
    open(12, file='params', status='old', action='read')
  
    read(12,*) StateKeep, TestFlag, RunType
    read(12,*) nmax, E_0, Kappa, N_A, Omega_r
    
    mmax = nmax
    varsize = (2*nmax+1)*(2*mmax+1)+1
    allocate(CurrentState(2*varsize))

    read(12,*) A_Pump, A_OmegaTil
    
    read(12,*) tinit, tend, tstep

    read(12,*) convergence_threshold


    close(12)

    Pump = A_Pump(1); OmegaTil = A_OmegaTil(1)
 

end subroutine

subroutine reset_params()
    
    Pump = A_Pump(1); OmegaTil = A_OmegaTil(1)
    call initial(CurrentState)
    
end subroutine

!function data_cluster(variety, t)

!    !Function output: an array containing data we want to write.
!    Complex(kind=kind(1.0D0)), allocatable :: data_cluster(:)

!    !What kind of data we want in the araay.
!    character(len=*)                   :: variety

!    !The time; only the evolving datd requires this.
!    real(kind=kind(1.0D0))              :: t

!    !We are required to include holder matrix and scalar,
!    !as we wish to output quantities involving these.
!    complex(kind=kind(1.0D0)) :: Phi(-nmax:nmax, -mmax:mmax), Lam

!    call encode_state(Phi, Lam, CurrentState)

!    select case(variety)
!        case('sweep_data')
!            allocate(data_cluster(4))
!            data_cluster = (/ Pump, Abs(Lam)**2, 2*Conjg(Phi(1,1))*Phi(0,0), &
!                           & 0.5*(Abs(2*Phi(1,1))**2 - Abs(Phi(0,0))**2) /)

!        case('evolving_data')
!            allocate(data_cluster(4))
!            data_cluster = (/ t/1000.0, sqrt((N_A*Omega_r/E_0))*Lam, &
!                            & 0.5*(abs(2*Phi(1,1))**2 - &
!                            & abs(Phi(0,0))**2), abs(Phi(1,0))**2 /)

!    end select
!    
!end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!This routine will change PUMP by one 
!increment.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine change_pump(direction)

    character(len=*) :: direction

    select case(direction)
        case('up') 
            Pump = Pump + (A_Pump(2) - A_Pump(1))/A_Pump(3)
        
        case('down')
            Pump = Pump - (A_Pump(2) - A_Pump(1))/A_Pump(3)
    end select

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!This routine will change OMEGATIL by one 
!increment.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine change_omegatil(direction)

    character(len=*) :: direction

    select case(direction)
        case('up') 
            OmegaTil = &
            & OmegaTil + (A_OmegaTil(2) - A_OmegaTil(1))/A_OmegaTil(3)
        
        case('down')
            OmegaTil = &
            & OmegaTil - (A_OmegaTil(2) - A_OmegaTil(1))/A_OmegaTil(3)
        end select 

end subroutine


subroutine write_output(t, variety)
    implicit none
    
    real(kind=kind(1.0D0)) :: t
    character(len=*) :: variety
    complex(kind=kind(1.0D0)) :: Phi(-nmax:nmax, -mmax:mmax), Lam
    
    real (KIND=kind(1.0D0)) :: sz, ipr
    complex (KIND=kind(1.0D0)) :: sp

    call decode_state(CurrentState, Phi, Lam)
    
    ! Spin components.
    sz=0.5*(abs(2*Phi(1,1))**2 - abs(Phi(0,0))**2)
    sp=2*Conjg(Phi(1,1))*Phi(0,0)

    ! Inverse participation ratio
    ipr = SUM(ABS(phi)**4) / ( (SUM(ABS(phi)**2))**2 )

    select case(variety)
        case('time')
           write(51,'(X,6(D14.7,X))') t/1000.0, Lam, sz, &
                            & ipr, abs(Phi(1,0))**2
        case('sweep')
            write(52,'(X,8(D12.5,X))') Pump, Abs(Lam)**2, &
                               & Lam, sp, sz, ipr
    end select
    
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Flattens matrix into a state vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine encode_state(Phi, Lam, state)
    implicit none
    
    complex(kind=kind(1.0D0)), intent(in) :: Phi(-nmax:nmax,-mmax:mmax), Lam
    real(kind=kind(1.0D0)), intent(out)   :: state(2*varsize)
    integer                               :: siz
    
    siz = (2*nmax+1)*(2*mmax+1)
    
    state(1    :  siz) = reshape(real(Phi), (/ siz /))
    state(siz+1:2*siz) = reshape(aimag(Phi),(/ siz /))
    state(2*siz+1:2*varsize) = (/ real(Lam), aimag(Lam) /)
    
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Extracts scalar and matrix from state vector
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine decode_state(state, Phi, Lam)
    implicit none
    
    real(kind=kind(1.0D0)), intent(in)     :: state(2*varsize)
    complex(kind=kind(1.0D0)), intent(out) :: Phi(-nmax:nmax, -mmax:mmax), Lam
    integer                                :: siz
    
    siz = (2*nmax+1)*(2*mmax+1)
    
    Phi =    reshape(state(1    :  siz), (/ 2*nmax+1, 2*mmax+1 /)) + &
          ii*reshape(state(siz+1:2*siz), (/ 2*nmax+1, 2*mmax+1 /))
          
    Lam = state(2*siz+1) + ii*state(2*varsize)
    
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Sets 'state' to the standard initial state value
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initial(state)
    implicit none

    real(kind=kind(1.0D0)) :: state(2*varsize), norm
    complex(kind=kind(1.0D0)) :: Phi(-nmax:nmax, -mmax:mmax), Lam

    select case(StateKeep)
        case(1)
            Phi = (1-ii)/(2*Sqrt(2*N_A)); Phi(0,0) = 1.0
            Phi(1,0)=0.0; Phi(0,1)=0.0; Phi(-1,0)=0.0; Phi(0,-1)=0.0
            Lam = 0.0

        case(0)
            open(13, file='teststate.txt', status='old', action='read')
            read(13,*) state
            
            call decode_state(state, Phi, Lam)
    end select

    norm = SUM(ABS(Phi)**2)
    Phi=Phi/sqrt(norm)
    call encode_state(Phi,Lam,state)

end subroutine

end module
