!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Summer 2014 - C.Stevenson, St. Andrews 110012872
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This module contains routines that evolve a state
!according to the cavity eom and control the output
!of the resulting data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module evolution
    use params
    use eom_so
    implicit none
    
    integer, parameter      :: BackCheck = 40
    real(kind=kind(1.0D0))  :: StateStore(BackCheck, 1:2)
    
    ! WHether to dump info relating to convergence or not
    logical :: dbg_convergence

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!This subroutine evolves the state from xinit to xend according to 
!the dynamics of eom; the equations of motion of the system.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine evolve_state(state, xinit, xend, eom)
    use nag_library, only: d02bjf, d02bjw, d02bjx
    implicit none
    
    real(kind=kind(1.0D0)), intent(inout)    :: state(2*varsize)
    real(kind=kind(1.0D0)), intent(inout)    :: xinit, xend

    real(kind=kind(1.0D0)) :: tol, W(40*varsize)
    integer                :: i, NN, ifail = 0

    interface
        subroutine eom(t, s, ds_dt)
            real(kind=kind(1.0D0)), intent(in) :: t, s(*)
            real(kind=kind(1.0D0)), intent(out) :: ds_dt(*)
        end subroutine
    end interface

    !set up NAG routine:
    NN = 2*varsize; tol = 0.000001;

    !Now evolve the equation to the end point.
    !Output and stopping both suppressed.
    call d02bjf(xinit, xend, NN, state, eom, tol, 'D', &
              & d02bjx, d02bjw, W, ifail)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!This routine evolves the state acording to eom until a steady state 
!is reached, indicated by the modulus of the derivative being below a 
!certain tolerance.
!
!If the state does not become steady by xend then this routine 
!returns a stop error.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine next_steady(state, xinit, xend, eom, step, buffer)
    implicit none
    real(kind=kind(1.0D0))  :: state(2*varsize), dstate(2*varsize)
    real(kind=kind(1.0D0))  :: xinit, xend, tempinit, tempend, step, store
    integer                 :: i, buffer

    interface
        subroutine eom(t, s, ds_dt)
            real(kind=kind(1.0D0)), intent(in) :: t, s(*)
            real(kind=kind(1.0D0)), intent(out) :: ds_dt(*)
        end subroutine
    end interface
    
    tempinit=xinit; tempend=xend

    call write_output(tempinit, 'time')

    !This first loop doesn't check for convegence as
    !the initial state may change slowly.
    do i=1, buffer
        store = tempinit + step
        call evolve_state(state, tempinit, store, eom)
        call storage()
        tempinit = store
        call write_output(tempinit, 'time')
    end do

    !Now begin the check for convergence and stop the
    !the calculation once state is steady, or stop if
    !maximum time is reached.    
    time_evo_loop: do while(.not. converged())
        store = tempinit + step
        call evolve_state(state, tempinit, store, eom)
        call storage()
        tempinit = store
        call write_output(tempinit, 'time')
        if(tempinit .ge. tempend) then
           ! Write warning, but carry on
           write(51,'("# Warning did not reach a steady state")')
           write(52,'("# Warning, non-converged point")')
           exit time_evo_loop
           !stop 'Error: State has not converged within alloted time.'
        end if
     end do time_evo_loop
    
     if (dbg_convergence) then
        do i = 1, 10
           write(*,*) StateStore(i, 1), StateStore(i, 2)
        end do
     end IF
     
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This will store a chosen element of the state in
!an array which is used to determine if the state
!is steady, ie if all 10 previous elements'
!derivatives are below a certain tolerance, then 
!"converged" flags true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine storage()
    implicit none

    complex(kind=kind(1.0D0)) :: Phi(-nmax:nmax, -mmax:mmax), Lam
    integer                   :: i
    
    call decode_state(CurrentState, Phi, Lam)
    
    do i = 1, (BackCheck - 1)
        StateStore(i,:) = StateStore(i+1,:)
    end do
    
    StateStore(BackCheck,1) = Real(Lam)
    StateStore(BackCheck,2) = Aimag(Lam)
        
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This function checks that the previous 10 states
!are below some tolerance.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function converged()
    implicit none
    logical                 :: converged
    integer                 :: i
    real(kind=kind(1.0D0))  :: maximum_R, maximum_I, minimum_R, minimum_I
    
    converged = .true.
    maximum_R = maxval(StateStore(:,1)); maximum_I = maxval(StateStore(:,2))
    minimum_R = minval(StateStore(:,1)); minimum_I = minval(StateStore(:,2))
    
    if( ( (maximum_R - minimum_R) .gt. convergence_threshold) .or. &
       &( (maximum_I - minimum_I) .gt. convergence_threshold)) then
            converged = .false.
    end if
    
    !insert decreasing criterion here.
!    if(maxval(StateStore(:BackCheck/2,1)) .10)
       
end function

end module
