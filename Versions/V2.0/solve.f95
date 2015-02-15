!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Summer 2014 - C.Stevenson, St. Andrews 110012872
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This module contains routines that evolve a state
!according to the cavity eom and control the output
!of the resulting data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module solve
    use params
    use state_control
    implicit none

    real(kind=kind(1.0D0)), save  :: prevstate(1:10)
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine evolves the state from xinit to 
!xend according to the dynamics of eom, the 
!equations of motion of the system.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine evolve_state(state, xinit, xend, eom)
    use params
    use state_control
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
    NN = 2*varsize; tol = 0.00001;

    !Now evolve the equation to the end point.
    !Output and stopping both suppressed.
    call d02bjf(xinit, xend, NN, state, eom, tol, 'D', &
                d02bjx, d02bjw, W, ifail)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This routine evolves the state acording to eom
!until a steady state is reached, indicated by 
!the modulus of the derivative being below a 
!certain tolerance.
!
!If the state does not become steady by xend then
!this routine returns a stop error.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine next_steady(state, xinit, xend, eom, step, output)

    implicit none
    real(kind=kind(1.0D0))  :: state(2*varsize), dstate(2*varsize)
    real(kind=kind(1.0D0))  :: xinit, xend, step, store
    integer                 :: i
    character(len=*)        :: output 

    interface
        subroutine eom(t, s, ds_dt)
            real(kind=kind(1.0D0)), intent(in) :: t, s(*)
            real(kind=kind(1.0D0)), intent(out) :: ds_dt(*)
        end subroutine
    end interface

    if(output .ne. '') then
     call write_state(xinit, state, output)
    end if

    !This first loop doesn't check for convegence as
    !the initial state may change slowly.
    do i=1, 10000
        store = xinit + step
        call evolve_state(state, xinit, store, eom)
        call storage(state, 2*varsize-1, prevstate)
        xinit = store

        if(output .ne. '') then
            call write_state(xinit, state, output)
        end if
    end do

    !Now begin the check for convergence and stop the
    !the calculation once state is steady, or stop if
    !maximum time is reached.    
    do while(.not. converged())
        store = xinit + step
        call evolve_state(state, xinit, store, eom)
        call storage(state, 2*varsize-1, prevstate)
        xinit = store

        if(output .ne. '') then
            call write_state(xinit, state, output)
        end if

        if(xinit .ge. xend) then
            stop 'Error: State has not converged within alloted time.'
        end if
    end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This will store a chosen element of the state in
!an array which is used to determine if the state
!is steady, ie if all 10 previous elements'
!derivatives are below a certain tolerance, then 
!"converged" flags true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine storage(state, element, list)
    implicit none
    
    real(kind=kind(1.0D0)) :: state(2*varsize), dstate(2*varsize), &
                              list(10)
    integer                :: i, element

    call diffeqn(state(1), state, dstate)

    do i=0,8
        list(10-i) = list(9-i)
    end do
    
    list(1) = dstate(element)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This function checks that the previous 10 states
!are below some tolerance.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function converged()

    implicit none
    logical                 :: converged
    integer                 :: i
    
    converged = .true.

    do i=1, 10
        if(abs(prevstate(i)) .gt. 0.0001) then
            converged = .false.
            exit
        end if
    end do

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This routine dictates the where the state is 
!outputted; the file, output. It writes the time
!and the relevant variables in the problem such 
!as the light field and the momentum components
!and potentially combinations thereof. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_state(t, state, output)
    implicit none

    
    real(kind=kind(1.0D0))       :: t, state(2*varsize)
    character(len=*)             :: output
    
    complex(kind=kind(1.0D0))    :: mat(-nmax:nmax, -mmax:mmax), lam
    logical :: extant 

    character(len = 41)          :: datapath = "/home/conchops/fortranfiles/project/Data/"
    
    datapath = trim(datapath)
    output = trim(output)

    call decode_state(state, mat, lam)

    output = trim(output)

    !Opens data file and writes to it; creates a new file and
    !writes if it doesn't exist already.
    inquire(file=datapath//output, exist=extant)

    if(extant) then
        open(1, file=datapath//output, status='old', action='write')
    else
        open(1,file=datapath//output,status='new', action='write')
    end if
    
    !This part writes the data in columns for reading by GNUPlot.
    write(1,*) t/1000.0, (N_A*Omega_r/E_0)*abs(lam)**2, 0.5*(abs(2*mat(1,1))**2-abs(mat(0,0))**2), abs(mat(1,0))**2

end subroutine

end module
