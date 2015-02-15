module solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This module contains routines that evolve a state
!according to the cavity eom and control the output
!of the resulting data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use params
    use state_control
    implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine evolves the state from xinit to 
!xend according to the dynamics of eom, the 
!equations of motion of the system.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine evolve_state(state, xinit, xend, eom)
    use params
    use state_control
    use nag_library, only: d02bjf, d02bjw !bjw to remove condition for stopping.
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

    !Call routine, it gives intermediate output according to 
    !the specification of the "output" routine.
    call d02bjf(xinit, xend, NN, state, eom, tol, 'D',&
                output, d02bjw, W, ifail)

end subroutine

!This routine is called internally by D02BJF and handles the
!intermediate output options. The solution state, y, must notbe changed
!and xsol must be set to the next time we wish to call output.
!
!This version of output prints the state to a file at intervals of
!the independent variable.
subroutine output(xsol, y)

    logical :: extant
    
    real(kind=kind(1.0D0)), intent(inout) :: xsol
    real(kind=kind(1.0D0)), intent(in) :: y(*)
    real(kind=kind(1.0D0))             :: dy(2*varsize)    
    
    call diffeqn(xsol, y, dy)

    !Opens data file and writes to it; creates a new file and
    !writes if it doesn't exist already.
    inquire(file='data.txt', exist=extant)

    if(extant) then
        open(1, file='data.txt', status='old', action='write')
    else
        open(1,file='data.txt',status='new', action='write')
    end if
    
    !This part writes the data in columns for reading by GNUPlot.
    write(1,*) (xsol/Omega_r)/10E3, y(2*varsize-1) , dy(2*varsize-1)

    !Now set xsol to the next point we want output.
    xsol = xsol + 0.1
        
end subroutine

end module
