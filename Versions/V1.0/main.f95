program main
!Given initial state defined by mat and lam, this program time evolves
!the state of the atoms and the light using the cavity differential
!equation defined in state_control.f95.
    use utility
    use params
    use state_control
    use nag_library, only: d02bjf, d02bjw !bjw to remove condition for stopping.
    implicit none
   
    !parameters used
    real(kind=kind(1.0D0)) :: xi, xend, tol, W(40*varsize),&
                              state(2*varsize)
    integer                :: NN, ifail = 0

    !variables relevant to problem:
    complex(kind=kind(1.0D0)) :: mat(-nmax:nmax,-mmax:mmax), lam

    !interface with the routine that handles intermediate output from
    !NAG routine.
    interface
        subroutine output(xsol,y)
            real(kind=kind(1.0D0)), intent(inout) :: xsol
            real(kind=kind(1.0D0)), intent(in) :: y(*)
        end subroutine
    end interface

    !initial conditions:
    mat = 1.0E-2; mat(0,0) = 1.0; lam = 0.0

    !encode into initial state:
    call encode_state(mat, lam, state)

    !set up NAG routine:
    xi = 0.0; xend = 5000; NN = 2*varsize; tol = 0.0001;

    !Call routine, it gives intermediate output according to 
    !the specification of the "output" routine.
    call d02bjf(xi, xend, NN, state, diffeqn, tol, 'D',&
                output, d02bjw, W, ifail)

end program

!This routine is called internally by D02BJF and handles the
!intermediate output options. The solution state, y, must notbe changed
!and xsol must be set to the next time we wish to call output.

!This version of output prints the state to a file at 0.1 intervals of
!the independent variable.
subroutine output(xsol, y)
    use state_control
    implicit none
    
    real(kind=kind(1.0D0)), intent(inout) :: xsol
    real(kind=kind(1.0D0)), intent(in) :: y(*)
    
    !Use these for any data manipulation we may want.
    complex(kind=kind(1.D0)) :: mat(-nmax:nmax, -mmax:mmax), lam
    logical :: extant

    call decode_state(y, mat, lam)

    !Opens data file and writes to it; creates a new file and
    !writes if it doesn't exist already.
    inquire(file='data1.txt', exist=extant)

    if(extant) then
        open(1, file='data1.txt', status='old', action='write')
    else
        open(1,file='data1.txt',status='new', action='write')
    end if

    inquire(file='data2.txt', exist=extant)

    if(extant) then
        open(2, file='data2.txt', status='old', action='write')
    else
        open(2,file='data2.txt',status='new', action='write')
    end if
    
    !This part writes the data in columns for reading by GNUPlot.
    write(1,*) xsol, abs(lam)**2
    write(2,*) xsol, real(mat(0,0)), aimag(mat(0,0)), abs(mat(0,0))
    !Now set xsol to the next point we want output.
    xsol = xsol + 0.01
        
end subroutine
