!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Summer 2014 - C.Stevenson, St. Andrews 110012872
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Given initial state defined by Phi and Lam, this
!program time evolves the state of the atoms and 
!the light using the cavity differential equation
!defined in state_control.f95.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
    use utility
    use params
    use state_control
    use solve
    implicit none

    character(len=30) :: outputfile
    real(kind=kind(1.0D0)) :: xinit, xend, step

    call setup()

    
    xinit = 0.0; xend = 500000.0; step = 2.0

!    write(*,*) gsqrtN, Omega

    call fileset(outputfile, gsqrtN, Omega)

    write(*,*) (N_A*Omega_r/E_0)*Abs(CurrentState(2*varsize-1)+ii*CurrentState(2*varsize))**2

    call next_steady(CurrentState, xinit, xend, diffeqn, step, outputfile)

    write(*,*) (N_A*Omega_r/E_0)*Abs(CurrentState(2*varsize-1)+ii*CurrentState(2*varsize))**2

    open(21, file='longstate.txt', status='old', access='append')

    write(21,*) gsqrtN, (N_A*Omega_r/E_0)*Abs(CurrentState(2*varsize-1)+ii*CurrentState(2*varsize))**2
    
    close(21)

    call make_continue()

!!!!!!!!!!!!!!!!!
!WORKING VERSION
!!!!!!!!!!!!!!!!!
!    call setup()

!    Pump = A_Pump(1); OmegaTil = A_OmegaTil(1); Omega_r = A_Omega_r(1)
!    allocate(Phi(-nmax:nmax, -mmax:mmax)); allocate(state(2*varsize))
!    call fileset(outfile, Pump, OmegaTil)

!    !initial conditions:
!    Phi = (1-ii)/(2*Sqrt(2*N_A)); Phi(0,0) = 1.0
!    Phi(1,0)=0.0; Phi(0,1)=0.0; Phi(-1,0)=0.0; Phi(0,-1)=0.0
!    Lam = 0.0
!    xinit = 0.0; xend = 10000; step = 2


!    !encode into our initial state
!    call encode_state(Phi, Lam, state)

!    !Evolve the state:
!    call next_steady(state, xinit, xend, diffeqn, step, outfile)


!Check which procedure to run, ie Multiple, or Not and Then 

end program

subroutine fileset(variable, par1, par2)
    implicit none

    character(len=*)      :: variable
    real(kind=kind(1.0D0)) :: par1, par2

    write(variable,"(a8, f5.3, a6, f6.3, a4)") "CS-gsqrtN_", par1, "Omega_", par2, ".txt"
!    write(*,*) variable
end subroutine
