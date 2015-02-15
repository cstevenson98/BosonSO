program main
!Given initial state defined by Phi and Lam, this program time evolves
!the state of the atoms and the light using the cavity differential
!equation defined in state_control.f95.
    use utility
    use params
    use state_control
    use solve

    !Light field.
    complex(kind=kind(1.0D0))              :: Lam
    !Momentum matrix.
    complex(kind=kind(1.0D0)), allocatable :: Phi(:,:)
    !Flattened, encoded state.
    real(kind=kind(1.0D0)), allocatable    :: state(:)
    !Initial and final time units.
    real(kind=kind(1.0D0))                 :: xinit, xend

    call init_params_a()

    allocate(Phi(-nmax:nmax, -mmax:mmax)); allocate(state(2*varsize))

    !initial conditions:
    Phi = (1-ii)/Sqrt(10E5); Phi(0,0) = 1.0
    Phi(1,0)=0.0; Phi(0,1)=0.0; Phi(-1,0)=0.0; Phi(0,-1)=0.0
    Lam = 0.0
    xinit = 0.0; xend = 20000

    !encode into our initial state
    call encode_state(Phi, Lam, state)

    !Evolve the state:
    call evolve_state(state, xinit, xend, diffeqn)

end program
