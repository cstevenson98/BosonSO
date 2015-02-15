!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Summer 2014 - C.Stevenson, St. Andrews 110012872
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!The program was written to simulate the boson self-organanisation 
!system studied in PRA Bhaseen, et.al. This will perform sweep in 
!various directions: "up", "down", "updown" and "downup".
!
!It outputs into several files: the sweep file, which contains the
!steady state solutions a varying pump strength, with each run's 
!initial conditions being the steady state from the previous run, 
!and the time evolution of each indulvidual run also.
!
!This program is easily extensible to include additional possible run
!configurations as well as extending the problem to include mutliple
!momentum states.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program boson_self_organisation
    use utility
    use params
    use evolution
    implicit none

    !The names of files where each individual run, dataoutfile;
    !the swept pump results, sweepresults.
    character(len=100)        :: DataOutFile, SweepResults

    !Momentum matrix and light field varioables Phi, and Lam.
    complex(kind=kind(1.0D0)), allocatable :: Phi(:,:), Lam

    !CALL ALL SETUP ROUTINES -- SETTING INITIAL CONDITION!
  !========================================================!

    call setup()
    call initial(CurrentState)


    !MAIN LOOP OF PROGRAM -- DOING EACH SET OF DATA WANTED!
  !=========================================================!

        !FIRST -- CHECK WHAT TYPE OF RUN TO DO:
    !===========================================!

    select case(RunType)
        case('sweep_up')
            call execute_sweep('up')

        case('sweep_down')
            call execute_sweep('down')

        case('sweep_updown')
            call execute_sweep('updown')

        case('sweep_downup')
            call execute_sweep('downup')

        case('varying')
            call varying_pump()

        case('single')
            call execute_single()
                      
        case default
            call execute_single()

    end select

contains 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!This routine hold the various permutations desired for the possible 
!sweeps one may wish to perform.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine execute_sweep(SweepType)
    use params
    use evolution
    implicit none
    character(len=*), intent(in) :: SweepType

    select case(SweepType)
        case('updown') 
            
            call pump_sweep_block('up')
            call pump_sweep_block('down')
            call change_omegatil('up')
            call initial(CurrentState)

            do while(OmegaTil .lt. A_OmegaTil(2))
                call pump_sweep_block('up')
                call pump_sweep_block('down')
                call change_omegatil('up')
                call initial(CurrentState)
            end do

        case('up')

            call pump_sweep_block('up')
            call change_omegatil('up')
            call initial(CurrentState)

            do while(OmegaTil .lt. A_OmegaTil(2))
                call pump_sweep_block('up')
                call change_omegatil('up')
                call initial(CurrentState)
            end do

        case('down')

            Pump = A_Pump(2); OmegaTil = A_OmegaTil(2)
            call pump_sweep_block('down')
            call change_omegatil('up')
            call initial(CurrentState)

            do while(OmegaTil .lt. A_OmegaTil(2))
                call pump_sweep_block('down')
                call change_omegatil('up')
                call initial(CurrentState)
            end do

        case('downup')

            Pump = A_Pump(2); OmegaTil = A_OmegaTil(2)
            call pump_sweep_block('down')
            call pump_sweep_block('up')
            call change_omegatil('up')
            call initial(CurrentState)

            do while(OmegaTil .lt. A_OmegaTil(2))
                call pump_sweep_block('down')
                call pump_sweep_block('up')
                call change_omegatil('up')
                call initial(CurrentState)
            end do

    end select
    
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!This contains the basic operations to carry out one sweep of the pump
!either up or down.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pump_sweep_block(BlockType)
    use params
    use evolution
    implicit none

    character(len=*)  :: BlockType
    complex(kind=kind(1.0D0)) :: Phi(-nmax:nmax, -mmax:mmax), Lam

    do while(cont_iteration(Pump, A_Pump, Blocktype))
        write(*,*) Pump, OmegaTil
        call fileset(BlockType)
        call next_steady(CurrentState,tinit,tend,diffeqn,tstep, 1000)
        call decode_state(CurrentState, Phi, Lam)
        call write_output(0.0D0, 'sweep')
        close(51); close(52)
        call change_pump(BlockType)
    end do

end subroutine
    !Function to determine the continuation condition required above.
function cont_iteration(Pump, A_Pump, BlockType)

    logical :: cont_iteration
    real(kind=kind(1.0D0)) :: Pump, A_Pump(3)
    character(len=*)       :: BlockType

    select case(blocktype)
    case('up')
        cont_iteration=(Pump .lt. A_Pump(2))
    case('down')
        cont_iteration=(Pump .gt. A_Pump(1))
    case default
        cont_iteration=.false.
    end select 

end function

subroutine varying_pump() !################IMPLEMENT
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!This routine simply performs an evolution for the first value of P
!in the array provided to the program.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine execute_single()
    use params
    use evolution    
    implicit none
                    
    write(*,*) Pump, OmegaTil
    call fileset('single')
    call next_steady(CurrentState,tinit,tend,diffeqn,tstep, 1000)
    close(51); close(52)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Calling this routine will set the names of the output files, both
!for the state evolution and the sweep evolution data.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fileset(direction) !REQUIRES CHANGES ##################
    implicit none

    character(len=*)       :: direction
    logical                :: extant

    ! Set up file names
    select case(direction)
        case('up')
            write(TimeOutput,"(a8, f6.3, a9, f6.3, a4)") &
                 & "UP-Pump_", Pump, "OmegaTil_", OmegaTil, ".dat"
        case('down')
            write(TimeOutput,"(a8, f6.3, a9, f6.3, a4)") &
                 & "DN-Pump_", Pump, "OmegaTil_", OmegaTil, ".dat"
        case('single')
            write(TimeOutput,"(a8, f6.3, a9, f6.3, a4)") &
                 & "SI-Pump_", Pump, "OmegaTil_", OmegaTil, ".dat"
    end select

    write(SweepOutput,"(a9, f6.3, a3, I2, a4)") &
         & "OmegaTil_", OmegaTil,"-n_", nmax, ".dat"
    
    TimeOutput = trim(TimeOutput); SweepOutput = trim(SweepOutput)
    
    ! Open files, print headers
    inquire(file=TimeOutput, exist=extant)
    if(extant) then
        open(51, file=TimeOutput, status='old', access='append')
    else
        open(51, file=TimeOutput, status='new', action='write')
    end if
    write(51,'("#",6(A14,X))') "Time(ms)", "R[Lam]", "I[Lam]", &
         &   "Sz", "IPR", "|Phi(1,0)|^2"


    inquire(file=SweepOutput, exist=extant)
    if(extant) then
        open(52, file=SweepOutput, status='old', access='append')
    else
        open(52, file=SweepOutput, status='new', action='write')
    end if
    write(52,'("#",8(A14,X))') "Pump", "|Lam|^2", "R[Lam]", "I[Lam]", &
         &   "R[Sp]", "I[Sp]", "Sz", "IPR"
    
end subroutine
end program
