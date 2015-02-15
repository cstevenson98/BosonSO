module params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This module contains all of the physical and 
!mathematical contants related to the problem.
!These are then able to be used globally.
!
!   nmax, mmax - x, y extent of momentum-space
!                matrix.        
!
!   varsize    - the number of complex variables
!                in problem.
!
!   other variable pending understanding...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    implicit none

    integer, public         :: nmax, mmax, varsize
    real(kind=kind(1.0D0)), public   :: UN, Omega_r, Kappa
    real(kind=kind(1.0D0)), public   :: OmegaTil, Pump
     
contains

subroutine init_params()
    
    write(*,*) 'Size of matrix, nmax:'
    read(*,*) nmax
    mmax = nmax
    varsize = (2*nmax+1)*(2*mmax+1)+1

    write(*,*) 'UN, Omega_r, Kappa:'
    read(*,*) UN, Omega_r, Kappa

    write(*,*) 'OmegaTil, Pump:'
    read(*,*) OmegaTil, Pump
    
end subroutine

subroutine init_params_a()
    implicit none

    real(kind=kind(1.0D0)) :: gsqrtN, Omega
  
    write(*,*) 'Size of matrix, nmax:'
    read(*,*) nmax
    mmax = nmax
    varsize = (2*nmax+1)*(2*mmax+1)+1

    write(*,*) 'UN, Omega_r, Kappa:'
    read(*,*) UN, Omega_r, Kappa

    write(*,*) 'Omega, gsqrt(N):'
    read(*,*) Omega, gsqrtN

    Pump = 0.5*gsqrtN/(Omega_r*abs(UN))
    OmegaTil = Omega/UN + 0.5
end subroutine
end module
