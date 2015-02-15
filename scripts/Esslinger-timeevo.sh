#!/bin/bash

# This script explores the irregular behaviour seen for some of the
# parameters based on Esslinger


E_0=40
Kappa=8.1 
N_A=100000 
Omega_r=0.05

nmax=10
. common-time-evo-params.sh

main()
{
    NEWDIR="Esslinger_Chaos_Omtil=${OmegaTil}_P=${P}_n=${nmax}"

    setup_new_directory

    setup_other_params
    create_and_run_script

    

    cd ..
}


# Flickering values
OmegaTil=1.9
P=0.469

main

# Very irregular
OmegaTil=0.5
P=0.261

