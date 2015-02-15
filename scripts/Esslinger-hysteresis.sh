#!/bin/bash

# This script explores hysteresis for Esslinger's experimental parameters
. common-hystersis-params.sh

E_0=40
Kappa=8.1 
N_A=100000 
Omega_r=0.05

main()
{
    NEWDIR="Esslinger_OmTil_$1-n_$2"
    
    setup_new_directory

    nmax=$2
    OmegaTil=$1

    echo $OmegaTil $nmax

    setup_other_params
    create_and_run_script

    cleanup_data

    cd ..
}

for omtil in 3.5 4.0 4.5
do
    for n in 1 2 3 5 10 
    do
	main $omtil $n 
    done
done

