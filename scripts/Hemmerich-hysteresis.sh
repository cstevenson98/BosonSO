#!/bin/bash

#This script explores hystersis for Hemmerich's parameters
. common-hysteresis-params.sh

E_0=0.025
Kappa=0.005 
N_A=100000 
Omega_r=0.003

convergence_threshold=0.002

main()
{
    NEWDIR="Hemmerich_OmTil_$1-n_$2"

    setup_new_directory

    nmax=$2
    OmegaTil=$1

    echo $OmegaTil $nmax

    setup_other_params
    create_and_run_script

    cleanup_data

    cd ..
}

for omtil in 2.1 2.5 3.0 1.9 1.5
do
    for n in 1 2 3 5 10
    do
	main $omtil $n 
    done
done

