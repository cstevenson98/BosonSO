#!/bin/bash

#This script will run various runs to produce large amounts of data.

main()
{
    NEWDIR="OmTil_$1-n_$2"
    mkdir $NEWDIR
    cp main $NEWDIR/ 
    cd $NEWDIR/
    
    CONTROL=params
    
    StateKeep=1
    TestFlag=1
    RunType='sweep_updown'

    nmax=$2
    E_0=40
    Kappa=8.1 
    N_A=100000 
    Omega_r=0.05

    OmegaTil_init=$1
    OmegaTil_end=1
    OmegaTil_no=1

    Pump_init=0.1
    Pump_end=1.0
    Pump_no=900

    cat << EOF > $CONTROL
$StateKeep, $TestFlag, $RunType
$nmax, $E_0, $Kappa, $N_A, $Omega_r
$Pump_init, $Pump_end, $Pump_no
$OmegaTil_init, $OmegaTil_end, $OmegaTil_no
EOF

    echo $OmegaTil_init $nmax

    ./main

    tar -cvzf "$NEWDIR-Evol.tar.gz" UP* DN*
    rm UP* DN*

    cd ..
}

for COUNTER in 6 7 8 9 10
do
for n in 1 2 3 5 10 
do
main $(echo 0.1*$COUNTER+2.000 | bc) $n
done
done

