#!/bin/bash

CONTROL=params

StateKeep=1
TestFlag=1
RunType='sweep_updown'

nmax=1

E_0=40
Kappa=8.1 
N_A=100000 
Omega_r=0.05

OmegaTil_init=0.5
OmegaTil_end=1
OmegaTil_no=40

Pump_init=0.1
Pump_end=1.0
Pump_no=600

    cat << EOF > $CONTROL
$StateKeep, $TestFlag, $RunType
$nmax, $E_0, $Kappa, $N_A, $Omega_r
$Pump_init, $Pump_end, $Pump_no
$OmegaTil_init, $OmegaTil_end, $OmegaTil_no
EOF

    ./main
    
rm params

