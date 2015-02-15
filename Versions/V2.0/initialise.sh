#!/bin/bash

CONTROL=params

run() {

    cat << EOF > $CONTROL
$Multiple, $StateKeep, $PumpEvolve, $Up
$nmax, $E_0, $Omega_0, $Kappa,$N_A
$Omega_init, $Omega_end, $Omega_step
$gsqrtN_init, $gsqrtN_end, $gsqrtN_step
EOF

continuefile="/home/conchops/fortranfiles/project/continue.dat"

for x in `seq 1 1 3`;
#while [ -e $continuefile ];
do
    ./main < $CONTROL
    #Something here to decide loading bar.
done

}

Multiple=1
StateKeep=1
PumpEvolve=0
Up=1

nmax=1

E_0=40
Omega_0=0.05
Kappa=8.1
N_A=100000 #1E5??

Omega_init=10
Omega_end=10
Omega_step=100

gsqrtN_init=0
gsqrtN_end=2
gsqrtN_step=3

run
