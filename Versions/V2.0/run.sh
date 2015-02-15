#!/bin/bash

CONTROL=params

run() {

    cat << EOF > $CONTROL
$nmax
$UN, $Omega_0, $Kappa
$Omega,$gsqrtN
EOF


    ./main < $CONTROL

}

nmax=1
UN=40
Omega_0=0.05
Kappa=8.1

Omega=-30
gsqrtN=0.9

run
