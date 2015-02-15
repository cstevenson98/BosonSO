# Parameters common to the time evo cases pictures

setup_other_params() {

CONTROL=params
    
StateKeep=1
TestFlag=1
RunType='single'

OmegaTil_init=$OmegaTil
OmegaTil_end=$OmegaTil
OmegaTil_no=1

Pump_init=$P
Pump_end=$P
Pump_no=900

tmin=0.0
tmax=5.0E6
tstep=2

convergence_threshold=0.002
}
. common-functions.sh

