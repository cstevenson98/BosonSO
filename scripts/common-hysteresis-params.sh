# Parameters common to the hysteresis pictures

setup_other_params() {

CONTROL=params
    
StateKeep=1
TestFlag=1
RunType='sweep_updown'

OmegaTil_init=$OmegaTil
OmegaTil_end=$OmegaTil
OmegaTil_no=1

Pump_init=0.1
Pump_end=1.0
Pump_no=900

tmin=0.0
tmax=5.0E6
tstep=2

convergence_threshold=0.002

}

. common-functions.sh

