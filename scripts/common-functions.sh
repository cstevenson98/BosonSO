setup_new_directory() {

    mkdir -p $NEWDIR
    cd $NEWDIR/

    ln -sf ../../main ./


}

create_and_run_script() {
    cat << EOF > $CONTROL
$StateKeep, $TestFlag, $RunType
$nmax, $E_0, $Kappa, $N_A, $Omega_r
$Pump_init, $Pump_end, $Pump_no
$OmegaTil_init, $OmegaTil_end, $OmegaTil_no
$tmin, $tmax, $tstep
$convergence_threshold
EOF


    ./main
}

cleanup_data() {
    tar -czf "$NEWDIR-Evol.tar.gz" UP* DN*
    rm UP* DN*
}
