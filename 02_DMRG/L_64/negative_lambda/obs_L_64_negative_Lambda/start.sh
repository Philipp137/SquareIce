for L in 64
do
    mkdir L$L
    cd  L$L
    for mu in 1.91 1.93 1.925 1.92 1.29 1.3 1.28  0.291 0.313 0.319
        do
        mkdir mu$mu
        cd mu$mu
        #cp /lustre/fs23/group/nic/stornati/bash_script .
        cp ../../bash_script .
        touch scaling
        echo " julia ../../run_dmrg_and_compute_obs.jl  $L  $mu  " >> bash_script 
        sbatch bash_script
        cd ..
        sleep 0.01
    done
    cd ..
done
