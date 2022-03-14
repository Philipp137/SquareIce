#!/bin/bash
for L in 25
do
    #mkdir L$L
    cd  L$L
    for mu in  0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.41 0.42 0.43
    do
        cd mu$mu
        cp /lustre/fs23/group/nic/stornati/bash_script .
        touch scaling
        echo " /lustre/fs23/group/nic/stornati/julia/julia /lustre/fs23/group/nic/stornati/SquareIce_new/02_DMRG/measure.jl  $L $L 2 -1 $mu 'SQ' >> scaling " >> bash_script
        sbatch bash_script
        cd ..
        sleep 0.05
    done
    cd ..
done
