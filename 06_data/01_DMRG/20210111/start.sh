#!/bin/bash
for L in 6 11 21 41 81 161
do
    #mkdir L$L
    #cd  L$L
    for Lambda in -10.0 -5.0 -1.0 -0.5 -0.2 -0.1 0 0.1 0.2 0.5
        do
        #mkdir Lambda$Lambda
        #cd Lambda$Lambda
        #cp /lustre/fs23/group/nic/stornati/bash_script .
        #cp ../../bash_script .
	    #touch scaling
        echo " julia main.jl  $L  160 $Lambda  >> dmrg.log " > bash_script
        sh bash_script
        #cd ..
        sleep 1
    done
    #cd ..
done
