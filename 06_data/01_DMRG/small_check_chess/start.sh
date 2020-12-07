#!/bin/bash
for L in 21 
do
  #  mkdir L$L
    cd  L$L
    for Lambda in +1.0 +0.75 +0.5  +0.25 -0.1  -0.3 -0.5 -0.7 -0.9 -1.0 -2.0 -3.0 -4.0 -5.0  
        do
        mkdir Lambda$Lambda
        cd Lambda$Lambda 
        #cp /lustre/fs23/group/nic/stornati/bash_script .
        cp ../../bash_script .
	touch scaling
        echo " /lustre/fs23/group/nic/stornati/julia/julia /lustre/fs23/group/nic/stornati/big_fss/fss.jl  $L  40 $Lambda  >> scaling " >> bash_script
        sbatch bash_script
        cd ..
        sleep 1
    done
    cd ..
done



