#!/bin/bash
for L in 8 
do
    #mkdir L$L
    #cd  L$L
    for Lambda in -2 0 1
        do
          for i in {0..30}
            do
             
              mu=$(bc <<< "$i*0.05")
          
        #mkdir Lambda$Lambda
        #cd Lambda$Lambda
        #cp /lustre/fs23/group/nic/stornati/bash_script .
        #cp ../../bash_script .
	    #touch scaling
            logfile="dmrg_Lx_$L.log"
            command="julia ../../../02_DMRG/main.jl  $L  $((L*2)) 2 $Lambda $mu"
            echo -e "\n$command" 
            $command >> $logfile 

        
        #sh bash_script
        #cd ..
        sleep 1
      done
    done
    #cd ..
done
