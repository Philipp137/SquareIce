#!/bin/bash
for L in 8 16 24 32 40 48 56 64 72 80 96 120 
do
 mkdir L$L
 cd L$L
 L2=$(expr $L / 2)
 cp ../xylargeQ3 .
 cp ../submitQ3 . 
 cp ../XY.job .
 sed -i "s/LS/${L}/g" submitQ3
 sed -i "s/QQC/0/g" submitQ3
 sed -i "s/LL2/${L2}/g" submitQ3
 mv XY.job XYL$L.job
 qsub XYL$L.job
 sleep 5
 cd .. 
done
