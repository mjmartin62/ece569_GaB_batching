#!/bin/bash

###  Runs nested loops for multiple CUDA configurations
# i = number of CUDA STREAMS
# j = BLOCK SIZE
# k = number of words per CUDA STREAM

### run sweep
for i in {5,20,35}
do
    for j in {32,128,512,1024}
    do
        for k in {5,20,35}
        do
            ./bin/main ~/ece569/proj/source/GaB/IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu $i $j $k 0
        done
    done
done



