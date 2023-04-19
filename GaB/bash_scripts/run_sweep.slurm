#!/bin/bash


### Specify "wallclock time", hhh:mm:ss. Required field
#SBATCH --time 00:03:00

###



### run sweep
for i in {2..4..2}
do
  ./bin/main ~/ece569/proj/source/GaB/IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu $i 32 100 0

done



