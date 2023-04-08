#!/bin/bash


### Specify "wallclock time", hhh:mm:ss. Required field
#SBATCH --time 00:03:00

### set directory for job execution, replace the following path with your own path
bin_path=~matthewmartin/ece569/proj/bin

cd $bin_path
###



### run sweep
for i in {0..10..2}
do
  time ./main ~/ece569/proj/source/GaB/IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu $i
done



