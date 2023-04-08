#!/bin/bash


### Specify "wallclock time", hhh:mm:ss. Required field
#SBATCH --time 00:03:00

### set directory for job execution, replace the following path with your own path
exec_path=~matthewmartin/ece569/proj/bin

### run profile
/usr/bin/time mpirun -np 1 nvprof --log-file profile.txt ./IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu 10 32 0 0
/usr/bin/time mpirun -np 1 nvprof --metrics all -o nvvp_profile.nvvp ./IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu 10 32 0 0
/usr/bin/time mpirun -n 1 nvprof -o profile.timeline ./IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu 10 32 0 0
/usr/bin/time mpirun -n 1 nvprof –o profile.metrics --analysis-metrics ./IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu 10 32 0 0

