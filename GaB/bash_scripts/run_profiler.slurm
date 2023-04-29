#!/bin/bash


### Specify "wallclock time", hhh:mm:ss. Required field
#SBATCH --time 00:03:00


### Launch cuda module
### module load cuda11/11.0

### collect profile data
nvprof --log-file profile.txt ./bin/main ~/ece569/proj/source/GaB/IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu 10 128 100 0
nvprof -f --metrics all -o nvvp_profile.nvvp ./bin/main ~/ece569/proj/source/GaB/IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu 10 128 100 0
nvprof -f -o profile.timeline ./bin/main ~/ece569/proj/source/GaB/IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu 10 128 100 0
nvprof -f -o profile.metrics --analysis-metrics ./bin/main ~/ece569/proj/source/GaB/IRISC_dv4_R050_L54_N1296_Dform IRISC_dv4_R050_L54_N1296_Dform_Res_cu 10 128 100 0

### Launch profiler
nvvp $