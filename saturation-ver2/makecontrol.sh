#!/usr/bin/env bash

mkdir ./BGKSTest

gcc ./Utilities/Auto_control.c -o Auto-Control -lm
gcc ./Utilities/Auto_main.c -o Auto_Run -lm
gcc ./Utilities/append-control.c -o Append -lm

./Auto-Control -dir ./BGKSTest -model 0 1 2 3 22 -qup 500 -lmass 0.0 -sudakov 1 -rfix 0


./Append "#define MASS_C2 1.69" ./BGKSTest/M*

./Append "#define STAR 1" ./BGKSTest/M*

./Append "#define N_SIMPS_R 100" ./BGKSTest/M*

./Append "#define PRINT_PROGRESS 1" ./BGKSTest/*Model3* 
