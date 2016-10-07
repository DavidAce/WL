#!/bin/bash

#SBATCH -J WL
#SBATCH -A 2016-10-17
#SBATCH -t 0-0:01:00
#SBATCH -N 1
#SBATCH --exclusive
export OMP_NUM_THREADS=1
aprun -n 32 ./build/Release/WL > OutputFiles/out.o
