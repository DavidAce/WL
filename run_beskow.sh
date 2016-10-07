#!/bin/bash

#SBATCH -J WL
#SBATCH -A 2016-10-17
#SBATCH -t 0-6:00:00
#SBATCH -N 6
#SBATCH --exclusive
export OMP_NUM_THREADS=1
aprun ./build/Release/WL
