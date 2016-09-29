#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-0:45:00
#SBATCH -N 1
#SBATCH --exclusive
export OMP_NUM_THREADS=1
mpprun ./build/Release/WL
