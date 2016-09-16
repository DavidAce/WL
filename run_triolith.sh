#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-01:00:00
#SBATCH -N 2
#SBATCH --exclusive
export OMP_NUM_THREADS=1
mpprun ./build/Release/WL
