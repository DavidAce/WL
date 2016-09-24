#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-03:00:00
#SBATCH -N 4
#SBATCH --exclusive
export OMP_NUM_THREADS=1
mpprun ./build/Release/WL
