#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-02:00:00
#SBATCH -N 4
#SBATCH --exclusive
export OMP_NUM_THREADS=1
mpprun ./build/Release/WL
