#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-00:05:00
#SBATCH -N 1
#SBATCH --exclusive
export OMP_NUM_THREADS=4
mpprun --pass="--bind-to-core --bysocket" ./build/Release/WL
