#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-2:00:00
#SBATCH -N 6
#SBATCH -e error_file.e
#SBATCH -o cout_file.o
#SBATCH --exclusive
export OMP_NUM_THREADS=1
mpprun ./build/Release/WL
