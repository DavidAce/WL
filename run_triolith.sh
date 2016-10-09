#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-0:40:00
#SBATCH -N 6
#SBATCH -e error_file.e
#SBATCH -o output_file.o
#SBATCH --exclusive
export OMP_NUM_THREADS=1
mpprun ./build/Release/WL
