#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-00:05:00
#SBATCH -N 1
#SBATCH --exclusive
export OMP_NUM_THREADS=4
valgrind --tool=memcheck --leak-check=full -v mpprun -l ./build/Debug/WL
