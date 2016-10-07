#!/bin/bash

#SBATCH -J WL
#SBATCH -A 2016-10-17
#SBATCH -t 0-00:05:00
#SBATCH -N 1
#SBATCH -e error_file.e
#SBATCH -o output_file.o
#SBATCH --exclusive
export OMP_NUM_THREADS=4
valgrind --tool=memcheck --leak-check=full -v aprun ./build/Debug/WL
