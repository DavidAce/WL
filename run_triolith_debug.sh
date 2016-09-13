#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-00:05:00
#SBATCH -N 1
#SBATCH --exclusive

valgrind --tool=memcheck --leak-check=full -v mpprun ./build/Debug/WL
