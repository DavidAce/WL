#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-00:05:00
#SBATCH -N 1
#SBATCH --exclusive

mpprun ./build/Release/WL
