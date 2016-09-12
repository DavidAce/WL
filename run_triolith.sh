#!/bin/bash

#SBATCH -J WL
#SBATCH -t 0-00:15:00
#SBATCH -N 2
#SBATCH --exclusive

mpprun ./build/Release/WL
