#!/bin/bash

echo "Starting job"

if [[ "$HOSTNAME" == *"triolith"* ]]
then
    echo "We're on triolith!";
    sbatch run_triolith.sh
else
    echo "We're on my pc!"
    ./run_my_pc.sh
fi