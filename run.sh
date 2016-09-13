#!/bin/bash
option=$1
echo "Starting job"

if [[ "$HOSTNAME" == *"triolith"* ]]
then
    echo "We're on triolith!";
    sbatch run_triolith.sh ${option}
else
    echo "We're on my pc!"
    ./run_my_pc.sh ${option}
fi