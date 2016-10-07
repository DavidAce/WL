#!/bin/bash
option=$1
echo "Starting job"
option=$1
if [[ "$HOSTNAME" == *"triolith"* ]]
then
    echo "We're on triolith!";
    if [[ "${option}" == *"valgrind"* ]]
    then
        sbatch run_triolith_debug.sh
    else
        sbatch run_triolith.sh
    fi
elif [[ "$HOSTNAME" == *"beskow"* ]]
then
    echo "We're on beskow!";
    if [[ "${option}" == *"valgrind"* ]]
    then
        sbatch run_beskow_debug.sh
    else
        sbatch run_beskow.sh
    fi
else
    echo "We're on my pc!"
    ./run_my_pc.sh ${option}
fi