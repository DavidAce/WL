#!/bin/bash
echo "Running"
# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

if [[ "$@" == *"valgrind"* ]]
then
    valgrind --tool=memcheck --leak-check=full -v mpirun -n 4 ./build/Debug/WL
elif [[ "$@" == *"gdb"* ]]
then
    mpirun -n 4 xterm gdb build/Debug/WL
#    mpirun -n 4 gdb build/Debug/WL
else
    mpirun -n 4  -bind-to core:overload-allowed build/Release/WL
fi
function ctrl_c() {
        echo "** Trapped CTRL-C"
}



