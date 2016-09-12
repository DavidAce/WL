#!/bin/bash
echo "Running"
# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
        echo "** Trapped CTRL-C"
}

mpirun -n 6  -bind-to core:overload-allowed bin/Release/WL


