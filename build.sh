#!/bin/bash
echo "Starting Build"
mkdir build
cd build
#rm -rf *

if [[ "$@" == "Debug" ]]
then
	buildtype="Debug"
else
    buildtype="Release"
fi

if [[ "$@" == "intel" ]]
then
	compiler="mpiicpc"
else
    compiler="mpic++"
fi

mkdir ${buildtype}
cd ${buildtype}

cmake -DCMAKE_BUILD_TYPE=${buildtype} -DCMAKE_CXX_COMPILER=${compiler} ../../
make