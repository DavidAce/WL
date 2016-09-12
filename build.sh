#!/bin/bash
echo "Starting Build"
mkdir build
cd build
#rm -rf *
buildtype="Release"
#CC="mpicc"
#CXX="mpic++"
if [[ "$@" == "Debug" ]]
then
	buildtype="Debug"
else
    buildtype="Release"
fi

#if [[ "$@" == "intel" ]]
#then
#	CC="mpiicc"
#	CXX="mpiicpc"
#else
#    CC="mpicc"
#    CXX="mpic++"
#fi

mkdir ${buildtype}
cd ${buildtype}

#cmake -DCMAKE_BUILD_TYPE=${buildtype} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX}  ../../
cmake -DCMAKE_BUILD_TYPE=${buildtype}  ../../
make