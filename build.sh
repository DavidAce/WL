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

#if [[ "$@" == "intel" ]]
#then
#	CCcompiler="mpiicc"
#	CXXcompiler="mpiicpc"
#else
#    CCcompiler="mpicc"
#    CXXcompiler="mpic++"
#fi
CCcompiler="mpiicc"
CXXcompiler="mpiicpc"
mkdir ${buildtype}
cd ${buildtype}

cmake -DCMAKE_BUILD_TYPE=${buildtype}  ../../
make