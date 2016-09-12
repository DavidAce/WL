#!/bin/bash
echo "Starting Build"
mkdir build
cd build
#rm -rf *
buildtype="Release"

if [[ "$@" == "Debug" ]]
then
	buildtype="Debug"
fi

if [[ "$@" == "Clean" ]]
then
	buildtype="Clean"
fi

mkdir ${buildtype}
cd ${buildtype}

cmake -DCMAKE_BUILD_TYPE=${buildtype}  ../../
make