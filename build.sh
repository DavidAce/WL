#!/bin/bash
echo "Starting Build"
mkdir build
cd build
#rm -rf *
buildtype="Release"
if [[ $# -eq 1 ]] ; then
	echo "Argument supplied: " $1
	buildtype=$1
fi
mkdir ${buildtype}
cd ${buildtype}
cmake -DCMAKE_BUILD_TYPE=${buildtype} ../../
make
