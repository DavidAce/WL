#!/bin/bash
echo "Starting Build"
mkdir build
cd build
#rm -rf *
buildtype="Release"

if [[ "$@" == *"ebug"* ]]
then
	buildtype="Debug"
fi

if [[ "$@" == *"lean"* ]]
then
	buildtype="Clean"
fi


if [[ "$HOSTNAME" == *"triolith"* ]]
then
    echo "We're on triolith!";
    export CC=/software/apps/gcc/5.3.0/build01/bin/gcc
    export CXX=/software/apps/gcc/5.3.0/build01/bin/g++
else
    echo "We're on my pc!"
fi


mkdir ${buildtype}
cd ${buildtype}

cmake -DCMAKE_BUILD_TYPE=${buildtype}  ../../
make