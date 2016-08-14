#!/bin/bash
echo "Starting Build"
mkdir build
cd build
#rm -rf *
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ../../
make



