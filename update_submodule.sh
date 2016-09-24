#!/bin/bash
echo "Updating submodule EMC-Lib"
cd EMC-Lib
git pull
cd ..
git commit -am 'Updated EMC-Lib'
echo "Finished updating"
