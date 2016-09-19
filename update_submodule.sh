#!/bin/bash
echo "Updating submodule EMC-Lib"
cd EMC-Lib
git pull origin master
cd ..
git commit -am 'Updated EMC-Lib'
echo "Finished updating"
