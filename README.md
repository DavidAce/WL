# Wang-Landau Algorithm

This program estimates the Density of States (DOS) of discrete models.

## Requirements
* Eigen library. If not found it will be automatically downloaded to `./libs`. It is safe to remove.
* OpenMPI.  `sudo apt install openmpi-bin libopenmpi-dev`
  
## Before install
Get the git submodule EMC-Lib

```
git submodule update --init
```
  
## Usage

Two scripts in the project root folder facilitate usage: `build.sh` and `run.sh`.
To learn more about optional parameters, run these with flag  `-h`.