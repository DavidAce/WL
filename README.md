# Wang-Landau Algorithm

The Wang-Landau algorithm estimates the Density of States (DOS) of discrete models as function of energy and (optionally) magnetization,
often denoted `g(E)` (or `g(E,M)`) in physics litterature.

The DOS can give us all the thermodynamic averages by differentiating the partition function `Z= Σ_E g(E,M)exp(-βE)`.

By default, this program will generate the DOS for a 10x10 2D Ising model. 

### In short
Independent "random walkers" build up the DOS by sampling the phase space while keeping a histogram of visited energies. The probability of 
taking the next step in phase space depends on the current value of the DOS and the value at the next point. 

At each step of the random walk the DOS at the current energy and magnetization is incremented
by a modification factor `f`, which starts out as `f=exp(1)`. A stage of the simulation runs until the convergence condition is met, which is when the random walker visits
energies with uniform probability, meaning that the histogram is roughly flat. When that happens, the next stage starts with a smaller modification factor,
`f=exp(0.5)` and runs until the convergence condition is met again. 

Reducing the modification factor throughout many stages allows the DOS to become
smoother and smoother. The simulation ends when the modification factor reaches a very small value, such as `f = f_min = 1e-7`.

 
### Parallelization 
The algorithm is parallelized by splitting the DOS into overlapping energy windows, and each MPI is assigned a window. The processes perform 
independent random walks within their own window, but can with a small probability swap windows with other threads. This is to avoid threads becoming
stuck for too long in one place, which would ruin the smoothness of the DOS.

Some energy windows may reach `f_min` quicker than others, in which case the "finished" processes cease their random walks and go
join another window, to help that part of the DOS converge quicker. 

## Minimum Requirements

The following software is required to build the project:

- C++17 compiler. Tested with:
    * *To build dependencies*: Fortran compiler
- CMake version >= 3.20.

## Dependencies

- [**Eigen**](http://eigen.tuxfamily.org) for tensor and matrix and linear algebra (tested with version >= 3.3).
- [**h5pp**](https://github.com/DavidAce/h5pp) a wrapper for HDF5.
- [**fmt**](https://github.com/fmtlib/fmt) for formatting strings (bundled with h5pp).
- [**spdlog**](https://github.com/gabime/spdlog) for logging (bundled with h5pp).
- [**CLI11**](https://github.com/CLIUtils/CLI11) For parsing input arguments (WIP).
- [**OpenMPI**](https://www.open-mpi.org/) For parallelization. On Ubuntu/Debian systems, install with  `sudo apt install openmpi-bin libopenmpi-dev` 


## Build
The simplest way to build is to use a CMake preset. Some common presets are found in `CMakePresets.json`, 
which can serve as a starting point for generating your own in `CMakeUserPresets.json`. 
[Read more](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html).

To list available presets, run

```
cmake --list-presets
```

Selecting a preset, for example `release-gcc-11-native-cmake`, which sets the `-march=native` compiler flag and uses
`WL_PACKAGE_MANAGER=cmake`, to use CMake exclusively to build the dependency tree. To configure and build, run

```
cmake --preset=release-gcc-11-native-cmake
cmake --build --preset=release-gcc-11-native-cmake
```

Use the CMake argument `WL_PACKAGE_MANAGER=<find|cmake|cpm|conan>` to choose how WL handles dependencies. 
The option `find` (default) will simply call `find_package(...)` for all dependencies, and leave their installation
up to the user. 

The CMake flag `WL_PACKAGE_MANAGER` controls the automated behavior for finding or installing dependencies. It can
take one of these string values:

| Value                | Description                                                                                                                                         |
|----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
| `find` **(default)** | Use CMake's `find_package` to find dependencies                                                                                                     |
| `cmake` **¹**        | Use isolated CMake instances to install dependencies during configure. Disregards pre-installed dependencies on your system             |
| `cpm`                | Use [CPM](https://github.com/cpm-cmake/CPM.cmake) to install dependencies. Disregards pre-installed dependencies on your system         |
| `conan` **²**        | Use the [Conan package manager](https://conan.io/) to install dependencies automatically. Disregards libraries elsewhere on your system |

There are several variables you can pass to CMake to guide `find_package` calls and install location,
see [CMake options](#cmake-options) below.

**¹** Dependencies are installed into `${WL_PKG_INSTALL_DIR}`.

**²** Conan is guided by `conanfile.txt` found in this project's root directory. This method requires conan to be
installed prior (for instance through `pip`, `conda`, `apt`, etc). To let CMake find conan you have three options:
* 
* Add Conan install (or bin) directory to the environment variable `PATH`.
* Export Conan install (or bin) directory in the environment variable `CONAN_PREFIX`, i.e. from command
  line: `export CONAN_PREFIX=<path-to-conan>`
* Give the variable `CONAN_PREFIX` directly to CMake, i.e. from command
  line: `cmake -DCONAN_PREFIX:PATH=<path-to-conan> ...`

## Run
The scripts `run*.sh` show example usage. By default `./run.sh` will run a simulation with `mpirun`
using 4 cores.

All parameters are defined at compile time, by modifying them in the header file `source/params/nmspc_WL_constants.h`.
Please refer to the comments in that file to learn more.


## Output data

All output data is placed under the `outdata/` folder as text files.
 
In subfolders numbered `outdata/0/` to `outdata/N-1/`, where `N` is the number of MPI threads that were run. 
These subfolders contain their corresponding part of the DOS and some thermodynamic data
used for the final averaged results.  

The folder `outdata/final/`  contains the main results, in the form of thermodynamic averages and corresponding error
estimates from made from bootstrapping portions of the DOS. For the bootstrap to be meaningful, you need to run multiple simulations. You can
control the number of independent simulations with the parameter `simulation_reps` found in the header file `source/params/nmspc_WL_constants.h`.
 
Enable `constants::collect_samples` to obtain a subfolder `outdata/samples/` containing lattice samples that are uniformly
distributed in energy and magnetization. This is useful for generating training data for machine learning algorithms, for instance.


