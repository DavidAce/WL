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



## Requirements
* Eigen library. If not found it will be automatically downloaded to `./libs`. It is safe to remove.
* OpenMPI.  `sudo apt install openmpi-bin libopenmpi-dev`
  
  
## Usage

Two scripts in the project root folder facilitate usage: `build.sh` and `run.sh`.
To learn more about optional parameters, run these with flag  `-h`.
By default, `./build.sh` will build in release mode, and `./run.sh` will run a simulation with `mpirun`
using 4 cores.
 
The default setting will run 2 independent and consecutive simulations, to generate 
good quality thermodynamic averages using bootstrap. See below on how to control the number of simulations.

## Run time parameters

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


