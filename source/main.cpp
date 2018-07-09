#include <iostream>
#include <mpi.h>
#include <params/nmspc_WL_constants.h>
#include <algorithm/class_WL_worker.h>
#include "algorithm/simulation.h"
#include "algorithm/bootstrap.h"


using namespace std;

int main() {

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    int world_ID,world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);           //Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);         //Get total number of threads
    constants::num_teams = min(world_size , constants::num_teams);
//    constants::team_size = world_size / constants::num_teams;
    class_worker worker(world_ID, world_size);
    do_simulations(worker);
    do_bootstrap  (worker);
    do_sampling   (worker);
    if(world_ID == 0){cout << "Finished successfully" << endl;}

    MPI_Finalize();
    return 0;
}