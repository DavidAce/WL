#include <iostream>
#include <mpi.h>
#include "source/simulation.h"
#include "source/bootstrap.h"


using namespace std;

int main() {

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    int world_ID,world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);           //Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);         //Get total number of threads
    class_worker worker(world_ID, world_size);
    do_simulations(worker);
    do_bootstrap  (worker);
    MPI_Finalize();
    return 0;
}