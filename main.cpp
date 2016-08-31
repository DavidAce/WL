#include <iostream>
#include <mpi.h>
#include "source/simulation.h"
using namespace std;

int main() {

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    class_worker worker;
    do_simulations(worker);

    MPI_Barrier(MPI_COMM_WORLD);

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}