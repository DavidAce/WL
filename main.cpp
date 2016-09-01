#include <iostream>
#include <mpi.h>
#include "source/simulation.h"
#include "source/bootstrap.h"
using namespace std;

int main() {

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    class_worker worker;
    do_simulations(worker);
    do_bootstrap  (worker);
    MPI_Finalize();
    return 0;
}