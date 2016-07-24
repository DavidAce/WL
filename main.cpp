#include <iostream>
#include <Eigen/Dense>
#include <mpi.h>
#include "source/class_worker.h"
using namespace std;

int main() {

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Get the number of processes

    // Get the rank of the process
    class_worker worker();
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
                   " out of %d processors\n",
           processor_name, w.world_ID, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}