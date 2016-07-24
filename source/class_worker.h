//
// Created by david on 2016-07-24.
//

#ifndef WL_CLASS_WORKER_H
#define WL_CLASS_WORKER_H
#include <mpi.h>
#include <random>
#include <fstream>
#include "class_lattice.h"
#include "constants.h"
#include "randomNumbers.h"

class class_worker {
private:

public:
    class_worker(){
        MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);       //Establish thread number of this worker
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);     //Get total number of threads
        //rn::rng.seed((unsigned long) world_ID);        //Seed the random number generator
    };

    //MPI Communicator
    int world_ID;                                       //Thread number
    int world_size;                                     //Total threads

    //Lattice
    class_lattice lattice;

};


#endif //WL_CLASS_WORKER_H
