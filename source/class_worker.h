//
// Created by david on 2016-07-24.
//

#ifndef WL_CLASS_WORKER_H
#define WL_CLASS_WORKER_H
#include <mpi.h>
#include <random>
#include "class_lattice.h"
#include "constants.h"


class class_worker {
private:
    typedef std::mt19937 RNGType;
    RNGType rng;
public:
    class_worker(){
        MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        rng.seed(world_ID);
    };
    int world_ID;
    int world_size;
    class_lattice lattice(constants::L);

    //Random functions
    double  uniform_double(RNGType *, const double &, const double &);
    int     uniform_integer(RNGType *, const int &, const int &);
    double  gaussian_truncated(RNGType *, const double &, const double &, const double &, const double &);


};


#endif //WL_CLASS_WORKER_H
