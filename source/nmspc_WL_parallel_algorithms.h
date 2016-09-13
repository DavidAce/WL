//
// Created by david on 2016-08-11.
//

#ifndef CLASS_WL_PARALLEL_ALGORITHMS_H
#define CLASS_WL_PARALLEL_ALGORITHMS_H
#include <thread>
#include <chrono>
#include "class_tic_toc.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#define debug_swap 1
#define debug_merge 1
#define debug_bcast 1
#define debug_divide 1

#include "class_WL_worker.h"
//class class_worker;
namespace mpi {
    extern void swap                           (class_worker &) ;
    extern void merge                          (class_worker &) ;
    extern void broadcast                      (class_worker &) ;
    extern void divide_global_range_energy     (class_worker &) ;
    extern void divide_global_range_dos_volume (class_worker &) ;

    template <typename Derived, typename mpitype>
    void bcast_dynamic (ArrayBase<Derived> &arr, mpitype MPI_TYPE , int master_id, int id){
        int rows = (int) arr.rows();
        int cols = (int) arr.cols();
        MPI_Bcast(&rows, 1, MPI_INT, master_id, MPI_COMM_WORLD);
        MPI_Bcast(&cols, 1, MPI_INT, master_id, MPI_COMM_WORLD);
        arr.derived().conservativeResize(rows, cols);
        MPI_Bcast(arr.derived().data(), cols*rows, MPI_TYPE, master_id, MPI_COMM_WORLD);
    }
};


#endif //CLASS_WL_PARALLEL_ALGORITHMS_H
