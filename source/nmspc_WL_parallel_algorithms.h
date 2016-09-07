//
// Created by david on 2016-08-11.
//

#ifndef CLASS_WL_PARALLEL_ALGORITHMS_H
#define CLASS_WL_PARALLEL_ALGORITHMS_H
#include <thread>
#include <chrono>
#include "class_tic_toc.h"
#include <Eigen/Dense>
#define debug_swap 0
#define debug_merge 1
#define debug_bcast 1
#define debug_divide 1

#include "class_WL_worker.h"
namespace mpi {
    extern void swap                           (class_worker &) ;
    extern void merge                          (class_worker &) ;
    extern void broadcast                      (class_worker &) ;
    extern void divide_global_range_dos_volume (class_worker &) ;

    template <typename Derived, typename mpitype>
    void bcast_dynamic (ArrayBase<Derived> &arr, mpitype MPI_TYPE , int master_id, int id){
        long int rows = arr.rows();
        long int cols = arr.cols();
//        if (id == master_id){
//            cout << "Master size: " << rows << " x " << cols<<endl;
//        }
        MPI_Bcast(&rows, 1, MPI_LONG_INT, master_id, MPI_COMM_WORLD);
        MPI_Bcast(&cols, 1, MPI_LONG_INT, master_id, MPI_COMM_WORLD);
        if (id != master_id){
            arr.resize(rows, cols);
        }
        MPI_Bcast(arr.derived().data(), (int)(cols*rows), MPI_TYPE, master_id, MPI_COMM_WORLD);
//        cout << "All Size: " << rows << " x " << cols<<endl;

    }
};


#endif //CLASS_WL_PARALLEL_ALGORITHMS_H
