//
// Created by david on 2016-08-11.
//

#ifndef CLASS_WL_PARALLELIZATION_H
#define CLASS_WL_PARALLELIZATION_H
#include <thread>
#include <chrono>
#include "general/class_tic_toc.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include "general/nmspc_mpi_extensions.h"
#include "class_WL_worker.h"


namespace parallel {
    extern void swap                           (class_worker &) ;
    extern void swap2                          (class_worker &) ;
    extern void merge                          (class_worker &, bool broadcast, bool setNaN) ;
    extern void broadcast_merger               (class_worker &);
    extern void divide_global_range_uniform    (class_worker &);
    extern void divide_global_range_dos_area   (class_worker &);
    extern void divide_global_range_dos_volume (class_worker &);
    extern void resize_global_range            (class_worker &);
    extern void synchronize_sets               (class_worker &);
    extern void adjust_local_bins              (class_worker &);
//    extern void sync_team                      (class_worker &)__attribute__((hot)) ;
//    extern void setup_team                     (class_worker &);
//    extern void setup_comm                     (class_worker &);
//    template <typename T>
//    void debug_print        (class_worker &worker, T input){
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (worker.world_ID == 0) {
//            cout << input;
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::microseconds(100));
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }




};


#endif //CLASS_WL_PARALLELIZATION_H
