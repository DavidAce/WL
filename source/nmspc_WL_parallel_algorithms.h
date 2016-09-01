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
};


#endif //CLASS_WL_PARALLEL_ALGORITHMS_H
