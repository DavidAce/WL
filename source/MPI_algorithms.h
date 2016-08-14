//
// Created by david on 2016-08-11.
//

#ifndef WL_CLASS_MPI_ALGORITHMS_H
#define WL_CLASS_MPI_ALGORITHMS_H
#include <thread>
#include <chrono>
#include "class_profiling.h"
#define debug_swap 0
#define debug_merge 1
#define debug_divide 1

#include "class_worker.h"
namespace mpi {
    extern void swap(class_worker &, class_profiling &t_swap);
    extern void merge(class_worker &);
    extern void merge2(class_worker &);
    extern void merge3(class_worker &);
    extern void divide_global_range_dos_volume(class_worker &);
};


#endif //WL_CLASS_MPI_ALGORITHMS_H
