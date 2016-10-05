//
// Created by david on 2016-07-24.
//

#ifndef WL_SIMULATION_H
#define WL_SIMULATION_H
#include <iomanip>
#include <thread>
#include <chrono>
#include <ratio>
#include <numeric>
#include "class_tic_toc.h"
#include "class_WL_print_data.h"
#include "nmspc_math_algorithms.h"
#include "nmspc_WL_constants.h"
#include "nmspc_WL_parallel_algorithms.h"
#include "nmspc_WL_counters_timers.h"

void do_simulations     (class_worker &);
void find_global_range  (class_worker &);
void wanglandau         (class_worker &);
void sweep              (class_worker &)                            __attribute__((hot));
void check_finish_line(class_worker &, outdata &out, int &)       __attribute__((hot));
void divide_range_find  (class_worker &)                            __attribute__((hot));
void divide_range       (class_worker &, class_backup &, outdata &) __attribute__((hot));
void backup_to_file(class_worker &, outdata &);
void print_status       (class_worker &, bool force)                __attribute__((hot));


template <typename T>
void debug_print        (class_worker &worker, T input){
    MPI_Barrier(MPI_COMM_WORLD);
    if (worker.world_ID == 0) {
        cout << input;
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(10));
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


#endif //WL_SIMULATION_H
