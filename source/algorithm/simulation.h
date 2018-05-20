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
#include "general/class_tic_toc.h"
#include "IO/class_WL_print_data.h"
#include "params/nmspc_WL_constants.h"
#include "algorithm/nmspc_WL_parallelization.h"
#include "nmspc_WL_counters_timers.h"

void do_simulations     (class_worker &);
void wanglandau         (class_worker &);
void check_finish_line  (class_worker &, class_backup &, int &);
void divide_range       (class_worker &, class_backup &);
void print_status       (class_worker &, bool force);

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
