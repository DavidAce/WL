//
// Created by david on 2016-07-24.
//

#ifndef WL_SIMULATION_H
#define WL_SIMULATION_H
#include <iomanip>
#include <thread>
#include <chrono>
#include "class_tic_toc.h"
#include "class_WL_print_data.h"
#include "nmspc_math_algorithms.h"
#include "nmspc_WL_constants.h"
#include "nmspc_WL_parallel_algorithms.h"
#include "nmspc_WL_counters_timers.h"

void do_simulations     (class_worker &);
void wanglandau         (class_worker &);
void sweep              (class_worker &);
void help_out           (class_worker &, class_worker &);
void check_convergence  (class_worker &, int &);
void add_hist_volume    (class_worker &);
void check_saturation   (class_worker &);
void check_one_over_t   (class_worker &);
void check_global_limits(class_worker &);
void divide_range       (class_worker &);
void debug_print        (class_worker &, string);
void print_status       (class_worker &);
void backup_data        (class_worker &, outdata &);
#endif //WL_SIMULATION_H
