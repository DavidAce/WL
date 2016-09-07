//
// Created by david on 8/31/16.
//

#ifndef WL_BOOTSTRAP_H
#define WL_BOOTSTRAP_H

#include "class_WL_worker.h"
#include "class_WL_print_data.h"
#include "class_WL_read_data.h"
#include "nmspc_math_algorithms.h"
#include "nmspc_WL_parallel_algorithms.h"
#include "class_WL_thermo.h"
static const int debug_boot                 =	1;
static const int debug_thermo                =	1;
static const int debug_stats                 =	1;
void do_bootstrap(class_worker &worker);
void do_thermodynamics(class_worker &worker);
void do_statistics(class_worker &worker);

#endif //WL_BOOTSTRAP_H
