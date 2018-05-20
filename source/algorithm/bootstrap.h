//
// Created by david on 8/31/16.
//

#ifndef WL_BOOTSTRAP_H
#define WL_BOOTSTRAP_H

#include "class_WL_worker.h"
#include "IO/class_WL_print_data.h"
#include "IO/class_WL_read_data.h"
#include "general/nmspc_math_algorithms.h"
#include "nmspc_WL_parallelization.h"
#include "class_WL_thermo.h"

void do_bootstrap(class_worker &worker);
void do_thermodynamics(class_worker &worker);
void do_statistics(class_worker &worker);

#endif //WL_BOOTSTRAP_H
