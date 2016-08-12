//
// Created by david on 2016-07-24.
//

#ifndef WL_WANGLANDAU_H
#define WL_WANGLANDAU_H
#include "class_data.h"
#include "constants.h"
#include "math_algorithms.h"
#include "counters_timers.h"
void WangLandau         (class_worker &);
void sweep              (class_worker &); //Function that accepts/rejects new configurations
void check_convergence  (class_worker &, int &);
void add_hist_volume    (class_worker &);
void check_saturation   (class_worker &);
void check_one_over_t   (class_worker &);
void check_global_limits(class_worker &);
void divide_range       (class_worker &);
void print_status       (class_worker &);
void backup_data        (class_worker &, outdata &);
int find_min_positive   (MatrixXi &);
#endif //WL_WANGLANDAU_H
