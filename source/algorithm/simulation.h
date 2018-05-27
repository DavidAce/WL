//
// Created by david on 2016-07-24.
//

#ifndef WL_SIMULATION_H
#define WL_SIMULATION_H
#include <mpi.h>
class class_worker;
class class_backup;


void do_simulations     (class_worker &);
void do_sampling        (class_worker &);
void wanglandau         (class_worker &);
void check_finish_line  (class_worker &, class_backup &, int &);
void divide_range       (class_worker &, class_backup &);
void sync_teams         (class_worker &);
void setup_teams        (class_worker &);
void print_status       (class_worker &, bool force);




#endif //WL_SIMULATION_H
