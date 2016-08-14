//
// Created by david on 2016-08-14.
//

#include "class_profiling.h"

class_profiling::class_profiling(const int & p): profiling(p) {
    if (profiling) {
        total_time = 0;
    }
}

void class_profiling::tic() {
    if (profiling) {
        delta_tic = high_resolution_clock::now();
    }
}

void class_profiling::toc() {
    if (profiling) {
        delta_toc   = high_resolution_clock::now();
        delta_time  = delta_toc - delta_tic;
        total_time += delta_time.count();
    }
}

void class_profiling::reset() {
    if (profiling) {
        total_time = 0;
    }
}


//Function for printing the lattice. Very easy if L is an Eigen matrix.
std::ostream &operator<<(std::ostream &os, const class_profiling &t) {
    if (t.profiling) {
        os << setprecision(5);
        os << t.total_time;
    }
    return os;
}