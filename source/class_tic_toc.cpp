//
// Created by david on 2016-08-14.
//

#include "class_tic_toc.h"

class_profiling::class_profiling(const int & p): profiling(p) {
    if (profiling) {
        total_time = 0;
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
