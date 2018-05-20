//
// Created by david on 2016-08-14.
//

#include "class_tic_toc.h"

class_profiling::class_profiling(const int & on_off, int prec, const string output_text)
        : profiling(on_off),
          print_precision(prec),
          name(output_text)
{
    if (profiling) {
        if (!name.empty()){
            name = " " + name + ": ";
        }
        total_time = total_time.zero();
        delta_time = delta_time.zero();
    }
}

void class_profiling::reset() {
    if (profiling) {
        total_time = total_time.zero();
    }
}


std::ostream &operator<<(std::ostream &os, const class_profiling &t) {
    if (t.profiling) {
        os  << setprecision(5)
            << t.name
            << t.total_time.count();
    }
    return os;
}

