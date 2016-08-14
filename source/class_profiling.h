//
// Created by david on 2016-08-14.
//

#ifndef WL_CLASS_PROFILING_H
#define WL_CLASS_PROFILING_H

#include <chrono>
#include <iostream>
#include <iomanip>
using namespace std;
using namespace chrono;

class class_profiling {
private:
    high_resolution_clock::time_point delta_tic;
    high_resolution_clock::time_point delta_toc;
    const int profiling;        //Whether we are profiling or not.
public:
    class_profiling(const int &);                 //Constructor
    void tic();
    void toc();
    void reset();
    double total_time;
    duration<double> delta_time;
    friend std::ostream &operator<<(std::ostream &, const class_profiling &);

};


#endif //WL_CLASS_PROFILING_H
