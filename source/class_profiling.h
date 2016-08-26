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
    inline void tic(){
        if (profiling) {
            delta_tic = high_resolution_clock::now();
        }
    }
    
    inline void toc(){
        if (profiling) {
            delta_toc   = high_resolution_clock::now();
            delta_time  = delta_toc - delta_tic;
            total_time += delta_time.count();
        }
    }
    void reset();
    double total_time;
    duration<double> delta_time;
    friend std::ostream &operator<<(std::ostream &, const class_profiling &);

};


#endif //WL_CLASS_PROFILING_H
