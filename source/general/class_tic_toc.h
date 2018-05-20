//
// Created by david on 2016-08-14.
//

#ifndef WL_CLASS_TIC_TOC_H
#define WL_CLASS_TIC_TOC_H


#include <chrono>
#include <ratio>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace chrono;
class class_profiling {
private:
    high_resolution_clock::time_point delta_tic;
    high_resolution_clock::time_point delta_toc;
    const int profiling;           //Whether we are profiling or not.
    int print_precision;
    string name;
public:
    class_profiling(const int & on_off, int prec, const string output_text);                 //Constructor
    high_resolution_clock::duration delta_time;
    high_resolution_clock::duration total_time;

    inline void tic(){
        if (profiling) {
            delta_tic = high_resolution_clock::now();
        }
    }
    
    inline void toc(){
        if (profiling) {
            delta_toc   = high_resolution_clock::now();
            delta_time  = delta_toc - delta_tic;
            total_time += delta_time;
        }
    }
    inline void print_delta(){
        if (profiling) {
            cout << setprecision(print_precision)
                 << name
                 << duration_cast<duration<double>>(delta_time).count();
        }
    }
    template<typename rep=double, typename time_type=std::ratio<1,1>>
    inline void print_delta(){
        if (profiling) {
            cout << name
                 << duration_cast<duration<rep,time_type>>(delta_time).count();
        }
    }

    inline void print_total(){
        if (profiling) {
            cout << setprecision(print_precision)
                 << name
                 << duration_cast<duration<double>>(total_time).count();
        }
    }

    template<typename rep=double, typename time_type=std::ratio<1,1>>
    inline void print_total(){
        if (profiling) {
            cout << name
                 << duration_cast<duration<rep,time_type>>(total_time).count();
        }
    }

    inline void print_total_reset(){
        if (profiling) {
            cout << setprecision(print_precision)
                 << name
                 << duration_cast<duration<double>>(total_time).count();
            reset();
        }
    }

    template<typename rep=double, typename time_type=std::ratio<1,1>>
    inline void print_total_reset(){
        if (profiling) {
            cout << name
                 << duration_cast<duration<rep,time_type>>(total_time).count();
            reset();
        }
    }

    void reset();
    friend std::ostream &operator<<(std::ostream &, const class_profiling &);
};


#endif //WL_CLASS_TIC_TOC_H
