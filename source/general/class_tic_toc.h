//
// Created by david on 2016-08-14.
//

#ifndef WL_CLASS_TIC_TOC_H
#define WL_CLASS_TIC_TOC_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <ratio>

class class_profiling {
    private:
    std::chrono::high_resolution_clock::time_point delta_tic;
    std::chrono::high_resolution_clock::time_point delta_toc;
    const int                                      profiling; // Whether we are profiling or not.
    int                                            print_precision;
    std::string                                    name;

    public:
    class_profiling(const int &on_off, int prec, const std::string &output_text); // Constructor
    std::chrono::high_resolution_clock::duration delta_time;
    std::chrono::high_resolution_clock::duration total_time;

    inline void tic() {
        if(profiling) { delta_tic = std::chrono::high_resolution_clock::now(); }
    }

    inline void toc() {
        if(profiling) {
            delta_toc  = std::chrono::high_resolution_clock::now();
            delta_time = delta_toc - delta_tic;
            total_time += delta_time;
        }
    }
    inline void print_delta() {
        if(profiling) {
            std::cout << std::setprecision(print_precision) << name << std::chrono::duration_cast<std::chrono::duration<double>>(delta_time).count();
        }
    }
    template<typename rep = double, typename time_type = std::ratio<1, 1>>
    inline void print_delta() {
        if(profiling) { std::cout << name << std::chrono::duration_cast<std::chrono::duration<rep, time_type>>(delta_time).count(); }
    }

    inline void print_total() {
        if(profiling) {
            std::cout << std::setprecision(print_precision) << name << std::chrono::duration_cast<std::chrono::duration<double>>(total_time).count();
        }
    }

    template<typename rep = double, typename time_type = std::ratio<1, 1>>
    inline void print_total() {
        if(profiling) { std::cout << name << std::chrono::duration_cast<std::chrono::duration<rep, time_type>>(total_time).count(); }
    }

    inline void print_total_reset() {
        if(profiling) {
            std::cout << std::setprecision(print_precision) << name << std::chrono::duration_cast<std::chrono::duration<double>>(total_time).count();
            reset();
        }
    }

    template<typename rep = double, typename time_type = std::ratio<1, 1>>
    inline void print_total_reset() {
        if(profiling) {
            std::cout << name << std::chrono::duration_cast<std::chrono::duration<rep, time_type>>(total_time).count();
            reset();
        }
    }

    void                 reset();
    friend std::ostream &operator<<(std::ostream &, const class_profiling &);
};

#endif // WL_CLASS_TIC_TOC_H
