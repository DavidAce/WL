//
// Created by david on 2016-07-24.
//

#ifndef WL_COUNTERS_TIMERS_H
#define WL_COUNTERS_TIMERS_H
#include <chrono>

namespace counter {
    //Counters
    extern int MCS;
    extern int walks;
    extern int swaps;
    extern int swap_accepts;
    extern int merges;
	
};

namespace timer {
    extern int add_hist_volume;
    extern int check_finish_line;
    extern int check_saturation;
    extern int backup;
    extern int print;
    extern int swap;
    extern int take_help;
    extern int setup_help;
    extern int divide_range;


    extern std::chrono::duration<double> elapsed_time_total;
    extern std::chrono::duration<double> elapsed_time_print;
    extern std::chrono::high_resolution_clock::time_point total_tic;
    extern std::chrono::high_resolution_clock::time_point total_toc;
    extern std::chrono::high_resolution_clock::time_point print_tic;
    extern std::chrono::high_resolution_clock::time_point print_toc;
}

namespace profiling {
    extern std::chrono::duration<double> elapsed_time_total;
    extern std::chrono::duration<double> elapsed_time_print;
    extern std::chrono::high_resolution_clock::time_point total_tic;
    extern std::chrono::high_resolution_clock::time_point total_toc;
    extern std::chrono::high_resolution_clock::time_point print_tic;
    extern std::chrono::high_resolution_clock::time_point print_toc;
}

#endif //WL_COUNTERS_TIMERS_H
