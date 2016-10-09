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
    extern int area_merges;
    extern int vol_merges;
    extern int no_global_change;

};

namespace timer {
    extern int increment;
    extern int add_dos;
    extern int add_hist_volume;
    extern int check_finish_line;
    extern int check_saturation;
    extern int backup;
    extern int print;
    extern int swap;
    extern int check_help;
    extern int take_help;
    extern int setup_help;
    extern int divide_range;
    extern int divide_range_find;

}


#endif //WL_COUNTERS_TIMERS_H
