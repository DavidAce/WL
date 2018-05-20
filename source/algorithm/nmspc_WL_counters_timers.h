//
// Created by david on 2016-07-24.
//

#ifndef WL_COUNTERS_TIMERS_H
#define WL_COUNTERS_TIMERS_H

namespace counter {
    //Counters
    extern int MCS;
    extern int walks;
    extern int swaps;
    extern int swap_accepts;
    extern int merges;
};

namespace timer {
    extern int increment;
    extern int add_hist_volume;
    extern int check_finish_line;
    extern int check_saturation;
    extern int backup;
    extern int print;
    extern int swap;
    extern int sync_team;
    extern int setup_team;
    extern int divide_range;
    extern int sampling;
}


#endif //WL_COUNTERS_TIMERS_H
