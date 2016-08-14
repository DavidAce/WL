//
// Created by david on 2016-07-24.
//

#ifndef WL_CONSTANTS_H
#define WL_CONSTANTS_H
#include <math.h>

namespace constants{

    //Lattice Properties
    static const int d = 2;         //Dimension
    static const int L = 32;        //Linear size
    static const int N = (int)pow(L,d);  //Number of spins/particles

    //DOS and Histogram properties
    static const int rw_dims    = 2;      //Dimension of random walks (1D or 2D WL)
    static const int bins       = 3;      //No lower than 3!

    //Rates for checking and printing (MCS units)
    static const int    rate_add_hist_volume   = 500;
    static const int    rate_check_finish_line = 5000;
    static const int    rate_check_saturation  = 5000;
    static const int    rate_check_limits      = 1000;
    static const int    rate_split_windows     = 50000;
    static const int    rate_swap              = 1000;
    static const int    rate_backup_data       = 500000;
    static const int    rate_print_status      = 10000;

    //Wang-Landau convergence criteria
    static const double minimum_lnf            = 1e-6;
    static const double check_saturation_from  = 0.9;   // From which fraction to check saturation convergence
    static const double reduce_factor_lnf      = 0.5;
    static const double overlap_factor         = 0.2;
    static const double one_over_t_factor      = 1.0;
    static const double one_over_t_exponent    = 1.0;

    //Parameters for sub-window splitting
    static const int     min_walks             = 5;
    static const int     max_merges            = 2;
}


#endif //WL_CONSTANTS_H
