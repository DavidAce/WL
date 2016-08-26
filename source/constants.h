//
// Created by david on 2016-07-24.
//

#ifndef WL_CONSTANTS_H
#define WL_CONSTANTS_H
#include <math.h>

namespace constants{

    //Lattice Properties
    static const int d = 2;         //Dimension
    static const int L = 16;        //Linear size
    static const int N = (int)pow(L,d);  //Number of spins/particles

    //DOS and Histogram properties
    static const int rw_dims    = 2;      //Dimension of random walks (1D or 2D WL)
    static const int bins       = 3;      //No lower than 3!

    //Rates for checking and printing (MCS units)
    static const int    rate_add_hist_volume   = 500;
    static const int    rate_check_finish_line = 5000;
    static const int    rate_check_saturation  = 5000;
    static const int    rate_check_limits      = 500;
    static const int    rate_split_windows     = 50000;
    static const int    rate_swap              = 500;
    static const int    rate_backup_data       = 500000;
    static const int    rate_print_status      = 10000;

    //Wang-Landau convergence criteria
    static const double minimum_lnf            = 1e-6;
    static const double check_saturation_from  = 0.95;
//    static const double reduce_factor_lnf      = 1-exp(-1.0); // 656 s (check from 0.75)
//    static const double reduce_factor_lnf      = 1-exp(-1.0); //  641 s (check from 0.75)
    static const double reduce_factor_lnf      = 0.5;           // 131 s (check from 0.9
//    static const double reduce_factor_lnf      = 0.5;           // 133 s (check from 0.95
//    static const double reduce_factor_lnf      = 0.25;           //  over 177 s
    static const double overlap_factor         = 0.25;
    static const double one_over_t_factor      = 1.0;
    static const double one_over_t_exponent    = 1.0;

    //Parameters for sub-window splitting
    static const int     min_walks             = 5;
    static const int     max_merges            = 2;
}


#endif //WL_CONSTANTS_H
