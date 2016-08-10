//
// Created by david on 2016-07-24.
//

#ifndef WL_CONSTANTS_H
#define WL_CONSTANTS_H
#include <math.h>

namespace constants{

    //Lattice Properties
    static const int d = 2;         //Dimension
    static const int L = 10;        //Linear size
    static const int N = (int)pow(L,d);  //Number of spins/particles

    //DOS and Histogram properties
    static const int rw_dims = 2;      //Dimension of random walks (1D or 2D WL)
    static const int bins = 50;

    //Rates for checking and printing (MCS units)
    static const int    rate_add_hist_volume        = 100;
    static const int    rate_check_finish_line      = 5000;
    static const int    rate_check_saturation       = 2000;
    static const int    rate_resize_global_range    = 1000;
    static const int    rate_swap                   = 500;
    static const int    rate_backup_data            = 50000;
    static const int    rate_print_status           = 10000;

    //Wang-Landau convergence criteria
    static const double minimum_lnf            = 1e-4;
    static const double check_saturation_from  = 0.7;   // From which fraction to check saturation convergence
    static const double reduce_factor_lnf      = 0.5;
    static const double overlap_factor         = 0.25;
}


#endif //WL_CONSTANTS_H
