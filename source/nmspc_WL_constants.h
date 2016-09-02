//
// Created by david on 2016-07-24.
//

#ifndef WL_NMSPC_WL_CONSTANTS_H
#define WL_NMSPC_WL_CONSTANTS_H
#include <math.h>

namespace constants{

    //WL boostrap properties
    static const int simulation_reps = 2;      //Number of independent do_simulations
    static const int bootstrap_reps  = 2;

    //WL thermodynamics properties
    static const int    T_num = 500;             //Number of temperatures for thermodynamic quantities
    static const double T_min = 0;               //Minimum temperature for thermodynamic quantities
    static const double T_max = 6;               //Maximum temperature for thermodynamic quantities

    //Lattice Properties
    static const int d = 2;         //Dimension
    static const int L = 4;        //Linear size
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
    static const double minimum_lnf            = 1e-4;
    static const double check_saturation_from  = 0.9;
    static const double reduce_factor_lnf      = 0.5;           // 131 s (check from 0.9
    static const double overlap_factor         = 0.5;
    static const double one_over_t_factor      = 1.0;
    static const double one_over_t_exponent    = 1.0;

    //Parameters for sub-window splitting
    static const int     min_walks             = 5;
    static const int     max_merges            = 1;
}


#endif //WL_NMSPC_WL_CONSTANTS_H
