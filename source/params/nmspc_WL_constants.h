#pragma once
#include <cmath>
#include <spdlog/common.h>
namespace constants {
    /* clang-format off */
    //WL boostrap properties
    inline constexpr int simulation_reps = 1;      //Number of independent do_simulations
    inline constexpr int bootstrap_reps  = 0;      //Number of times to bootstrap-sample existing simulation data

    inline constexpr bool collect_samples     = false;    //Collect energy/magnetization and lattice samples after the simulation has converged.
    inline constexpr int samples_to_collect   = 0;  // How many samples to collect per thread.
    inline constexpr int rate_sampling        = 10; // How often to sample a lattice with respective energies and magnetization

    // Console properties
    inline spdlog::level::level_enum loglevel = spdlog::level::info; // Log level during simulation
    inline bool                      logstamp = false; // prepend time stamp to each log message

    //WL thermodynamics properties
    inline constexpr int    T_num = 500;             //Number of temperatures for thermodynamic quantities
    inline constexpr double T_min = 0.01;            //Minimum temperature for thermodynamic quantities
    inline constexpr double T_max = 6;               //Maximum temperature for thermodynamic quantities

    //Lattice Properties
    inline constexpr int d = 2;                               //Dimension
    inline int L = 10;                                        //Linear size
    inline int N(){return static_cast<int>(std::pow(L,d));}   //Number of spins/particles

    //DOS and Histogram properties
    inline constexpr int rw_dims                   = 1;        //Dimension of random walks (1D or 2D WL)
    inline constexpr int bins                      = 10;       //No lower than 10! (per worker)

    //Rates for checking and printing (MCS units)
    inline constexpr int    rate_swap              = 100;      //How often to swap walkers
    inline constexpr int    rate_add_hist_volume   = 500;      //How often to append reduced volume to an array called "saturation", which indicates if the current walk has converged when it flattens out.
    inline constexpr int    rate_check_finish_line = 5000;     //Check if everybody's modification factor is below minimum_lnf
    inline constexpr int    rate_check_saturation  = 5000;     //How often to check if saturation has flattened out
    inline constexpr int    rate_divide_range      = 5000;     //How often to check if we can merge all dos and split energy subwindows in a smarter way.
    inline constexpr int    rate_sync_team         = 5;
    inline constexpr int    rate_setup_team        = 5000;
    inline constexpr int    rate_print_status      = 1000;    //How often to print in terminal

    //Wang-Landau convergence criteria
    inline double           minimum_lnf            = 1e-6;  // Finish when the natural logarithm of the Wang-Landau modification factor "ln(f)" has reached this value. Small values take longer but give more precise results.
    inline constexpr double check_saturation_from  = 0.9;
    inline constexpr double reduce_factor_lnf      = 0.5;  // Reduce ln(f) by this factor after each histogram saturation
    inline constexpr double overlap_factor_energy  = 0.75; // Controls how much each energy window should overlap with the next. 0 is no overlap, 0.5 , and 1 is fully overlapping.
    inline constexpr double overlap_factor_dos_vol = 1.0;
    inline constexpr double overlap_factor_dos_area= 1.0;

    //Parameter for initial finding global range
//    inline constexpr int max_no_global_change       = 3;

    //Parameters for parallelization
//    extern int           team_size; //See main.cpp!
    inline int               num_teams           = 8; // Default number of random walker teams. In main.cpp this is adjusted to min(num_teams, MPI WORLD SIZE)
    //Parameters for sub-window splitting
    inline constexpr int     max_merges          = 4; // Maximum number of energy window mergers (as the walkers find new energy limits, the windows need to shift)
/* clang-format ON */

}


