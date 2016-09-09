//
// Created by david on 2016-07-24.
//

#ifndef WL_CLASS_WL_WORKER_H
#define WL_CLASS_WL_WORKER_H
#include <mpi.h>
#include <random>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <set>
#include "class_model.h"
#include "class_tic_toc.h"
#include "nmspc_WL_constants.h"
#include "nmspc_WL_counters_timers.h"
#include "nmspc_math_algorithms.h"
#include "nmspc_random_numbers.h"
static const int profiling_sweep                =	1;
static const int profiling_swap                 =	0;
static const int profiling_check_global_limits  =	0;
static const int profiling_check_convergence	=   0;
static const int profiling_make_MC_trial 		=	0;
static const int profiling_acceptance_criterion =	0;
static const int debug_comp_numb_bins           =   0;
static const int debug_divide_dos_vol           =   0;
using namespace Eigen;
class class_worker {
private:

public:
    class_worker();                 //Constructor
    //Main data structures of the WL algorithm. Needed very often.
    double   lnf;                         //Modification factor of WL-algorithm
    //WL DOS and Histograms
    ArrayXXd dos;
    ArrayXXi histogram;
    ArrayXd E_bins, M_bins;
    //Model with lattice etc
    class_model lattice;

    //MPI Communicator
    int world_ID;                   //Thread number
    int world_size;                 //Total number of threads

    //Lattice
    class_model model;
    //WL Energy and Order parameter and their limits
    double E,M;                         //Current Energy and Order parameter
    double E_trial, M_trial;                 //Proposed
    int E_idx, M_idx;
    int E_idx_trial, M_idx_trial;         //Position of trial values
    double E_min_global, M_min_global;    //Global minimum
    double E_max_global, M_max_global;    //Global maximum
    double E_min_local , M_min_local ;    //Local minimum
    double E_max_local , M_max_local ;    //Local maximum
    //Sets containing discrete spectrums
    std::set<double> E_set;              //Set of found energies, used in discrete do_simulations.
    std::set<double> M_set;              //Set of found energies, used in discrete do_simulations.


    //WL acceptance criterion
    bool accept;
    bool in_window;
    int  need_to_resize_global;
    int  need_to_resize_local;
    //WL convergence parameters
    int     flag_one_over_t;             //turns to 1 when 1/t algorithm starts
    int     finish_line;                 //turns to 1 when converged

    ArrayXi saturation;                //Measures the histogram saturation

    //Holders for total, merged data
    ArrayXXd dos_total;
    ArrayXd E_bins_total, M_bins_total;

	//Used for profiling functions in worker
    class_profiling t_sweep 				,
					t_swap 					,
					t_check_global_limits 	,
					t_check_convergence 	,
					t_make_MC_trial 		,
					t_acceptance_criterion 	;
    int iteration;
    //Functions
    void find_current_state();           //Compute current E and M (and their indices)
    void find_initial_limits();
    void start_counters();
    void set_initial_local_bins();
    void update_global_range();
    void resize_global_range() __attribute__((hot));
    void divide_global_range();
    void resize_local_bins();
    void resize_local_range();
    void compute_number_of_bins(int &, int &);
    bool check_in_window(const double &);
    void make_MC_trial() __attribute__((hot));
    void acceptance_criterion() __attribute__((hot));
    void accept_MC_trial() __attribute__((hot));
    void reject_MC_trial() __attribute__((hot));
    void next_WL_iteration();
    void prev_WL_iteration();
    void rewind_to_lowest_walk();
    void rewind_to_zero();
    friend std::ostream &operator<<(std::ostream &, const class_worker &);

};


#endif //WL_CLASS_WL_WORKER_H
