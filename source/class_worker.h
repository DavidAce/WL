//
// Created by david on 2016-07-24.
//

#ifndef WL_CLASS_WORKER_H
#define WL_CLASS_WORKER_H
#include <mpi.h>
#include <random>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <set>
#include "class_model.h"
#include "constants.h"
#include "counters_timers.h"
#include "math_algorithms.h"
#include "randomNumbers.h"

using namespace Eigen;
class class_worker {
private:

public:
    class_worker();                 //Constructor

    //Functions
    void find_current_state();           //Compute current E and M (and their indices)
    void find_initial_limits();
    void start_counters();
    void set_initial_local_bins();
    void update_global_range();    //
    void resize_global_range();
    void divide_global_range();     //

    void resize_local_bins();       //
    void compute_number_of_bins(int &, int &);
    bool check_in_window(const double &);
    void make_MC_trial();
    void acceptance_criterion();
    void accept_MC_trial();
    void reject_MC_trial();
    void next_WL_iteration();
    void prev_WL_iteration();


    //MPI Communicator
    int world_ID;                   //Thread number
    int world_size;                 //Total number of threads

    //Lattice
    class_model model;

    //WL acceptance criterion
    bool accept;
    bool in_window;
    int  need_to_resize;
    //WL DOS and Histograms
    MatrixXi histogram, histogram_temp;
    MatrixXd dos, dos_temp;
    VectorXi saturation;                //Measures the histogram saturation


    //WL convergence parameters
    int     flag_one_over_t;             //turns to 1 when 1/t algorithm starts
    int     finish_line;                 //turns to 1 when converged
    double  lnf;                         //Modification factor of WL-algorithm

    //WL Energy and Order parameter and their limits
    double E,M;                         //Current Energy and Order parameter
    double E_trial, M_trial;                 //Proposed
    int E_idx, M_idx;
    int E_idx_trial, M_idx_trial;         //Position of trial values
    double E_min_global, M_min_global;    //Global minimum
    double E_max_global, M_max_global;    //Global maximum
    double E_min_local , M_min_local ;    //Local minimum
    double E_max_local , M_max_local ;    //Local maximum



    VectorXd E_bins;                     //Energy spectrum
    VectorXd M_bins;                     //Order parameter spectrum
    VectorXd X_bins;                     //Auxiliary spectrum for calculations.
    VectorXd Y_bins;                     //Auxiliary spectrum for calculations.



    std::set<double> E_set;              //Set of found energies, used in discrete simulations.
    std::set<double> M_set;              //Set of found energies, used in discrete simulations.

    //Holders for total, merged data
    MatrixXd dos_total;
    VectorXd E_bins_total, M_bins_total;

    friend std::ostream &operator<<(std::ostream &, const class_worker &);

};


#endif //WL_CLASS_WORKER_H
