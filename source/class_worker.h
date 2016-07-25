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

#include "class_model.h"
#include "constants.h"
#include "randomNumbers.h"

using namespace Eigen;
class class_worker {
private:

public:
    class_worker();                 //Constructor

    //Functions
    void measure_state();           //Compute current E and M
    void initial_limits();
    void initial_spectrum();
    void update_limits();           //
    void update_Eidx();
    void reset_timers();
    void make_MC_trial();
    void accept_MC_trial();
    void reject_MC_trial();
    //int  find_idx(const double list[], const double &x, const int &idx_min, const int &idx_max);

   // template <typename List_type, typename T>
    //int const& binary_search(const List_type (&)[], const T&, const long int &size);
    template <typename List_type, typename T, typename size_type>
    int binary_search(const List_type &list , const T& x, const size_type &size){
       auto low = std::upper_bound (list, list + size, x); //          ^
       return  (low - list - 1);
    }


    //MPI Communicator
    int world_ID;                   //Thread number
    int world_size;                 //Total threads

    //Lattice
    class_model model;

    //WL DOS and Histograms
    MatrixXi histogram, histogram_temp;
    MatrixXd dos, dos_temp;


    int finish_line;                    //turns to 1 when converged

    double lnf;
    //WL Energy and Order parameter limits
    double E,M;                         //Current Energy and Order parameter
    double E_trial, M_trial;                 //Proposed
    int Eidx, Midx;
    int Eidx_trial, Midx_trial;         //Position of trial values
    double Emin_global, Mmin_global;    //Global minimum
    double Emax_global, Mmax_global;    //Global maximum
    double Emin_local , Mmin_local ;    //Local minimum
    double Emax_local , Mmax_local ;    //Local maximum



    VectorXd E_bins;                     //Energy spectrum
    VectorXd M_bins;                     //Order parameter spectrum
    VectorXd X_bins;                     //Auxiliary spectrum for calculations.
    friend std::ostream &operator<<(std::ostream &, const class_worker &);

    //Counters and flags
    int MCS;
    int timer_append;
    int timer_check;
    int timer_print;
    int timer_swap;
    int timer_resize;

    int flag_foundE;
    int flag_foundM;

};


#endif //WL_CLASS_WORKER_H
