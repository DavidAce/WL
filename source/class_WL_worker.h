//
// Created by david on 2016-07-24.
//

#ifndef WL_CLASS_WL_WORKER_H
#define WL_CLASS_WL_WORKER_H
#include <mpi.h>
#include <random>
#include <fstream>
#include <thread>
#include <chrono>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <set>
#include <iterator>
#include "class_model.h"
#include "class_tic_toc.h"
#include "nmspc_WL_constants.h"
#include "nmspc_WL_counters_timers.h"
#include "nmspc_math_algorithms.h"
#include "nmspc_random_numbers.h"

using namespace Eigen;

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    if ( !v.empty() ) {
        out << "[ ";
        std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::set<T>& v) {
    if ( !v.empty() ) {
        out << "[ ";
        std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}

template <typename Derived, typename T>
typename Derived::PlainObject operator<< (ArrayBase<Derived> &in , const std::set<T>& v) {
    if (!v.empty()) {
        Array<typename Derived::Scalar,Dynamic, 1> temp;
        temp.resize(v.size());
        int k = 0;
        for (auto it = v.begin(); it != v.end(); it++) {
            temp(k) = *it;
            k++;
        }
        in = temp;
    }
    return in;
}





class class_worker {
private:

public:
    class_worker(int &, int &);                 //Constructor
    //Main data structures of the WL algorithm. Needed very often.
    double   lnf;       //Modification factor of WL-algorithm
    double P_increment;    //Increment probability should be proportional to number of bins

    //WL DOS and Histograms
    ArrayXXd dos;
    ArrayXXi histogram;
    ArrayXd E_bins, M_bins;

    //MPI Communicator
    int world_ID;                   //Thread number
    int world_size;                 //Total number of threads


    //Lattice
    class_model model; //Takes argument L, the linear size of your lattice
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
    bool state_in_window;
    bool state_is_valid;
    int  need_to_resize_global;
    //WL convergence parameters
    int     flag_one_over_t;             //turns to 1 when 1/t algorithm starts
    int     finish_line;                 //turns to 1 when converged
    vector<int> saturation;
    double slope;

    //Holders for total, merged data
    ArrayXXd dos_total;
    ArrayXd E_bins_total, M_bins_total;
    int iteration;

	//Used for profiling functions in worker
    class_profiling t_total                 ,
                    t_print                 ,
                    t_sweep 				,
                    t_swap 					,
                    t_merge 				,
                    t_help_setup            ,
                    t_help                  ,
					t_divide_range          ,
					t_check_convergence 	,
					t_make_MC_trial 		,
					t_acceptance_criterion 	;

    //Used when finished and helping others out
    class helper{
    private:
    public:
        helper(){
            reset();
        }
        void reset(){
            giving_help = false;
            getting_help = false;
            helping_id = -1;
            available  = 0;
            help_walks = 0;
            MPI_COMM_HELP = MPI_COMM_NULL;
        }
        ArrayXXi histogram_recv; //Receive histogram from helpers
        bool giving_help;
        bool getting_help;
        int  helping_id;
        int available;
        int color;
        int key;
        int help_rank;
        int help_size;
        int help_walks;
        MPI_Comm MPI_COMM_HELP;
    };
    helper help;

    //Functions
    void find_current_state();           //Compute current E and M (and their indices)
    void find_next_state_exact() __attribute__((hot));
    void find_next_state() __attribute__((hot));
    void find_next_state(bool) __attribute__((hot));
    void find_initial_limits();
    void start_counters();
    void rewind_timers();

    void set_initial_local_bins();
    void update_global_range();
    void resize_global_range();
    void divide_global_range_uniformly();
    void adjust_local_bins();
    void synchronize_sets();
    void resize_local_bins2();
    void compute_number_of_bins(int &, int &);
//    bool check_in_window(const double) __attribute__((always_inline));
    bool __attribute__((always_inline)) check_in_window(const double x) {
        return x >= E_min_local && x <= E_max_local;
    }
    void make_MC_trial() __attribute__((hot));
    void insert_state(double new_E,double new_M) __attribute__((hot));
    void walk_away_from_window();
    void walk_towards_window();
    void acceptance_criterion() __attribute__((hot));
    void acceptance_criterion2();
    void accept_MC_trial()  __attribute__((hot));
    void reject_MC_trial() __attribute__((hot));
    void set_P_increment();
    void next_WL_iteration();
    void prev_WL_iteration();
    void rewind_to_lowest_walk();
    void rewind_to_zero();
    void add_hist_volume() __attribute__((hot));
    void check_saturation();
    friend std::ostream &operator<<(std::ostream &, const class_worker &);
};


class class_backup{
private:
    bool backed_up;
public:
    class_backup(): model(constants::L){
        backed_up = false;
    }
    void backup_state(class_worker &worker){
        if (!backed_up) {
            lnf             = worker.lnf;
            P_increment     = worker.P_increment;
            dos             = worker.dos;
            histogram       = worker.histogram;
            E_bins          = worker.E_bins;
            M_bins          = worker.M_bins;
            model.lattice   = worker.model.lattice;
            E               = worker.E;
            M               = worker.M;
            E_idx           = worker.E_idx;
            M_idx           = worker.M_idx;
            E_min_local     = worker.E_min_local;
            E_max_local     = worker.E_max_local;
            M_min_local     = worker.M_min_local;
            M_max_local     = worker.M_max_local;
            E_set           = worker.E_set;
            M_set           = worker.M_set;
            in_window       = worker.state_in_window;
            slope           = worker.slope;
            MCS             = counter::MCS;
            walks           = counter::walks;
            backed_up       = true;
            cout << "ID: " << worker.world_ID << " Is backed up" << endl;
        }
    }

    void restore_state(class_worker &worker){
        if (backed_up) {
            worker.lnf          = lnf;
            worker.P_increment  = P_increment;
            worker.dos          = dos;
            worker.histogram    = histogram;
            worker.E_bins       = E_bins;
            worker.M_bins       = M_bins;
            worker.model.lattice= model.lattice;
            worker.E            = E;
            worker.M            = M;
            worker.E_idx        = E_idx;
            worker.M_idx        = M_idx;
            worker.E_min_local  = E_min_local;
            worker.E_max_local  = E_max_local;
            worker.M_min_local  = M_min_local;
            worker.M_max_local  = M_max_local;
            worker.E_set        = E_set;
            worker.M_set        = M_set;
            worker.state_in_window    = in_window;
            worker.slope        = slope;
            counter::MCS        = MCS;
            counter::walks      = walks;
            backed_up = false;
            cout << "ID: " << worker.world_ID << " Is now restored" << endl;
        }
    }

    //Main data structures of the WL algorithm. Needed very often.
    double lnf;             //Modification factor of WL-algorithm
    double P_increment;     //Increment probability should be proportional to number of bins

    //WL DOS and Histograms
    ArrayXXd dos;
    ArrayXXi histogram;
    ArrayXd E_bins, M_bins;

    //Lattice
    class_model model;
    //WL Energy and Order parameter and their limits
    double E,M;                         //Current Energy and Order parameter
    int E_idx, M_idx;
    double E_min_local , M_min_local ;    //Local minimum
    double E_max_local , M_max_local ;    //Local maximum
    //Sets containing discrete spectrums
    std::set<double> E_set;              //Set of found energies, used in discrete do_simulations.
    std::set<double> M_set;              //Set of found energies, used in discrete do_simulations.

    bool in_window;
    double slope;
    int MCS;
    int walks;

};


#endif //WL_CLASS_WL_WORKER_H
