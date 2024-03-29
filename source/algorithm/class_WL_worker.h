//
// Created by david on 2016-07-24.
//

#ifndef WL_CLASS_WL_WORKER_H
#define WL_CLASS_WL_WORKER_H
#include "general/class_tic_toc.h"
#include "general/nmspc_math_algorithms.h"
#include "general/nmspc_random_numbers.h"
#include "nmspc_WL_counters_timers.h"
#include "params/class_model.h"
#include "params/nmspc_WL_constants.h"
#include <chrono>
#include <Eigen/Core>
#include <exception>
#include <fstream>
#include <iterator>
#include <memory>
#include <mpi.h>
#include <random>
#include <set>
#include <thread>

struct state {
    long E_idx;
    long M_idx;
};

std::ostream &operator<<(std::ostream &out, const std::vector<state> &v);

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
    if(!v.empty()) {
        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::set<T> &v) {
    if(!v.empty()) {
        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}

template<typename Derived, typename T>
typename Derived::PlainObject operator<<(Eigen::ArrayBase<Derived> &in, const std::set<T> &v) {
    if(!v.empty()) {
        Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> temp;
        temp.resize(v.size());
        int k = 0;
        for(auto it = v.begin(); it != v.end(); it++) {
            temp(k) = *it;
            k++;
        }
        in = temp;
    }
    return in;
}

class class_WL_teams;

class class_worker {
    private:
    public:
    class_worker(int &, int &); // Constructor
    // Main data structures of the WL algorithm. Needed very often.
    double lnf;            // Modification factor of WL-algorithm
    int    rate_increment; // Increment rate should be proportional to number of bins

    // WL DOS and Histograms
    Eigen::ArrayXXd    dos;
    Eigen::ArrayXXi    histogram;
    std::vector<state> random_walk;
    Eigen::ArrayXd     E_bins, M_bins;

    // MPI Communicator
    int world_ID;   // Thread number
    int world_size; // Total number of threads

    // Lattice
    class_model model; // Takes argument L, the linear size of your lattice
    // WL Energy and Order parameter and their limits
    double E, M;             // Current Energy and Order parameter
    double E_trial, M_trial; // Proposed
    long   E_idx, M_idx;
    long   E_idx_trial, M_idx_trial;   // Position of trial values
    double E_min_global, M_min_global; // Global minimum
    double E_max_global, M_max_global; // Global maximum
    double E_min_local, M_min_local;   // Local minimum
    double E_max_local, M_max_local;   // Local maximum

    // Sets containing discrete spectrums
    std::set<double> E_set; // Set of found energies, used in discrete do_simulations.
    std::set<double> M_set; // Set of found energies, used in discrete do_simulations.

    // WL acceptance criterion
    bool accept;
    bool state_in_window, trial_in_window, state_is_valid, trial_is_valid;
    int  need_to_resize_global;
    // WL convergence parameters
    int              flag_one_over_t; // turns to 1 when 1/t algorithm starts
    int              finish_line;     // turns to 1 when converged
    std::vector<int> saturation;
    double           slope;

    // Holders for total, merged data
    Eigen::ArrayXXd dos_total;
    Eigen::ArrayXd  E_bins_total, M_bins_total;
    int             iteration;

    // Used for profiling functions in worker
    class_profiling t_total, t_print, t_sweep, t_swap, t_merge, t_setup_team, t_sync_team, t_divide_range, t_check_convergence, t_acceptance_criterion;

    // Used to keep track of your team
    std::shared_ptr<class_WL_teams> team;

    // Functions

    // Startup
    void find_initial_limits();
    void start_counters();
    void rewind_timers();
    void set_initial_local_bins();
    bool __attribute__((always_inline)) check_in_window(const double x) { return x >= E_min_local && x <= E_max_local; }

    // Monte Carlo Sweep
    void        sweep() __attribute__((hot));
    inline void acceptance_criterion() __attribute__((always_inline));
    inline void update_global_range() __attribute__((always_inline));
    void        walk_away_from_window();
    void        walk_towards_window();
    // Convergence
    void                 next_WL_iteration();
    void                 prev_WL_iteration();
    void                 rewind_to_lowest_walk();
    void                 rewind_to_zero();
    void                 set_rate_increment();
    void                 add_hist_volume();
    void                 check_saturation();
    friend std::ostream &operator<<(std::ostream &, const class_worker &);

    template<int on = 0, typename T>
    void debug_print(T input) {
        if constexpr(on == 1) {
            MPI_Barrier(MPI_COMM_WORLD);
            if(world_ID == 0) {
                std::cout << input << std::flush;
                std::this_thread::sleep_for(std::chrono::microseconds(10));
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    template<int on = 0, typename T>
    void debug_print_all(T input) {
        if constexpr(on == 1) {
            MPI_Barrier(MPI_COMM_WORLD);
            for(int i = 0; i < world_size; i++) {
                if(i == world_ID) {
                    std::cout << "ID " << world_ID << ": " << input << std::flush;
                    std::this_thread::sleep_for(std::chrono::microseconds(10));
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
};

class class_backup {
    private:
    bool backed_up;

    public:
    class_backup() { backed_up = false; }
    void backup_state(class_worker &worker) {
        if(!backed_up) {
            lnf            = worker.lnf;
            rate_increment = worker.rate_increment;
            dos            = worker.dos;
            E_bins         = worker.E_bins;
            M_bins         = worker.M_bins;
            E_min_local    = worker.E_min_local;
            E_max_local    = worker.E_max_local;
            M_min_local    = worker.M_min_local;
            M_max_local    = worker.M_max_local;
            E_set          = worker.E_set;
            M_set          = worker.M_set;
            MCS            = counter::MCS;
            walks          = counter::walks;
            backed_up      = true;
            std::cout << "ID: " << worker.world_ID << " Is backed up" << std::endl;
        }
    }

    void restore_state(class_worker &worker) {
        if(backed_up) {
            worker.lnf            = lnf;
            worker.rate_increment = rate_increment;
            worker.dos            = dos;
            worker.E_bins         = E_bins;
            worker.M_bins         = M_bins;
            worker.E_min_local    = E_min_local;
            worker.E_max_local    = E_max_local;
            worker.M_min_local    = M_min_local;
            worker.M_max_local    = M_max_local;
            worker.E_set          = E_set;
            worker.M_set          = M_set;
            counter::MCS          = MCS;
            counter::walks        = walks;
            worker.histogram      = Eigen::ArrayXXi::Zero(dos.rows(), dos.cols());
            worker.saturation.clear();
            worker.random_walk.clear();
            worker.state_is_valid = false;
            backed_up             = false;
            std::cout << "ID: " << worker.world_ID << " Is now restored" << std::endl;
        }
    }

    // Main data structures of the WL algorithm. Needed very often.
    double lnf;            // Modification factor of WL-algorithm
    int    rate_increment; // Increment probability should be proportional to number of bins

    // WL DOS and Histograms
    Eigen::ArrayXXd dos;
    Eigen::ArrayXd  E_bins, M_bins;

    // Lattice
    //    class_model model;
    // WL Energy and Order parameter and their limits
    double E_min_local, M_min_local; // Local minimum
    double E_max_local, M_max_local; // Local maximum
    // Sets containing discrete spectrums
    std::set<double> E_set; // Set of found energies, used in discrete do_simulations.
    std::set<double> M_set; // Set of found energies, used in discrete do_simulations.
    int              MCS;
    int              walks;
};

#endif // WL_CLASS_WL_WORKER_H
