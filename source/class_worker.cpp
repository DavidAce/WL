//
// Created by david on 2016-07-24.
//
#include <math.h>
#include <fstream>
#include <iomanip>
#include "class_worker.h"
#include "counters_timers.h"
#include "math_algorithms.h"
using namespace std;
using namespace Eigen;
int counter::MCS;
int counter::saturation;
int counter::walks;
int timer::add_hist_volume;
int timer::check_saturation;
int timer::check_finish_line;
int timer::split_windows;
int timer::print;
int timer::swap;
int timer::resize;

//Constructor
class_worker::class_worker(): model(), finish_line(false){
    MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);           //Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);         //Get total number of threads
    rn::rng.seed((unsigned long)world_ID);
    //Initialize matrices
    switch(constants::rw_dims){
        case 1:
            histogram.conservativeResizeLike(MatrixXi::Zero(constants::bins,1));
            dos.conservativeResizeLike(MatrixXd::Zero(constants::bins,1));

            break;
        case 2:
            histogram.conservativeResizeLike(MatrixXi::Zero(constants::bins,constants::bins));
            dos.conservativeResizeLike(MatrixXd::Zero(constants::bins,constants::bins));
            break;
        default:
            std::cout << "Error in class_worker constructor. Wrong dimension for WL-random walk (rw_dims = ?) " << std::endl;
    }
    dos_temp = dos;
    histogram_temp = histogram;
    lnf = 1;

    initial_limits();
    update_local_spectrum();
    measure_current_state();
    start_counters();
    cout << "ID: " << world_ID << " Started OK"<<endl;
}

void class_worker::measure_current_state(){
    E = model.get_E();
    //Only give value to M if we are doing 2D random walks
    switch(constants::rw_dims){
        case 1:
            M = 0;
            break;
        case 2:
            M = model.get_M();
            break;
        default:
            cout << "Error in class_worker::measure_current_state(). Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            exit(1);
    }
    E_idx = math::binary_search(E_bins.data(), E, E_bins.size());
    M_idx = math::binary_search(M_bins.data(), M, M_bins.size());
    in_window = check_in_window(E);
    cout << "ID: " << world_ID << " In window: " << in_window << endl;
    MPI_Barrier(MPI_COMM_WORLD);
}

void class_worker::initial_limits(){
    //Measure, randomize and measure again to get 2 states
    double E1 = model.get_E();
    double M1 = model.get_M();
    double E2 = E1;
    double M2 = M1;

    while (E2 == E1 || M2 == M1){
        model.randomize_lattice();
        E2 = model.get_E();
        M2 = model.get_M();
    }
    E_min_local = fmin(E1,E2);
    E_max_local = fmax(E1,E2);
    M_min_local = fmin(M1,M2);
    M_max_local = fmax(M1,M2);
    //Now compare to the other workers to find global limits
    switch(constants::rw_dims){
        case 1:
            MPI_Allreduce(&E_min_local, &E_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
            MPI_Allreduce(&E_max_local, &E_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
            M_min_global = 0;
            M_max_global = 0;
            break;
        case 2:
            MPI_Allreduce(&E_min_local, &E_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
            MPI_Allreduce(&E_max_local, &E_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
            MPI_Allreduce(&M_min_local, &M_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
            MPI_Allreduce(&M_max_local, &M_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
            break;
        default:
            cout << "Error in class_worker::initial_state(). Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            exit(1);
    }
    update_local_limits();
}

void class_worker::start_counters() {
    counter::MCS = 0;
    counter::saturation = 0;
    counter::walks = 0;

    timer::add_hist_volume      = 0;
    timer::check_saturation     = 0;
    timer::check_finish_line    = 0;
    timer::split_windows        = 0;
    timer::print                = 0;
    timer::swap                 = 0;
    timer::resize               = 0;

    flag_one_over_t             = 0;
}

void class_worker::update_local_spectrum(){
    switch(constants::rw_dims){
        case 1:
            E_bins  = VectorXd::LinSpaced(constants::bins, E_min_local, E_max_local);
            M_bins  = VectorXd::Zero(1);
            break;
        case 2:
            E_bins  = VectorXd::LinSpaced(constants::bins, E_min_local, E_max_local);
            M_bins  = VectorXd::LinSpaced(constants::bins, M_min_local, M_max_local);
            break;
        default:
            cout << "Error in class_worker::update_local_spectrum. Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            exit(1);
    }
}

void class_worker::update_global_limits() {
    E_min_global     = fmin(E_trial, E_min_global);
    E_max_global     = fmax(E_trial, E_max_global);
    M_min_global     = fmin(M_trial, M_min_global);
    M_max_global     = fmax(M_trial, M_max_global);
}


// This function does rebinning of dos and histograms.
// If E_set contains more than the default number of bins, then enlarge E_bins, otherwise shrink it!
// If M_set contains more ----" " ---
void class_worker::update_local_bins() {
    int x, y, i, j, k;
    double weight, weight_sum, dE, dM, dR, dx, dy;
    bool zero = false;
    bool flag = false;
    update_local_limits();

    if (class_model::discrete_model){
        //Purge elements from E_set and M_set outside of window
        for (auto it = E_set.begin(); it != E_set.end(); ) {
            if (!check_in_window(*it)){
                it = E_set.erase(it);
            }
            else {
                ++it;
            }
        }
    }


    //Check if we need more bins
    int E_old_size = (int) E_bins.size();
    int M_old_size = (int) M_bins.size();
    int E_new_size = max(E_old_size, (int) E_set.size());
    int M_new_size = max(M_old_size, (int) M_set.size());

    histogram_temp  = MatrixXi::Zero(E_new_size, M_new_size);
    dos_temp        = MatrixXd::Zero(E_new_size, M_new_size);
    dE              = (E_max_local - E_min_local) / (max(E_new_size,E_old_size) * 2.0);
    dM              = (M_max_local - M_min_local) / (max(M_new_size,E_old_size) * 2.0);
    dR              = sqrt(dE * dE + dM * dM);
    X_bins          = E_bins;
    Y_bins          = M_bins;
    E_bins          = VectorXd::LinSpaced(E_new_size, E_min_local, E_max_local);
    M_bins          = VectorXd::LinSpaced(M_new_size, M_min_local, M_max_local);
    //Coarsen the histogram and dos.
    for (j = 0; j < M_new_size; j++) {
        for (i = 0; i < E_new_size; i++) {
            k = 0;
            weight_sum  = 0;
            zero        = false;
            for (y = 0; y < M_old_size; y++) {
                for (x = 0; x < E_old_size; x++) {
                    if (dos(x, y) == 0) {
                        zero = true;
                        break;
                    }
                    dx = fabs(E_bins(i) - X_bins(x));
                    dy = fabs(M_bins(j) - Y_bins(y));
                    if (dx >= dE) { continue; }
                    if (dy >= dM) { continue; }
                    else {
                        weight                   = fabs(1.0 - sqrt(dx * dx + dy * dy) / dR);
                        dos_temp(i, j)          += dos(x, y) * weight;
                        histogram_temp(i, j)    += histogram(x, y);
                        weight_sum              += weight;
                        k++;
                    }
                }
                if (zero) { break; }
            }
            if (weight_sum > 0 && !zero) {
                dos_temp(i, j)          /= weight_sum;
                histogram_temp(i, j)    /= k;
            } else {
                dos_temp(i, j)          = 0;
                histogram_temp(i, j)    = 0;
            }
        }
    }
    dos             = dos_temp;
    histogram       = histogram_temp;
    E_idx           = math::binary_search(E_bins.data(), E, E_new_size);
    M_idx           = math::binary_search(M_bins.data(), M, M_new_size);
}

void class_worker::update_local_limits() {
    //Update limits
    double global_range     = E_max_global - E_min_global;
    double local_range      = global_range/world_size;
    E_min_local             = fmax(E_min_global,
                                   E_min_global + local_range * (world_ID - 0.5));
    E_max_local             = fmin(E_max_global,
                                   E_max_global - local_range * (world_size - (world_ID + 1) - 0.5));
}

void class_worker::split_global_spectrum(){
    //Compare to the other workers to find global limits
    if (timer::split_windows >= constants::rate_split_windows) {
        timer::split_windows = 0;
        //merge_windows(worker);
        switch (constants::rw_dims) {
            case 1:
                MPI_Allreduce(&E_min_global, &E_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
                MPI_Allreduce(&E_max_global, &E_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                M_min_global = 0;
                M_max_global = 0;
                break;
            case 2:
                MPI_Allreduce(&E_min_global, &E_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
                MPI_Allreduce(&E_max_global, &E_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                MPI_Allreduce(&M_min_global, &M_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
                MPI_Allreduce(&M_max_global, &M_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                break;
            default:
                cout << "Error in check_windows(). Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
                exit(1);
        }
        update_local_limits();
        update_local_bins();


    }else{
        timer::split_windows++;
    }
}


bool class_worker::check_in_window(const double &x) {
    return x >= E_min_local && x <= E_max_local;
}

void class_worker::make_MC_trial(){
    model.make_new_state(E,M, E_trial, M_trial);
    update_global_limits();
}

void class_worker::acceptance_criterion(){
    //Scenarios:
    //  1) in_window = true, E_trial is in_window ->  Find idx -> MC-test
    //  2) in_window = true, E_trial not in_window->  Reject
    //  3) in_window = false,E_trial is in_window ->  Accept -> Find idx -> set in_window true
    //  4) in_window = false,E_trial not in_window->  Accept (without updating dos)

    if (in_window){
        E_set.insert(E);
        M_set.insert(M);
        E_bins(E_idx) = E;
        M_bins(M_idx) = M;
        if (check_in_window(E_trial)){
            E_idx_trial = math::binary_search(E_bins.data(), E_trial, E_bins.size());
            M_idx_trial = math::binary_search(M_bins.data(), M_trial, M_bins.size());
            accept      = rn::uniform_double(0,1) < exp(dos(E_idx, M_idx) - dos(E_idx_trial, M_idx_trial));
        }else{
            accept      = false;
        }
    }else{
        if (check_in_window(E_trial)){
            accept      = true;
            in_window   = true;
            E_idx_trial = math::binary_search(E_bins.data(), E_trial, E_bins.size());
            M_idx_trial = math::binary_search(M_bins.data(), M_trial, M_bins.size());
        }else{
            accept      = true;
        }
    }
}

void class_worker::accept_MC_trial() {
    E = E_trial;
    M = M_trial;
    E_idx = E_idx_trial;
    M_idx = M_idx_trial;
    model.flip();
    if (in_window) {
        dos(E_idx, M_idx) += lnf;
        histogram(E_idx, M_idx) += 1;
    }
}

void class_worker::reject_MC_trial() {
    if (in_window) {
        dos(E_idx, M_idx) += lnf;
        histogram(E_idx, M_idx) += 1;
    }
}

void class_worker::next_WL_iteration() {
    lnf *= constants::reduce_factor_lnf;
    histogram.fill(0);
    saturation.fill(0);
    counter::saturation = 0;
    counter::walks++;
}


//Function for printing the lattice. Very easy if L is an Eigen matrix.
std::ostream &operator<<(std::ostream &os, const class_worker &worker) {
    os << "ID: " << worker.world_ID
       << " Current State: " << endl
       << "     E = " << worker.E << " (" << worker.E_idx << ")"
       << "     M = " << worker.M << " (" << worker.M_idx << ")" << endl
       << "     E_trial = " << worker.E_trial << " (" << worker.E_idx_trial << ")"
       << "     M_trial = " << worker.M_trial << " (" << worker.M_idx_trial << ")" << endl
       << "     Local E "
       << "["
       << worker.E_min_local << " "
       << worker.E_max_local << "]"
       << " M "
       << "["
       << worker.M_min_local << " "
       << worker.M_max_local << "]" << endl
       << "     Global E "
       << "["
       << worker.E_min_global << " "
       << worker.E_max_global << "]"
       << " M "
       << "["
       << worker.M_min_global << " "
       << worker.M_max_global << "]"
       << std::endl
       << "     E_bins: " << worker.E_bins.transpose() << endl
       << "     M_bins: " << worker.M_bins.transpose() << endl
       << "     E_set : ";
        for (auto it = worker.E_set.begin(); it != worker.E_set.end(); it++) {
            if (*it == worker.E){ os << "[";}
            os << *it ;
            if (*it == worker.E){ os << "]";}
            os << " ";

        }
        os << endl
        << "     M_set : ";
        for (auto it = worker.M_set.begin(); it != worker.M_set.end(); it++) {
            if (*it == worker.M){ os << "[";}
            os << *it;
            if (*it == worker.M){ os << "]";}
            os << " ";

        }
        os << endl;
    return os;
}



