//
// Created by david on 2016-07-24.
//
#include <math.h>
#include <fstream>
#include <iomanip>
#include "class_worker.h"
#include "counters_timers.h"
#include "math_algorithms.h"
#include <thread>
#include <chrono>
#define debug_divide 1
#define debug_resize_local_range 1
using namespace std;
using namespace Eigen;
int counter::MCS;
int counter::saturation;
int counter::walks;
int counter::swaps;
int counter::swap_accepts;
int counter::merges;
int timer::add_hist_volume;
int timer::check_saturation;
int timer::check_finish_line;
int timer::check_limits;
int timer::backup;
int timer::print;
int timer::swap;
int timer::split_windows;
std::chrono::duration<double> timer::elapsed_time_total;
std::chrono::duration<double> timer::elapsed_time_print;
std::chrono::high_resolution_clock::time_point timer::total_tic;
std::chrono::high_resolution_clock::time_point timer::total_toc;
std::chrono::high_resolution_clock::time_point timer::print_tic;
std::chrono::high_resolution_clock::time_point timer::print_toc;


//Constructor
class_worker::class_worker(): model(), finish_line(false){
    MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);           //Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);         //Get total number of threads
    rn::rng.seed((unsigned long)world_ID);
    lnf = 1.0;
    find_initial_limits();
    resize_global_range();
    divide_global_range();
    set_initial_local_bins();
    find_current_state();
    start_counters();
    need_to_resize_global = 0;
    cout << "ID: " << world_ID << " Started OK"<<endl;

}

void class_worker::find_current_state(){
    E_idx = math::binary_search(E_bins.data(), E, E_bins.size());
    M_idx = math::binary_search(M_bins.data(), M, M_bins.size());
	in_window = check_in_window(E);
}

void class_worker::find_initial_limits(){
    //Measure, randomize and measure again to get 2 states
    double E0 = model.get_E();
    double M0 = model.get_M();
    E = E0;
    M = M0;

    while (E == E0 || M == M0){
        model.randomize_lattice();
        E = model.get_E();
        M = model.get_M();
    }
    E_min_global = fmin(E0,E);
    E_max_global = fmax(E0,E);
    M_min_global = fmin(M0,M);
    M_max_global = fmax(M0,M);
}

void class_worker::start_counters() {
    counter::MCS                = 0;
    counter::saturation         = 0;
    counter::walks              = 0;
    counter::swaps              = 0;
    counter::swap_accepts       = 0;
    counter::merges             = 0;
    timer::add_hist_volume      = 0;
    timer::check_saturation     = 0;
    timer::check_finish_line    = 0;
    timer::check_limits         = 0;
    timer::backup               = 0;
    timer::print                = 0;
	timer::swap 				= 0;
    timer::split_windows        = 0;
    timer::total_tic            = std::chrono::high_resolution_clock::now();
    timer::print_tic            = std::chrono::high_resolution_clock::now();
    flag_one_over_t             = 0;
}

void class_worker::set_initial_local_bins(){
    switch(constants::rw_dims){
        case 1:
            E_bins  = VectorXd::LinSpaced(constants::bins, E_min_local, E_max_local);
            M_bins  = VectorXd::Zero(1);
            histogram.conservativeResizeLike(MatrixXi::Zero(constants::bins,1));
            dos.conservativeResizeLike(MatrixXd::Zero(constants::bins,1));
            break;
        case 2:
            E_bins  = VectorXd::LinSpaced(constants::bins, E_min_local, E_max_local);
            M_bins  = VectorXd::LinSpaced(constants::bins, M_min_local, M_max_local);
            histogram.conservativeResizeLike(MatrixXi::Zero(constants::bins,constants::bins));
            dos.conservativeResizeLike(MatrixXd::Zero(constants::bins,constants::bins));
            break;
        default:
            cout << "Error in class_worker::set_initial_local_bins. Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            exit(1);
    }
}

void class_worker::update_global_range() {
    if (E_trial < E_min_global){
        E_min_global  = E_trial;	
		need_to_resize_global = 1;
    }
    if (E_trial > E_max_global){
        E_max_global  = E_trial;
        need_to_resize_global = 1;
    }
    if (M_trial < M_min_global){
        M_min_global = M_trial;
        need_to_resize_global = 1;
    }
    if (M_trial > M_max_global){
        M_max_global = M_trial;
        need_to_resize_global = 1;
    }
}


void class_worker::resize_global_range() {
    //Compare to the other workers to find global limits
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
}

void class_worker::divide_global_range(){
    //Update limits
    double global_range     = fabs(E_max_global - E_min_global);
    double local_range      = global_range / (world_size);
    double overlap_range    = local_range * constants::overlap_factor;

    if (world_ID == 0){
        E_min_local = E_min_global;
        E_max_local = E_min_local + local_range + overlap_range;
    }else if (world_ID == world_size - 1){
        E_max_local = E_max_global;
        E_min_local = E_max_local - local_range - overlap_range;
    }else{
        E_min_local = E_min_global + world_ID*local_range - overlap_range;
        E_max_local = E_min_global + (world_ID +1)*local_range + overlap_range;
    }
    M_min_local = M_min_global;
    M_max_local = M_max_global;
    in_window = check_in_window(E);
    if (in_window){
        E_idx = math::binary_search(E_bins.data(), E, E_bins.size());
        M_idx = math::binary_search(M_bins.data(), M, M_bins.size());
    }
}

void class_worker::resize_local_range(){
    //Scenarios
    //  1) E_bins.maxCoeff() < E_max_local, E_max_local = *E_set.end()
    //  2) E_bins.minCoeff() > E_min_local, E_min_local = *E_set.begin()
    need_to_resize_local = 0;
    if (model.discrete_model) {
        for (int w = 0; w < world_size; w++) {
            if (w == world_ID) {
                if (E_bins.maxCoeff() != E_max_local && E_bins.size() <= E_set.size()) {
                    if (debug_resize_local_range) {
                        cout << "E_max_local = " << E_max_local << " adjusted to " << fmax(E_bins.maxCoeff(),*E_set.rbegin())
                             << endl;
                        cout << "E_bins: " << E_bins.transpose() << endl;
                        cout << "E_set:  ";
                        for (auto it = E_set.begin(); it != E_set.end(); it++) {
                            if (*it == E) { cout << "["; }
                            cout << *it;
                            if (*it == E) { cout << "]"; }
                            cout << " ";
                        }
                        cout << endl;
                    }
                    E_max_local =  fmax(E_bins.maxCoeff(),*E_set.rbegin());
                    need_to_resize_local = 1;
                }
                if (E_bins.minCoeff() != E_min_local && E_bins.size() <= E_set.size()) {
                    if (debug_resize_local_range) {

                        cout << "E_min_local = " << E_min_local << " adjusted to " <<  fmin(E_bins.minCoeff(),*E_set.begin())
                             << endl;

                    }
                    E_min_local = fmin(E_bins.minCoeff(),*E_set.begin());
                    need_to_resize_local = 1;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

}

void class_worker::resize_local_bins() {
    // This function does rebinning of dos and histograms.
    // If E_set contains more than the default number of bins, then enlarge E_bins, otherwise shrink it!
    // If M_set contains more ----" " ---
    int x, y, i, j, k;
    double weight, weight_sum, dE, dM, dR, dx, dy;
    bool zero;



    //Check if we need more bins
    int E_old_size = (int) E_bins.size();
    int M_old_size = (int) M_bins.size();
    int E_new_size;
    int M_new_size;

    compute_number_of_bins(E_new_size, M_new_size);

    histogram_temp  = MatrixXi::Zero(E_new_size, M_new_size);
    dos_temp        = MatrixXd::Zero(E_new_size, M_new_size);
    dE              = (E_max_local - E_min_local) / (max(E_new_size,E_old_size) * 2.0);
    dM              = (M_max_local - M_min_local) / (max(M_new_size,M_old_size) * 2.0);
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
                    }else{zero = false;}
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
    find_current_state();
}

void class_worker::compute_number_of_bins(int & E_new_size, int & M_new_size) {
    //Scenarios:
    //  1) E_set contains more elements than the current E_bins.size() -> E_new_size = E_set.size()
    //  2) E_set contains less elements than the current E_bins.size() -> E


    if (model.discrete_model){
        //Purge elements from E_set outside of window
        for (auto it = E_set.begin(); it != E_set.end(); ) {
            if (!check_in_window(*it)){
                it = E_set.erase(it);
            }
            else {
                ++it;
            }
        }
        int spacing = (int) ceil((E_bins.tail(E_bins.size() - 1) - E_bins.head(E_bins.size()-1)).mean());
        int maxElem = (int) ceil((E_max_local - E_min_local) / spacing);
        E_new_size = max(constants::bins, (int) E_set.size()); //What if there is a jump in E_set?
        E_new_size = max(E_new_size, maxElem);
        M_new_size = max(constants::bins, (int) M_set.size());
        MPI_Allreduce(&M_new_size, &M_new_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }else{
        //Modify this for continuous models later
        E_new_size = (int) E_bins.size();
        M_new_size = (int) M_bins.size();
    }

}

bool class_worker::check_in_window(const double &x) {
    return x >= E_min_local && x <= E_max_local;
}

void class_worker::make_MC_trial()  {
    model.make_new_state(E,M, E_trial, M_trial);
    update_global_range();
}

void class_worker::acceptance_criterion(){
    //Scenarios:
    //  1) in_window = true, E_trial is in_window ->  Find idx -> MC-test
    //  2) in_window = true, E_trial not in_window->  Reject
    //  3) in_window = false,E_trial is in_window ->  Accept -> Find idx -> set in_window true
    //  4) in_window = false,E_trial not in_window->  Accept (without updating dos)

    //  4a)in_window = false, E_trial not in_window, E_trial towards window = accept, otherwise accept 50% chance?
    //  5) need_to_resize_global = true (because E_trial out of global bounds) -> accept with 50% chance?
    if (in_window){

        if (check_in_window(E_trial)){
            E_idx_trial = math::binary_search(E_bins.data(), E_trial, E_bins.size());
            M_idx_trial = math::binary_search(M_bins.data(), M_trial, M_bins.size());
            accept      = rn::uniform_double_1() < exp(dos(E_idx, M_idx) - dos(E_idx_trial, M_idx_trial));
        }else{
            accept      = false;
        }
        if (model.discrete_model){
            E_bins(E_idx) = E;
            M_bins(M_idx) = M;
//            if (E_bins(0) < E_min_local || E_bins(E_bins.size()  - 1 ) > E_max_local){
//                need_to_resize_global = 1;
//            }
            E_set.insert(E);
            M_set.insert(M);
        }
    }else{
        if (check_in_window(E_trial)){
            //Reentering the window
            accept      = true;
            in_window   = true;
            E_idx_trial = math::binary_search(E_bins.data(), E_trial, E_bins.size());
            M_idx_trial = math::binary_search(M_bins.data(), M_trial, M_bins.size());
        }else{
            //Still out of window... prefer to move towards window.
            in_window   = false;
            if (E_trial >= E && E < E_min_local){
                accept = true;
            }else if (E_trial <= E && E > E_max_local){
                accept = true;
            }else if(E_trial < E && E < E_min_local){
                accept = rn::uniform_double_1() > 0.5;
            }else if (E_trial > E && E > E_max_local){
                accept = rn::uniform_double_1() > 0.5;
            }else{
                accept = true;
            }
        }

        if (need_to_resize_global == 1){
            //Broke through global limits. Might as well explore
            if (E_trial > E && E_trial >= E_max_local) {
                accept = true;
            }else if(E_trial < E && E_trial <= E_max_local){
                accept = true;
            }
            if (M_trial > M && M_trial >= M_max_global){
                accept = true;
            }else if(M_trial < M && M_trial <= M_min_global){
                accept = true;
            }
        }
    }
}

void class_worker::accept_MC_trial() {
    E                           = E_trial;
    M                           = M_trial;
    model.flip();
    if (in_window) {
        E_idx                    = E_idx_trial;
        M_idx                    = M_idx_trial;
        dos(E_idx, M_idx)       += lnf;
        histogram(E_idx, M_idx) += 1;
    }
}

void class_worker::reject_MC_trial() {
    if (in_window) {
        dos(E_idx, M_idx)       += lnf;
        histogram(E_idx, M_idx) += 1;
    }
}

void class_worker::next_WL_iteration() {
    lnf = fmax(1e-12, lnf*constants::reduce_factor_lnf);
    histogram.fill(0);
    saturation.fill(0);
    counter::saturation = 0;
    counter::walks++;
    finish_line = lnf > constants::minimum_lnf ? 0 : 1;

}

void class_worker::rewind_to_lowest_walk(){
    int min_walks;
    MPI_Allreduce(&counter::walks, &min_walks, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    counter::walks = min_walks;
    lnf = pow(constants::reduce_factor_lnf, min_walks);

    timer::add_hist_volume  = 0;
    timer::check_finish_line= 0;
    timer::check_saturation = 0;
//    counter::MCS            = (int) ( constants::one_over_t_factor / pow(lnf, constants::one_over_t_exponent));
    counter::MCS            = (int) 1.0/lnf;
    cout << "New MCS = " << counter::MCS << endl;
    finish_line = lnf > constants::minimum_lnf ? 0 : 1;


}

void class_worker::prev_WL_iteration() {
    if (counter::walks > 0) {
        lnf /= constants::reduce_factor_lnf;
        counter::walks--;
    }

    timer::add_hist_volume  = 0;
    timer::backup           = 0;
    timer::check_finish_line= 0;
    timer::check_saturation = 0;
    if (flag_one_over_t == 0) {
        //counter::MCS = 0;
        counter::MCS            = (int) 1.0 / lnf;

        histogram.fill(0);
        saturation.fill(0);
        counter::saturation = 0;
    }else {
//        counter::MCS = (int) (constants::one_over_t_factor / pow(lnf, constants::one_over_t_exponent));
        counter::MCS            = (int) 1.0 / lnf;

    }
    finish_line = lnf > constants::minimum_lnf ? 0 : 1;

}


//Function for printing the lattice. Very easy if L is an Eigen matrix.
std::ostream &operator<<(std::ostream &os, const class_worker &worker) {
    os << setprecision(2);
    os << "ID: " << worker.world_ID
       << " Current State: " << endl
       << "     E = " << worker.E << " (" << worker.E_idx << ")"
       << "     M = " << worker.M << " (" << worker.M_idx << ")" << endl
       << "     E_trial = " << worker.E_trial << " (" << worker.E_idx_trial << ")"
       << "     M_trial = " << worker.M_trial << " (" << worker.M_idx_trial << ")" << endl
       << "     E_size  = " << worker.E_bins.size()
       << "     M_size  = " << worker.M_bins.size() << endl
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
        //os << setprecision(0) << worker.dos << endl;

    return os;
}



