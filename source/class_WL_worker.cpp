//
// Created by david on 2016-07-24.
//
#include "class_WL_worker.h"

#define profiling_sweep                	0
#define profiling_swap                 	0
#define profiling_check_global_limits  	0
#define profiling_check_convergence	    0
#define profiling_make_MC_trial 		0
#define profiling_acceptance_criterion 	0
#define debug_comp_numb_bins            1
#define debug_divide_energy             0
#define debug_resize_local_bins         0
#define debug_resize_local_range        1
#define debug_insert_state              1
#define debug_accept_trial              1
using namespace std;
using namespace Eigen;
int counter::MCS;
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


//Constructors
class_worker::class_worker(): 	model(),
                                finish_line 			(false),
                                t_sweep                (profiling_sweep)                ,
                                t_swap                 (profiling_swap)                 ,
                                t_check_global_limits  (profiling_check_global_limits)  ,
                                t_check_convergence    (profiling_check_convergence)    ,
                                t_make_MC_trial        (profiling_make_MC_trial)        ,
                                t_acceptance_criterion (profiling_acceptance_criterion)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);           //Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);         //Get total number of threads
    rn::rng.seed((unsigned long)world_ID);
    lnf = 1.0;
    find_initial_limits();
    resize_global_range();
    divide_global_range_energy();
    set_initial_local_bins();
    find_current_state();
    in_window = check_in_window(E);
    need_to_resize_global = 0;
    update_global_range();
    start_counters();
    cout << "ID: " << world_ID << " Started OK"<<endl;
}

void class_worker::find_current_state(){
    switch(constants::rw_dims){
        case 1:
            E_idx = math::binary_search(E_bins, E);
            M_idx = 0;
            break;
        case 2:
            E_idx = math::binary_search(E_bins, E);
            M_idx = math::binary_search(M_bins, M);
            break;
        default:
            cout << "Wrong rw-dimension  constants::rw_dims" << endl;
    }

}

void class_worker::find_next_state(){
    switch(constants::rw_dims) {
        case 1:
            E_idx_trial = math::binary_search(E_bins, E_trial);
            M_idx_trial = 0;
            break;
        case 2:
            E_idx_trial = math::binary_search(E_bins, E_trial);
            M_idx_trial = math::binary_search(M_bins, M_trial);
            break;
        default:
            cout << "Wrong dimension:  constants::rw_dims" << endl;
    }
    if (E_idx_trial < 0 || E_idx_trial > E_bins.size()-1){
        cout << "E_idx_trial is out of bounds!" << endl;
        cout << "E_idx_trial = " << E_idx_trial << " Size = "<< E_bins.size() << endl;
        exit(1);
    }

}

void class_worker::find_next_state(bool & dummy){
    switch (constants::rw_dims) {
        case 1:
            E_idx_trial = math::binary_search(E_bins, E_trial, E, E_idx);
            M_idx_trial = 0;
            break;
        case 2:
            E_idx_trial = math::binary_search(E_bins, E_trial, E, E_idx);
            M_idx_trial = math::binary_search(M_bins, M_trial, M, M_idx);

            break;
        default:
            cout << "Wrong dimension:  constants::rw_dims" << endl;
    }
    if (E_idx_trial < 0 || E_idx_trial > E_bins.size()-1){
        cout << "E_idx_trial is out of bounds!" << endl;
        cout << "E_idx_trial = " << E_idx_trial << " Size = "<< E_bins.size() << endl;
        exit(1);
    }

}

void class_worker::find_initial_limits(){
    //Measure, randomize and measure again to get 2 states
    double E0 = model.get_E();
    double M0 = model.get_M();

    switch (constants::rw_dims){
        case 1:
            E = E0;
            M = 0;
            while (E == E0){
                model.randomize_lattice();
                E = model.get_E();
            }
            E_min_global = fmin(E0,E);
            E_max_global = fmax(E0,E);
            M_min_global = 0;
            M_max_global = 0;
            break;
        case 2:
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
            break;
        default:
            cout << "Error: Wrong dimension  constants::rw_dims" << endl;
            MPI_Finalize();
            exit(5);
    }

}

void class_worker::start_counters() {
    counter::MCS                = 0;
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
            E_bins  = ArrayXd::LinSpaced(constants::bins, E_min_local, E_max_local);
            M_bins  = ArrayXd::Zero(1);
            histogram.conservativeResizeLike(ArrayXXi::Zero(constants::bins,1));
            dos.conservativeResizeLike(ArrayXXd::Zero(constants::bins,1));
            break;
        case 2:
            E_bins  = ArrayXd::LinSpaced(constants::bins, E_min_local, E_max_local);
            M_bins  = ArrayXd::LinSpaced(constants::bins, M_min_local, M_max_local);
            histogram.conservativeResizeLike(ArrayXXi::Zero(constants::bins,constants::bins));
            dos.conservativeResizeLike(ArrayXXd::Zero(constants::bins,constants::bins));
            break;
        default:
            cout << "Error in class_worker::set_initial_local_bins. Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            MPI_Finalize();
            exit(1);
    }
}

void class_worker::update_global_range() {
    if (E_trial < E_min_global){
        E_min_global = E_trial;
		need_to_resize_global = 1;
    }
    if (E_trial > E_max_global){
        E_max_global  = E_trial;
        need_to_resize_global = 1;
    }
    if (M_trial < M_min_global){
        M_min_global = M_trial;
        M_min_local = M_trial;
        need_to_resize_global = 1;
    }
    if (M_trial > M_max_global){
        M_max_global = M_trial;
        M_max_local = M_trial;
        need_to_resize_global = 1;
    }
}

void class_worker::resize_global_range() {
    //Compare to the other workers to find global limits
    E_min_global = fmin(E_min_local, E_min_global);
    E_max_global = fmax(E_max_local, E_max_global);
    switch (constants::rw_dims) {
        case 1:
            MPI_Allreduce(MPI_IN_PLACE, &E_min_global,  1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &E_max_global,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            M_min_global = 0;
            M_max_global = 0;
            break;
        case 2:
            MPI_Allreduce(MPI_IN_PLACE, &E_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &E_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &M_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &M_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            break;
        default:
            cout << "Error in check_windows(). Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            MPI_Finalize();
            exit(1);
    }
}

void class_worker::divide_global_range_energy(){
    //Update limits
    double global_range     = fabs(E_max_global - E_min_global);
    double local_range      = global_range / (double)(world_size);
    double x                =  constants::overlap_factor_energy;
    //Add a little bit extra if there are too many workers (Add nothing if world_size == 2, and up to local_volume if world_size == inf)
    double overlap_range = local_range * x ;// * 2.0*(world_size - 2.0 + x)/world_size;
    //Make sure the overlap covers at least a hefty percent of the global range.
    //This is needed when you have too many workers on a small global energy range.
    if (overlap_range*2 + local_range < global_range * 0.2){
        if (world_ID == 0) {
            cout << "Force readjust of local energy range" << endl;
        }
        overlap_range = (global_range * 0.2 - local_range)/2;
    }
    if      (world_ID == 0)             {overlap_range = overlap_range/1;}
    else if (world_ID == world_size - 1){overlap_range = overlap_range/1;}
    else                                {overlap_range = overlap_range/2;}

    E_min_local = E_min_global + (world_ID     *local_range)  - overlap_range;
    E_max_local = E_min_global + (world_ID +1) *local_range   + overlap_range;

    E_max_local = fmin(E_max_local,  E_max_global);
    E_min_local = fmax(E_min_local,  E_min_global);
    M_min_local = M_min_global;
    M_max_local = M_max_global;
    if(debug_divide_energy){
        if (world_ID == 0){
            cout << "Dividing according to energy range" << endl;
            cout << flush;
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for (int w = 0; w < world_size ; w++){
            if (w == world_ID){
                cout << "ID: " << world_ID << endl;
                cout << "   Global Bound= " << E_min_global << " " << E_max_global << endl;
                cout << "   Local Range = " << local_range << " overlap_range = " << overlap_range << endl;
                cout << "   Bounds E: [" <<  E_min_local << " " << E_max_local << "]" << endl;
                cout << "          M: [" <<  M_min_local << " " << M_max_local << "]" << endl;
                std::this_thread::sleep_for(std::chrono::microseconds(10));

            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

}

//void class_worker::resize_local_bins() {
//    // This function does rebinning of dos and histograms.
//    // If E_set contains more than the default number of bins, then enlarge E_bins, otherwise shrink it!
//    // If M_set contains more ----" " ---
//    int x, y, i, j, k;
//    double weight, weight_sum, dE, dM, dR, dx, dy;
//    bool zero;
//
//
//    //Check if we need more bins
//    int E_old_size = (int) E_bins.size();
//    int M_old_size = (int) M_bins.size();
//    int E_new_size;
//    int M_new_size;
//    ArrayXXd dos_temp;
//    ArrayXXi histogram_temp;
//
//    compute_number_of_bins(E_new_size, M_new_size);
//    histogram_temp  = ArrayXXi::Zero(E_new_size, M_new_size);
//    dos_temp        = ArrayXXd::Zero(E_new_size, M_new_size);
//    dE              = (E_max_local - E_min_local) / E_new_size;    //New spacing in E_bins
//    dM              = (M_max_local - M_min_local) / M_new_size;    //New spacing in
//    dR              = sqrt(dE * dE + dM * dM);
//    ArrayXd X_bins = E_bins;
//    ArrayXd Y_bins = M_bins;
//    E_bins          = ArrayXd::LinSpaced(E_new_size, E_min_local, E_max_local);
//    M_bins          = ArrayXd::LinSpaced(M_new_size, M_min_local, M_max_local);
//    //Coarsen the histogram and dos.
//    for (j = 0; j < M_new_size; j++) {
//        for (i = 0; i < E_new_size; i++) {
//            k = 0;
//            weight_sum  = 0;
//            zero        = false;
//            for (y = 0; y < M_old_size; y++) {
//                for (x = 0; x < E_old_size; x++) {
//                    if (dos(x, y) == 0) {zero = true; break;}else{zero = false;}
//                    dx = fabs(E_bins(i) - X_bins(x));             //Distance between old and new bins
//                    dy = fabs(M_bins(j) - Y_bins(y));
//                    if (dx >= dE/2) { continue; }
//                    if (dy >= dM/2) { continue; }
//                    else {
//                        weight                   = fabs(1.0 - sqrt(dx * dx + dy * dy) / dR);
//                        dos_temp(i, j)          += dos(x, y) * weight;
//                        histogram_temp(i, j)    += histogram(x, y);
//                        weight_sum              += weight;
//                        k++;
//                    }
//                }
//                if (zero) { break; }
//            }
//            if (weight_sum > 0 && !zero) {
//                dos_temp(i, j)          /= weight_sum;
//                histogram_temp(i, j)    /= k;
//            } else {
//                dos_temp(i, j)          = 0;
//                histogram_temp(i, j)    = 0;
//            }
//        }
//    }
//   if (world_ID == 0){
//        cout << dos_temp << endl << endl;
//        cout << histogram_temp << endl << endl << endl << endl;
//    }
//    dos             = dos_temp;
//    histogram       = histogram_temp;
//    math::subtract_min_nonzero(histogram);
//    if (debug_resize_local_bins)
//        for (int w = 0 ; w < world_size; w++){
//            if (w == world_ID){
//                cout << *this << endl;
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//}

void class_worker::resize_local_bins() {
    // This function does rebinning of dos and histograms.
    // If E_set contains more than the default number of bins, then enlarge E_bins, otherwise shrink it!
    // If M_set contains more ----" " ---
    int x, y, i, j;
    double dE, dM, dR, dx, dy;
    bool zero;


    //Check if we need more bins
    int E_old_size = (int) E_bins.size();
    int M_old_size = (int) M_bins.size();
    int E_new_size;
    int M_new_size;


    compute_number_of_bins(E_new_size, M_new_size);
    ArrayXXi histogram_new  = ArrayXXi::Zero(E_new_size, M_new_size);
    ArrayXXd dos_new        = ArrayXXd::Zero(E_new_size, M_new_size);
    dE              = fabs(E_max_local - E_min_local) / E_new_size;    //New spacing in E_bins
    dM              = fabs(M_max_local - M_min_local) / M_new_size;    //New spacing in M_bins
    dR              = sqrt(dE * dE + dM * dM);                         //Diagonal distance ?
    ArrayXd E_old = E_bins;
    ArrayXd M_old = M_bins;
    E_bins          = ArrayXd::LinSpaced(E_new_size, E_min_local, E_max_local);
    M_bins          = ArrayXd::LinSpaced(M_new_size, M_min_local, M_max_local);
    ArrayXXd weight = ArrayXXd::Zero(E_new_size, M_new_size);
    ArrayXXi count  = ArrayXXi::Zero(E_new_size, M_new_size);
    //Coarsen the histogram and dos.
    for (y = 0; y < M_old_size; y++) {
        for (x = 0; x < E_old_size; x++) {
            for (j = 0; j < M_new_size; j++) {
                for (i = 0; i < E_new_size; i++) {
                    dx = fabs(E_bins(i) - E_old(x));   //Distance between old and new bins
                    dy = fabs(M_bins(j) - M_old(y));   //Distance between old and new bins
                    //Distance between old and new bins should not exceed dE/2
                    //This is so to avoid double counting
                    if (dx >= dE) { continue; }
                    if (dy >= dM) { continue; }
                    double w                  = fabs(1.0 - sqrt(dx * dx + dy * dy) / dR);
//                    weight(i,j)              += fabs(1.0 - 4*dx*dy/dE/dM);
                    weight(i,j)              += w;
                    count (i,j)              += 1;
                    dos_new(i,j)             += w * dos(x,y);
                    histogram_new(i,j)       += histogram(x,y);
                }
            }
        }
    }
    //We have now inserted all old entries to the new dos and histogram, and we only need to divide by the weight.
    for (j = 0; j < M_new_size; j++) {
        for (i = 0; i < E_new_size; i++) {
            dos_new(i,j)        = weight(i,j) > 0  ? dos_new(i,j)/weight(i,j) : 0;
            histogram_new(i,j)  = count(i,j) > 0  ? histogram_new(i,j)  / count(i,j)  : 0;

        }
    }
    dos             = dos_new;
    histogram       = histogram_new;
    if (debug_resize_local_bins)
    for (int w = 0 ; w < world_size; w++){
        if (w == world_ID){
            cout << *this << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

}

void class_worker::compute_number_of_bins(int & E_new_size, int & M_new_size) {
    //Scenarios:
    //  1) E_set contains more elements than the current E_bins.size() -> E_new_size = E_set.size()
    //  2) E_set contains less elements than the current E_bins.size() -> E

    if (model.discrete_model){
        int i = 0;
        ArrayXd E_set_2_vector(E_set.size());
        ArrayXd M_set_2_vector(M_set.size());
        for (auto it = E_set.begin(); it != E_set.end(); it++) {
            E_set_2_vector(i) = *it;
            i++;
        }
        i = 0;
        for (auto it = M_set.begin(); it != M_set.end(); it++) {
            M_set_2_vector(i) = *it;
            i++;
        }

        MPI_Allreduce(MPI_IN_PLACE, &M_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &M_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        M_min_local = M_min_global;
        M_max_local = M_max_global;
        double E_spacing =  math::typical_spacing(E_set_2_vector);
        double M_spacing =  math::typical_spacing(M_set_2_vector);
//        double E_spacing     =  fmax(E_spacing_set,math::typical_spacing(E_bins));
//        double M_spacing     =  fmax(M_spacing_set,math::typical_spacing(M_bins));
        int E_maxElem     = (int) ceil((E_max_local - E_min_local + E_spacing) / E_spacing);
        int M_maxElem     = (int) ceil((M_max_global - M_min_global + M_spacing) / M_spacing);

        //Purge elements from E_set outside of window
        for (auto it = E_set.begin(); it != E_set.end(); ) {
            if (!check_in_window(*it)){
                it = E_set.erase(it);
            }
            else {
                ++it;
            }
        }
        E_new_size = max(constants::bins , E_maxElem); //What if there is a jump in E_set?
//        E_new_size = max(E_new_size      , E_maxElem);
        M_new_size = max(constants::bins , M_maxElem);
//        M_new_size = max(M_new_size      , M_maxElem);
        M_new_size = constants::rw_dims == 1 ? 1 : M_new_size;
        MPI_Allreduce( MPI_IN_PLACE, &M_new_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (debug_comp_numb_bins) {
            for (i = 0; i < world_size; i++) {
                if (i == world_ID) {
                    cout << "ID: " << world_ID << endl
                         << "    Spacing         = " << E_spacing        << " " << M_spacing     << endl
                         << "    MaxElems        = " << E_maxElem        << " " << M_maxElem     << endl
                         << "    Set sizes       = " << E_set.size()     << " " << M_set.size()  << endl
                         << "    E local bounds  = " << E_min_local      << " " << E_max_local   << endl
                         << "    E global bounds = " << E_min_global     << " " << E_max_global  << endl
                         << "    E_bins          = " << E_bins.transpose() << endl
                         << "    M_bins          = " << M_bins.transpose() << endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }

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
	t_make_MC_trial.tic();
    model.make_new_state(E,M, E_trial, M_trial);
    update_global_range();
	t_make_MC_trial.toc();
}

void class_worker::insert_state(){
    if (model.discrete_model && counter::merges == 0) {
        if (debug_insert_state && in_window) {
            int E_idx_temp = E_idx;
            find_current_state();
            if (E_idx_temp != E_idx) {
                cout << "ID: " << world_ID <<  " Insert state failed: State mismatch!" << endl;
                cout << "E_idx reported: " << E_idx_temp << endl;
                cout << "E_idx reality : " << E_idx << endl;
                cout << "Current E: " << E << "  E_bins: " << E_bins.transpose() << endl;
                cout.flush();
                std::this_thread::sleep_for(std::chrono::microseconds(1000));
                exit(1);
            }
        }

        E_bins(E_idx) = E;
        M_bins(M_idx) = M;
        E_set.insert(E);
        M_set.insert(M);
    }
}

void class_worker::walk_away_from_window(){
//    if ((E_trial >= E && E_trial >= E_max_local) || (constants::rw_dims == 2  && M_trial >= M && M_trial >= M_max_local)) {
    if ((E_trial >= E && E_trial >= E_max_local) ) {
        accept = true;
//    } else if ((E_trial <= E && E_trial <= E_min_local) || (constants::rw_dims == 2  && M_trial <= M && M_trial <= M_min_local)) {
    } else if ((E_trial <= E && E_trial <= E_min_local) ) {
        accept = true;
    }else {
        //Either you're taking a step towards your window, in which case you can't accept,
        //...or you're already inside your window, in which case you need to leave. Go to the nearest boundary!
        if (check_in_window(E)){
            bool move_up = fabs(E - E_min_local ) >= fabs(E - E_max_local);
            if (move_up){
                accept = E_trial >= E ;
            }else{
                accept = E_trial <= E ;
            }
        }else{
            accept = false;
        }
    }
}

void class_worker::walk_towards_window(){
    //If E_trial not in window, and neither is E. Then we'd rather find our window.
    if        (E_trial >= E && E < E_min_local) {
        //Step towards window from below
        accept = true;
    } else if (E_trial <= E && E > E_max_local) {
        //Step towards window from above
        accept = true;
    } else{
        accept = rn::uniform_double_1() < 0.1;

//        accept = false;

    }
}

//void class_worker::acceptance_criterion2() {
//    //Scenarios:
//    //  1) in_window = true, E_trial is in_window ->  Find idx -> MC-test
//    //  2) in_window = true, E_trial not in_window->  Reject
//    //  3) in_window = false,E_trial is in_window ->  Accept -> Find idx -> set in_window true
//    //  4) in_window = false,E_trial not in_window->  Accept (without updating dos)
//
//    //  4a)in_window = false, E_trial not in_window, E_trial towards window = accept, otherwise accept 50% chance?
//    //  5) need_to_resize_global = true (because E_trial out of global bounds) -> accept with 50% chance?
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (!need_to_resize_global) {
//        if (check_in_window(E_trial)) {
//            if (in_window) {
//                if (debug_accept_trial) {
//                    int E_idx_temp = E_idx;
//                    find_current_state();
//                    if (E_idx_temp != E_idx) {
//                        cout << "ID: " << world_ID <<  " Accept trial failed: State mismatch!" << endl;
//                        cout << "E_idx reported: " << E_idx_temp << endl;
//                        cout << "E_idx reality : " << E_idx << endl;
//                        cout << "Current E: " << E << "  E_bins: " << E_bins.transpose() << endl;
//                        cout.flush();
//                        std::this_thread::sleep_for(std::chrono::microseconds(1000));
////                        exit(1);
//                    }
//                }
////                insert_state();
//                find_next_state(in_window);
//                accept = rn::uniform_double_1() < exp(dos(E_idx, M_idx) - dos(E_idx_trial, M_idx_trial));
////                cout << fixed << showpoint << setprecision(1) << *this << endl;
//
//
//                if (!std::is_sorted(E_bins.data(),E_bins.data()+(int)E_bins.size())){
//                    cout << "Unsorted!" << endl;
//                    cout << E_bins.transpose() << endl;
//                    exit(0);
//                }
//
//            } else { //If  we are currently not in window, but E_trial is! then reenter
//                in_window = true;
//                accept = true;
//                if (debug_accept_trial) {
//                    int E_idx_temp = E_idx;
//                    find_current_state();
//                    if (E_idx_temp != E_idx) {
//                        cout << "ID: " << world_ID <<  " Reenter failed: State mismatch!" << endl;
//                        cout << "E_idx reported: " << E_idx_temp << endl;
//                        cout << "E_idx reality : " << E_idx << endl;
//                        cout << "Current E: " << E << "  E_bins: " << E_bins.transpose() << endl;
//                        cout.flush();
//                        std::this_thread::sleep_for(std::chrono::microseconds(1000));
////                        exit(1);
//                    }
//                }
//                find_next_state();
//            }
//        } else {
//            if (in_window){ //If E_trial not in window, but E is already, then reject
////                insert_state();
//                accept == false;
//            }else{          //If E_trial not in window, Then we'd rather find our window.
//                walk_towards_window();
//            }
//
//        }
//
//    } else { //If need to resize we might as well explore the E-M landscape, so use different strategy
//        //Windows should no longer be relevant, everybody should just try to go outwards
//        //Accept if taking a step outward
//        in_window = false;
////        insert_state();
//        walk_away_from_window();
//
//    }
//}


void class_worker::acceptance_criterion(){
    //Scenarios:
    //  1) in_window = true, E_trial is in_window ->  Find idx -> MC-test
    //  2) in_window = true, E_trial not in_window->  Reject
    //  3) in_window = false,E_trial is in_window ->  Accept -> Find idx -> set in_window true
    //  4) in_window = false,E_trial not in_window->  Accept (without updating dos)

    //  4a)in_window = false, E_trial not in_window, E_trial towards window = accept, otherwise accept 50% chance?
    //  5) need_to_resize_global = true (because E_trial out of global bounds) -> accept with 50% chance?
    if (!need_to_resize_global) {
        if (in_window) {
            t_acceptance_criterion.tic();
            insert_state();
            if (check_in_window(E_trial)) {
                find_next_state(in_window);
                P_accept  =  exp(dos(E_idx, M_idx) - dos(E_idx_trial, M_idx_trial));
                accept = rn::uniform_double_1() < P_accept;
            } else {
                accept = false;
            }
            t_acceptance_criterion.toc();
        } else {
            if (check_in_window(E_trial)) {
                //Reentering the window
                accept = true;
                in_window = true;
                find_next_state();

            } else {
                //Still out of window... prefer to move towards window.
                in_window = false;
                walk_towards_window();
            }
        }
    }else{
        //Broke through global limits. Might as well explore
        //The current state will probably be out of window, so E_idx points to either the first
        //or last element in E_bins if we try to find it. Even so, let's insert it anyway.
        in_window = false;
        find_current_state();
        insert_state();
        walk_away_from_window();
    }

}

void class_worker::accept_MC_trial() {
    E                           = E_trial;
    M                           = M_trial;
    model.flip();
    if (in_window) {
        E_idx                       = E_idx_trial;
        M_idx                       = M_idx_trial;
        dos(E_idx, M_idx)       += lnf*P_accept ;
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
    saturation.clear();
    counter::walks++;
    finish_line = lnf > constants::minimum_lnf ? 0 : 1;

}

void class_worker::rewind_to_lowest_walk(){
    int min_walks;
    MPI_Allreduce(&counter::walks, &min_walks, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    counter::walks = min_walks;
    if (flag_one_over_t){

    }
    lnf = pow(constants::reduce_factor_lnf, min_walks);

    timer::add_hist_volume  = 0;
    timer::check_finish_line= 0;
    timer::check_saturation = 0;
    saturation.clear();
    counter::MCS            = (int) (1.0/lnf);
    cout << "New MCS = " << counter::MCS << endl;
    flag_one_over_t = 0;
    finish_line = lnf > constants::minimum_lnf ? 0 : 1;
}


void class_worker::rewind_to_zero(){
    counter::walks = 0;
    lnf = 1;
    start_counters();
    counter::merges = 1;
    finish_line = 0;
    dos.fill(0);
    timer::check_saturation = 0;
    histogram = ArrayXXi::Zero(dos.rows(), dos.cols());
    saturation.clear();
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
        saturation.clear();
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
       << " Current State: "
       << " Need to resize: " << worker.need_to_resize_global << endl
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
       << "     E_bins " << "[" << worker.E_bins.size() << "] "<< worker.E_bins.transpose() << endl
       << "     M_bins " << "[" << worker.M_bins.size() << "] "<< worker.M_bins.transpose() << endl
       << "     E_set  " << "[" << worker.E_set.size() << "] ";
        for (auto it = worker.E_set.begin(); it != worker.E_set.end(); it++) {
            if (*it == worker.E){ os << "[";}
            os << *it ;
            if (*it == worker.E){ os << "]";}
            os << " ";

        }
        os << endl
        << "     M_set : " << "[" << worker.M_set.size() << "] ";
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



