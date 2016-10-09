//
// Created by david on 2016-07-24.
//
#include "class_WL_worker.h"

#define profiling_total                	1
#define profiling_print                	1
#define profiling_sweep                	1
#define profiling_swap                 	1
#define profiling_merge                	1
#define profiling_help_setup           	1
#define profiling_help                 	1
#define profiling_divide_range        	1
#define profiling_check_convergence	    1
#define profiling_make_MC_trial 		0
#define profiling_acceptance_criterion 	0

#define debug_comp_numb_bins            0
#define debug_divide_energy             0
#define debug_resize_local_bins         0
#define debug_resize_local_range        0
#define debug_insert_state              0
#define debug_accept_trial              0

using namespace std;
using namespace Eigen;
int counter::MCS;
int counter::walks;
int counter::swaps;
int counter::swap_accepts;
int counter::vol_merges;
int timer::increment;
int timer::add_dos;
int timer::add_hist_volume;
int timer::check_saturation;
int timer::check_finish_line;
int timer::backup;
int timer::print;
int timer::swap;
int timer::check_help;
int timer::take_help;
int timer::setup_help;
int timer::divide_range;


std::ostream& operator<< (std::ostream& out, const std::vector<state>& v) {
    if (!v.empty()) {
        out << "[ ";
        for (int i = 0; i < v.size(); i++) {
            out << "(" << v[i].E_idx << "," << v[i].M_idx << ") ";
        }
        out << "]";
    }
    return out;
}

//Constructors
class_worker::class_worker(int & id, int & size):
                                world_ID(id),
                                world_size(size),
                                model(constants::L),
                                finish_line 			(false),
                                t_total                (profiling_total,                3,"Total Time"),
                                t_print                (profiling_print,                3,"Time"),
                                t_sweep                (profiling_sweep,                3,"t_sweep"),
                                t_swap                 (profiling_swap,                 3,"t_swap" ),
                                t_merge                (profiling_merge,                3,"t_merge" ),
                                t_help_setup           (profiling_help_setup,           3,"t_help_s" ),
                                t_help                 (profiling_help,                 3,"t_help" ),
                                t_divide_range         (profiling_divide_range,         3,"t_divr" ),
                                t_check_convergence    (profiling_check_convergence,    3,"t_conv") ,
                                t_make_MC_trial        (profiling_make_MC_trial,        3,"t_mkMC") ,
                                t_acceptance_criterion (profiling_acceptance_criterion, 3,"t_accr")
{
    rn::rng.seed((unsigned long)world_ID);
    lnf = 1.0;
    find_initial_limits();
    resize_global_range();
    divide_global_range_uniformly();
    set_initial_local_bins();
    find_current_state();
    need_to_resize_global = 0;
    update_global_range();
    start_counters();
    rate_increment = 1;
    cout << "ID: " << world_ID << " Started OK"<<endl;
}

void class_worker::find_current_state(){
    state_in_window = check_in_window(E);
    if(state_in_window) {
        switch (constants::rw_dims) {
            case 1:
                E_idx = math::binary_search_nearest(E_bins, E);
                M_idx = 0;
                break;
            case 2:
                E_idx = math::binary_search_nearest(E_bins, E);
                M_idx = math::binary_search_nearest(M_bins, M);
                break;
            default:
                cout << "Wrong rw-dimension  constants::rw_dims" << endl;
        }
    }
}


void class_worker::find_next_state_exact(){
    switch(constants::rw_dims) {
        case 1:
            E_idx_trial = math::binary_search_exact(E_bins, E_trial);
            M_idx_trial = 0;
            state_is_valid = E_idx_trial != -1;
            break;
        case 2:
            E_idx_trial = math::binary_search_exact(E_bins, E_trial);
            M_idx_trial = math::binary_search_exact(M_bins, M_trial);
            state_is_valid = !(E_idx_trial == -1 || M_idx_trial == -1);
            break;
        default:
            cout << "Wrong dimension:  constants::rw_dims" << endl;
    }

}

void class_worker::find_next_state(){

    switch(constants::rw_dims) {
        case 1:
            E_idx_trial = math::binary_search_nearest(E_bins, E_trial);
            M_idx_trial = 0;
            break;
        case 2:
            E_idx_trial = math::binary_search_nearest(E_bins, E_trial);
            M_idx_trial = math::binary_search_nearest(M_bins, M_trial);
            break;
        default:
            cout << "Wrong dimension:  constants::rw_dims" << endl;
    }

}

void class_worker::find_next_state(bool dummy){
    switch (constants::rw_dims) {
        case 1:
            E_idx_trial = math::binary_search_nearest(E_bins, E_trial, E, E_idx);
            M_idx_trial = 0;
            break;
        case 2:
            E_idx_trial = math::binary_search_nearest(E_bins, E_trial, E, E_idx);
            M_idx_trial = math::binary_search_nearest(M_bins, M_trial, M, M_idx);
            break;
        default:
            cout << "Wrong dimension:  constants::rw_dims" << endl;
    }
}

void class_worker::find_initial_limits(){
    //Measure, randomize and measure again to get 2 states
    double E0 = model.get_E();
    double M0 = model.get_M();
    E_set.insert(E0);
    M_set.insert(M0);


    switch (constants::rw_dims){
        case 1:
            E = E0;
            M = 0;
            while (E == E0){
                model.randomize_lattice();
                E = model.get_E();
                E_set.insert(E);
            }
            E_min_global = *E_set.begin() ;// fmin(E0,E);
            E_max_global = *E_set.rbegin();//E0,E);
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
                E_set.insert(E);
                M_set.insert(M);

            }
            E_min_global = *E_set.begin() ;
            E_max_global = *E_set.rbegin();
            M_min_global = *M_set.begin() ;
            M_max_global = *M_set.rbegin();
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
    counter::vol_merges         = 0;
    timer::increment            = 0;
    timer::add_dos              = 0;
    timer::add_hist_volume      = 0;
    timer::check_saturation     = 0;
    timer::check_finish_line    = 0;
    timer::backup               = 0;
    timer::print                = 0;
    timer::swap 				= 0;
    timer::check_help        	= 0;
    timer::take_help        	= 0;
    timer::setup_help			= 0;
    timer::divide_range         = 0;
    flag_one_over_t             = 0;
}

void class_worker::rewind_timers(){
    timer::increment            = 0;
    timer::add_dos              = 0;// math::mod(counter::MCS, constants::rate_add_hist_volume  );
    timer::add_hist_volume      = 0;// math::mod(counter::MCS, constants::rate_add_hist_volume  );
    timer::check_saturation     = 0;// math::mod(counter::MCS, constants::rate_check_saturation );
    timer::check_finish_line    = 0;// math::mod(counter::MCS, constants::rate_check_finish_line);
    timer::backup               = 0;// math::mod(counter::MCS, constants::rate_backup_data      );
    timer::print                = 0;// math::mod(counter::MCS, constants::rate_print_status     );
    timer::swap 				= 0;// math::mod(counter::MCS, constants::rate_swap             );
    timer::check_help        	= 0;// math::mod(counter::MCS, constants::rate_take_help        );
    timer::take_help        	= 0;// math::mod(counter::MCS, constants::rate_take_help        );
    timer::setup_help			= 0;// math::mod(counter::MCS, constants::rate_setup_help       );
    timer::divide_range         = 0;// math::mod(counter::MCS, constants::rate_divide_range     );
}

void class_worker::set_initial_local_bins(){
    random_walk.clear();
    switch(constants::rw_dims){
        case 1:
            E_bins << E_set;
//            copy_set(E_bins ,E_set);
            M_bins  = ArrayXd::Zero(1);
            histogram.conservativeResizeLike(ArrayXXi::Zero(E_bins.size(),M_bins.size()));
            dos.conservativeResizeLike(ArrayXXd::Zero(E_bins.size(),M_bins.size()));
//            histogram_incr = histogram;
            break;
        case 2:
            E_bins << E_set;
            M_bins << M_set;
//            copy_set(E_bins, E_set);
//            copy_set(M_bins, M_set);
            histogram.conservativeResizeLike(ArrayXXi::Zero(E_bins.size(),M_bins.size()));
            dos.conservativeResizeLike(ArrayXXd::Zero(E_bins.size(),M_bins.size()));
//            histogram_incr = histogram;
//            E_bins  = ArrayXd::LinSpaced(constants::bins, E_min_local, E_max_local);
//            M_bins  = ArrayXd::LinSpaced(constants::bins, M_min_local, M_max_local);
//            histogram.conservativeResizeLike(ArrayXXi::Zero(constants::bins,constants::bins));
//            dos.conservativeResizeLike(ArrayXXd::Zero(constants::bins,constants::bins));
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

void class_worker::divide_global_range_uniformly(){
    //Update limits
    double global_range     = fabs(E_max_global - E_min_global);
    double local_range      = global_range / (double)(world_size);
    double x                = constants::overlap_factor_energy;
    //Add a little bit extra if there are too many workers (Add nothing if world_size == 2, and up to local_volume if world_size == inf)
    double overlap_range = local_range * x ;// * 2.0*(world_size - 2.0 + x)/world_size;
    //Make sure the overlap covers at least a hefty percent of the global range.
    //This is needed when you have too many workers on a small global energy range.
    if (overlap_range*2 + local_range < global_range * 0.2){
        if (world_ID == 0) {
            cout << "Forced readjustment of local energy range" << endl;
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


void class_worker::synchronize_sets(){
    //This function collects all known energies from all sets
    vector<double> E_vector(E_set.begin(),E_set.end());
    vector<double> M_vector(M_set.begin(),M_set.end());
    vector<int> E_sizes((unsigned long)world_size);
    vector<int> M_sizes((unsigned long)world_size);
    int Esize = (int) E_vector.size();
    int Msize = (int) M_vector.size();
    MPI_Allgather(&Esize,1,MPI_INT, E_sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&Msize,1,MPI_INT, M_sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
    vector<double> E_recv;
    vector<double> M_recv;
    for (int w = 0; w < world_size; w++){
        if (w == world_ID){
            E_recv = E_vector;
            M_recv = M_vector;
        }else{
            E_recv.resize((unsigned long) E_sizes[w]);
            M_recv.resize((unsigned long) M_sizes[w]);
        }
        MPI_Bcast(E_recv.data(),(int)E_sizes[w], MPI_DOUBLE,w,MPI_COMM_WORLD);
        MPI_Bcast(M_recv.data(),(int)M_sizes[w], MPI_DOUBLE,w,MPI_COMM_WORLD);
        E_set.insert(E_recv.begin(),E_recv.end());
        M_set.insert(M_recv.begin(),M_recv.end());
    }

}

void class_worker::adjust_local_bins() {
    // This function does rebinning of dos and histograms.

    //Snap the local boundaries to existing energies
    ArrayXd E_set_to_array;
    ArrayXd M_set_to_array;
    E_set_to_array << E_set;
    M_set_to_array << M_set;

    int E_idx_min = math::binary_search_nearest(E_set_to_array, E_min_local);
    int E_idx_max = math::binary_search_nearest(E_set_to_array, E_max_local);
    E_min_local   = E_set_to_array(E_idx_min);
    E_max_local   = E_set_to_array(E_idx_max);
    int M_idx_min = math::binary_search_nearest(M_set_to_array, M_min_local);
    int M_idx_max = math::binary_search_nearest(M_set_to_array, M_max_local);
    M_min_local   = M_set_to_array(M_idx_min);
    M_max_local   = M_set_to_array(M_idx_max);


    ArrayXd E_old = E_bins;
    ArrayXd M_old = M_bins;

    E_bins                  = E_set_to_array.segment(E_idx_min, E_idx_max-E_idx_min+1);
    M_bins                  = M_set_to_array.segment(M_idx_min, M_idx_max-M_idx_min+1);
    histogram               = ArrayXXi::Zero(E_bins.size(), M_bins.size());
    dos                     = ArrayXXd::Zero(E_bins.size(), M_bins.size());


//    cout << "Hej 4" << endl;
//    int x, y, i, j;
//    double dR, dx, dy;
//    int E_old_size = (int) E_bins.size();
//    int M_old_size = (int) M_bins.size();
//    ArrayXXi histogram_new  = ArrayXXi::Zero(E_bins.size(), M_bins.size());
//    ArrayXXd dos_new        = ArrayXXd::Zero(E_bins.size(), M_bins.size());
//
//    ArrayXXd weight = ArrayXXd::Zero(E_bins.size(), M_bins.size());
//    ArrayXXi count  = ArrayXXi::Zero(E_bins.size(), M_bins.size());

//    double dE_dn,dE_up,dM_dn,dM_up;
//    for (y = 0; y < M_old_size; y++) {
//        for (x = 0; x < E_old_size; x++) {
//            if (dos(x,y) == 0){continue;}
//            dE_dn = x > 0               ? E_old(x-1) - E_old(x) : 0;
//            dE_up = x < E_old_size-1    ? E_old(x+1) - E_old(x) : 0;
//            dM_dn = y > 0               ? M_old(y-1) - M_old(y) : 0;
//            dM_up = y < M_old_size-1    ? M_old(y+1) - M_old(y) : 0;
//            dR = sqrt(fmax(dE_dn*dE_dn, dE_up*dE_up) + fmax(dM_dn*dM_dn, dM_up*dM_up));
//            for (j = 0; j < M_new_size; j++) {
//                for (i = 0; i < E_new_size; i++) {
//                    dx = E_bins(i) - E_old(x);   //Distance between old and new bins
//                    dy = M_bins(j) - M_old(y);   //Distance between old and new bins
//                    //Distance between old and new bins should not exceed dE
//                    //This is so to avoid double counting
//                    if (dx <= dE_dn || dx >= dE_up) { continue; }
//                    if (dy <= dM_dn || dy >= dM_up) { continue; }
//                    double w                  = fabs(1.0 - sqrt(dx * dx + dy * dy) / dR);
////                    weight(i,j)              += fabs(1.0 - 4*dx*dy/dE/dM);
//                    weight(i,j)              += w;
//                    count (i,j)              += 1;
//                    dos_new(i,j)             += w * dos(x,y);
//                    histogram_new(i,j)       += histogram(x,y);
//                }
//            }
//        }
//    }
//
//    //We have now inserted all old entries to the new dos and histogram, and we only need to divide by the weight.
//    for (j = 0; j < M_new_size; j++) {
//        for (i = 0; i < E_new_size; i++) {
//            dos_new(i,j)        = weight(i,j) > 0  ? dos_new(i,j)/weight(i,j) : 0;
//            histogram_new(i,j)  = count(i,j)  > 0  ? histogram_new(i,j)  / count(i,j)  : 0;
//
//        }
//    }
//    dos             = dos_new;
//    histogram       = histogram_new;
    if (debug_resize_local_bins) {
        for (int w = 0; w < world_size; w++) {
            if (w == world_ID) {
                cout << *this << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

//void class_worker::compute_number_of_bins(int & E_new_size, int & M_new_size) {
//    //Scenarios:
//    //  1) E_set contains more elements than the current E_bins.size() -> E_new_size = E_set.size()
//    //  2) E_set contains less elements than the current E_bins.size() -> E
//
//    if (model.discrete_model){
//        int i = 0;
//        ArrayXd E_set_2_vector(E_set.size());
//        ArrayXd M_set_2_vector(M_set.size());
//        for (auto it = E_set.begin(); it != E_set.end(); it++) {
//            E_set_2_vector(i) = *it;
//            i++;
//        }
//        i = 0;
//        for (auto it = M_set.begin(); it != M_set.end(); it++) {
//            M_set_2_vector(i) = *it;
//            i++;
//        }
//
//        MPI_Allreduce(MPI_IN_PLACE, &M_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//        MPI_Allreduce(MPI_IN_PLACE, &M_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//        M_min_local = M_min_global;
//        M_max_local = M_max_global;
//        double E_spacing =  math::typical_spacing(E_set_2_vector);
//        double M_spacing =  math::typical_spacing(M_set_2_vector);
//
//
//        MPI_Allreduce(MPI_IN_PLACE, &E_spacing, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        E_spacing = E_spacing/(world_size);
//        int E_global_size = 1 + (int)((E_max_global - E_min_global)/nearbyint(E_spacing));
//        //Now try estimating the global size (i.e. the size of the global spectrum)
//        //Adjust local boundaries
//        ArrayXd E_bins_global = ArrayXd::LinSpaced(E_global_size, E_min_global, E_max_global);
//        int E_idx_min = math::binary_search_nearest(E_bins_global, E_min_local);
//        int E_idx_max = math::binary_search_nearest(E_bins_global, E_max_local);
//        E_min_local   = E_bins_global(E_idx_min);
//        E_max_local   = E_bins_global(E_idx_max);
//
//        int E_maxElem     = (int) ceil((E_max_local - E_min_local + E_spacing) / E_spacing);
//        int M_maxElem     = (int) ceil((M_max_global - M_min_global + M_spacing) / M_spacing);
//
//        //Purge elements from E_set outside of window
//        for (auto it = E_set.begin(); it != E_set.end(); ) {
//            if (!check_in_window(*it)){
//                it = E_set.erase(it);
//            }
//            else {
//                ++it;
//            }
//        }
//        E_new_size = max(constants::bins , E_maxElem); //What if there is a jump in E_set?
//        M_new_size = max(constants::bins , M_maxElem);
//        M_new_size = constants::rw_dims == 1 ? 1 : M_new_size;
//        MPI_Allreduce( MPI_IN_PLACE, &M_new_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//        if (debug_comp_numb_bins) {
//            for (i = 0; i < world_size; i++) {
//                if (i == world_ID) {
//                    cout << "ID: " << world_ID << endl
//                         << "    Spacing         = " << E_spacing        << " " << M_spacing     << endl
//                         << "    MaxElems        = " << E_maxElem        << " " << M_maxElem     << endl
//                         << "    Set sizes       = " << E_set.size()     << " " << M_set.size()  << endl
//                         << "    E local bounds  = " << E_min_local      << " " << E_max_local   << endl
//                         << "    E global bounds = " << E_min_global     << " " << E_max_global  << endl
//                         << "    E_bins          = " << E_bins.transpose() << endl
//                         << "    M_bins          = " << M_bins.transpose() << endl;
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//            }
//        }
//
//
//
//
//    }else{
//        //Modify this for continuous models later
//        E_new_size = (int) E_bins.size();
//        M_new_size = (int) M_bins.size();
//    }
//}

void class_worker::make_MC_trial()  {
	t_make_MC_trial.tic();
    model.make_new_state(E,M, E_trial, M_trial);
    update_global_range();
	t_make_MC_trial.toc();
}

void class_worker::insert_state(double new_E, double new_M){
    if (model.discrete_model) {
        if (debug_insert_state && state_in_window) {
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
        E_set.insert(new_E);
        M_set.insert(new_M);
    }
}

void class_worker::walk_away_from_window(){
    if ((E_trial >= E && E_trial >= E_max_local) ) {
        accept = true;
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
        accept = rn::uniform_double_1() < 0.01;
    }
}

void class_worker::acceptance_criterion(){
    t_acceptance_criterion.tic();
    if (!need_to_resize_global) {
        if (state_in_window) {
            if (check_in_window(E_trial)) {
                find_next_state_exact();
                if(state_is_valid){
                    //The energy exists, proceed with normal WL
                    accept = rn::uniform_double_1() < fmin(1, exp(dos(E_idx, M_idx) - dos(E_idx_trial, M_idx_trial)));

                }
                else {
                    //A new energy has been discovered in window!
                    insert_state(E_trial, M_trial);
                    need_to_resize_global = 1;
                    accept          = true;
                  }
            } else {
                //The proposed E_trial is out of window, don't accept!
                accept = false;
            }
        } else {
            if (check_in_window(E_trial)) {
                //Has been outside, now reentering the window
                accept          = true;
                state_in_window = true;
//                cout << "ID: " << world_ID << " Reentering ("<< counter::MCS << "): " << random_walk <<endl;
                find_next_state_exact();
                if (!state_is_valid){
                    //New energy has been found, insert
                    insert_state(E_trial, M_trial);
                    need_to_resize_global = 1;
                }
            }else{
                //Still out of window... prefer to move towards window.
//                cout << "ID: " << world_ID << " Walking towards window("<< counter::MCS << "): " << random_walk << endl;

                state_in_window = false;
                walk_towards_window();
            }
        }
    }else{
        //Broke through global limits. Might as well explore
        //The current state will probably be out of window, so E_idx points to either the first
        //or last element in E_bins if we try to find it. Even so, let's insert it anyway.
        state_in_window = false;
        state_is_valid  = false;
        insert_state(E_trial, M_trial);
        walk_away_from_window();
    }
    t_acceptance_criterion.toc();


}

void class_worker::accept_MC_trial() {
    E                           = E_trial;
    M                           = M_trial;
    model.flip();
    if (state_in_window && state_is_valid) {
        E_idx                       = E_idx_trial;
        M_idx                       = M_idx_trial;
        if (++timer::increment >= rate_increment) {
            timer::increment = 0;
            random_walk.push_back({E_idx,M_idx});


        }
    }
}

void class_worker::reject_MC_trial() {
    if (state_in_window && state_is_valid) {
        if (++timer::increment >= rate_increment) {
            timer::increment = 0;
            random_walk.push_back({E_idx,M_idx});
        }
    }
}

void class_worker::set_rate_increment(){
    double dos_width = 0;
    double dos_height = 0;

    for (int i = 0; i <dos.rows(); i++){
        if((dos.row(i) != 0).any()){
            dos_height += 1;
        }
    }
    for (int i = 0; i <dos.cols(); i++){
        if((dos.col(i) != 0).any()){
            dos_width += 1;
        }
    }
//    rate_increment = max(1, (int)std::sqrt(fmax(dos_width,dos_height)));
    rate_increment = constants::N;
}

void class_worker::next_WL_iteration() {
    lnf = fmax(1e-12, lnf*constants::reduce_factor_lnf);
    histogram.setZero();
    saturation.clear();
    random_walk.clear();
    counter::walks++;
}

void class_worker::rewind_to_lowest_walk(){
    int min_walks;
    MPI_Allreduce(&counter::walks, &min_walks, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    counter::walks = min_walks;
    lnf = pow(constants::reduce_factor_lnf, min_walks);
    finish_line = 0;
    saturation.clear();
    counter::MCS            = (int) (1.0/lnf);
//    rewind_timers();
    random_walk.clear();
    help.reset();
}

void class_worker::rewind_to_zero(){
    counter::walks = 0;
    lnf = 1;
    int save_vol_merges  = counter::vol_merges;
    start_counters();
    counter::vol_merges  = save_vol_merges;
    finish_line = 0;
    dos.setZero();
    histogram.setZero();
    random_walk.clear();
    saturation.clear();
    rewind_timers();
    help.reset();

    flag_one_over_t = 0;
}

void class_worker::prev_WL_iteration() {
    if (counter::walks > 0) {
        lnf /= constants::reduce_factor_lnf;
        counter::walks--;
    }
    counter::MCS            = counter::walks == 0 ? 0 : (int) (1.0 / lnf);
    flag_one_over_t         = 0;
    histogram.setZero();
    saturation.clear();
    random_walk.clear();
}


void class_worker::add_dos() {
   if (help.MPI_COMM_HELP == MPI_COMM_NULL) {
       timer::add_dos = 0;
       for (int i = 0; i < random_walk.size(); i++) {
           histogram(random_walk[i].E_idx, random_walk[i].M_idx) += 1;
           dos      (random_walk[i].E_idx, random_walk[i].M_idx) += lnf;
       }
       random_walk.clear();
   }
}


void class_worker::add_hist_volume() {
    timer::add_hist_volume = 0;
    if (flag_one_over_t == 0) {
        t_check_convergence.tic();
        math::subtract_min_nonzero_one(histogram);
        saturation.push_back(histogram.sum());
        t_check_convergence.toc();
    }
}




void class_worker::check_saturation() {
    timer::check_saturation = 0;
    if (flag_one_over_t == 0) {
        int idx_to = (int) saturation.size() - 1;
        int idx_from = (int) (constants::check_saturation_from * idx_to);
        if (saturation.empty() || idx_to == idx_from || need_to_resize_global) {
            slope = 0;
            return;
        }
        vector<double> sat_double(saturation.begin() + idx_from, saturation.end()); //Cast to double
        ArrayXd Y = Map<ArrayXd>(sat_double.data(), sat_double.size());             //Cast to eigen array
        ArrayXd X = ArrayXd::LinSpaced(Y.size(), idx_from + 1, idx_to);
        slope = ((X - X.mean()).cwiseProduct(Y - Y.mean())).sum() / fmax(1, (X - X.mean()).cwiseAbs2().sum());
        if (slope < 0) {
            next_WL_iteration();
        }
        if (lnf < 1.0 / counter::MCS) {
            lnf = 1.0 / counter::MCS;
            flag_one_over_t = 1;
        }
    } else {
        counter::walks = (int) (-log(lnf) / log(2));
    }
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



