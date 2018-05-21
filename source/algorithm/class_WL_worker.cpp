//
// Created by david on 2016-07-24.
//
#include "class_WL_worker.h"

#define profiling_total                	1
#define profiling_print                	1
#define profiling_sweep                	1
#define profiling_swap                 	1
#define profiling_merge                	1
#define profiling_setup_team           	0
#define profiling_sync_team            	0
#define profiling_divide_range        	0
#define profiling_check_convergence	    0
#define profiling_make_MC_trial 		0
#define profiling_acceptance_criterion 	0
#define debug_sweep                     0

using namespace std;
using namespace Eigen;
int constants::team_size;
int constants::num_teams;
int counter::MCS;
int counter::walks;
int counter::swaps;
int counter::swap_accepts;
int counter::merges;
int timer::increment;
int timer::add_hist_volume;
int timer::check_saturation;
int timer::check_finish_line;
int timer::backup;
int timer::print;
int timer::swap;
int timer::sync_team;
int timer::setup_team;
int timer::divide_range;
int timer::sampling;

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
                                t_setup_team           (profiling_setup_team,           3,"t_setup" ),
                                t_sync_team            (profiling_sync_team,            3,"t_sync" ),
                                t_divide_range         (profiling_divide_range,         3,"t_divr" ),
                                t_check_convergence    (profiling_check_convergence,    3,"t_conv") ,
                                t_acceptance_criterion (profiling_acceptance_criterion, 3,"t_accr"),
                                team(id)
{
    rn::rng.seed((unsigned long)world_ID);
    lnf = 1.0;
    find_initial_limits();
    set_initial_local_bins();
    need_to_resize_global = 0;
    update_global_range();
    start_counters();
    rate_increment = 1;
    state_is_valid = false;
    trial_is_valid = false;
    cout << "ID: " << world_ID << " Started OK"<<endl;
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

void class_worker::start_counters(){
    counter::MCS                = 0;
    counter::walks              = 0;
    counter::swaps              = 0;
    counter::swap_accepts       = 0;
    counter::merges             = 0;
    timer::increment            = 0;
    timer::add_hist_volume      = 0;
    timer::check_saturation     = 0;
    timer::check_finish_line    = 0;
    timer::backup               = 0;
    timer::print                = 0;
    timer::swap 				= 0;
    timer::sync_team        	= 0;
    timer::setup_team			= 0;
    timer::divide_range         = 0;
    timer::sampling             = 0;
    flag_one_over_t             = 0;
}

void class_worker::rewind_timers(){
    timer::increment            = 0;
    timer::add_hist_volume      = 0;// math::mod(counter::MCS, constants::rate_add_hist_volume  );
    timer::check_saturation     = 0;// math::mod(counter::MCS, constants::rate_check_saturation );
    timer::check_finish_line    = 0;// math::mod(counter::MCS, constants::rate_check_finish_line);
    timer::backup               = 0;// math::mod(counter::MCS, constants::rate_backup_data      );
    timer::print                = 0;// math::mod(counter::MCS, constants::rate_print_status     );
    timer::swap 				= 0;// math::mod(counter::MCS, constants::rate_swap             );
    timer::sync_team        	= 0;// math::mod(counter::MCS, constants::rate_sync_team        );
    timer::setup_team			= 0;// math::mod(counter::MCS, constants::rate_setup_team       );
    timer::divide_range         = 0;// math::mod(counter::MCS, constants::rate_divide_range     );
    timer::sampling             = 0;
}

void class_worker::set_initial_local_bins(){
    random_walk.clear();
    switch(constants::rw_dims){
        case 1:
            E_bins << E_set;
            M_bins  = ArrayXd::Zero(1);
            histogram.conservativeResizeLike(ArrayXXi::Zero(E_bins.size(),M_bins.size()));
            dos.conservativeResizeLike(ArrayXXd::Zero(E_bins.size(),M_bins.size()));
            break;
        case 2:
            E_bins << E_set;
            M_bins << M_set;
            histogram.conservativeResizeLike(ArrayXXi::Zero(E_bins.size(),M_bins.size()));
            dos.conservativeResizeLike(ArrayXXd::Zero(E_bins.size(),M_bins.size()));
            break;
        default:
            cout << "Error in class_worker::set_initial_local_bins. Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            MPI_Finalize();
            exit(1);
    }
}

void class_worker::sweep(){
    t_sweep.tic();
    for (int i = 0; i < constants::N ; i++){
        model.make_new_state(E,M, E_trial, M_trial);
        update_global_range();
        acceptance_criterion();
        if(accept){
            E                               = E_trial;
            M                               = M_trial;
            model.flip();
            if (trial_is_valid) {
                E_idx                       = E_idx_trial;
                M_idx                       = M_idx_trial;
                if (++timer::increment >= rate_increment) {
                    timer::increment = 0;
                    random_walk.push_back({E_idx,M_idx});
                }
            }
        }
        else{
            if (state_is_valid) {
                if (++timer::increment >= rate_increment) {
                    timer::increment        = 0;
                    random_walk.push_back({E_idx,M_idx});

                }
            }
        }
    }

    if (flag_one_over_t) {
        lnf = 1.0 / max(1, counter::MCS);
    } else {
        flag_one_over_t = lnf < 1.0 / max(1, counter::MCS) ? 1 : 0;
    }
    t_sweep.toc();

}

void class_worker::acceptance_criterion(){
    t_acceptance_criterion.tic();
    state_in_window     = check_in_window(E);
    trial_in_window     = check_in_window(E_trial);
    if(state_is_valid){
        E_idx_trial         = math::binary_search_exact(E_bins, E_trial);
        M_idx_trial         = constants::rw_dims == 1 ? 0 : math::binary_search_exact(M_bins, M_trial);
    }else{
        E_idx               = math::binary_search_exact(E_bins, E);
        E_idx_trial         = math::binary_search_exact(E_bins, E_trial);
        M_idx               = constants::rw_dims == 1 ? 0 : math::binary_search_exact(M_bins, M);
        M_idx_trial         = constants::rw_dims == 1 ? 0 : math::binary_search_exact(M_bins, M_trial);
    }
    state_is_valid      = E_idx       != -1 && M_idx       != -1;
    trial_is_valid      = E_idx_trial != -1 && M_idx_trial != -1;
    if (!need_to_resize_global && state_is_valid && E_idx >= E_bins.size()){
        cout << "E_idx out of bounds" << endl;
    }
    if (!need_to_resize_global && state_is_valid && E != E_bins(E_idx)){
        cout << "Big error E: "<<E << " E: " << E_bins(E_idx) << endl;
    }

    bool normal_step = need_to_resize_global == 0 && state_is_valid; //Most often
    if(normal_step){
        if (trial_is_valid){
            accept = rn::uniform_double_1() < fmin(1, exp(dos(E_idx, M_idx) - dos(E_idx_trial, M_idx_trial)));
        }else{
            accept = false;
        }
        t_acceptance_criterion.toc();
        return;
    }

    bool go_away = need_to_resize_global == 1 || (state_in_window && trial_in_window && !trial_is_valid);
    bool go_home = !normal_step && !go_away;
//    if(world_ID==0) {
//        std::cout << " ID: " << world_ID
//                  << " go_away: " << go_away
//                  << " go_home: " << go_home
//                << " need_to_resize_global " << need_to_resize_global
//                << " state_in_window " << state_in_window
//                << " trial_in_window " << trial_in_window
//                << " trial_is_valid " << trial_is_valid
//                << " E_idx " << E_idx
//                << " M_idx " << M_idx
//                << " E_idx_trial " << E_idx_trial
//                << " M_idx_trial " << M_idx_trial
//                  << std::endl;
//    }
    if(go_away){
//        need_to_resize_global   = 1;
        E_set.insert(E_trial);
        M_set.insert(M_trial);
        walk_away_from_window();
        t_acceptance_criterion.toc();
        return;
    }

    if(go_home){
        walk_towards_window();
        t_acceptance_criterion.toc();
        return;
    }
    cout << "Undefined behavior! " << endl;
    exit(1);
    t_acceptance_criterion.toc();
}

void class_worker::update_global_range(){
    if (E_trial < E_min_global){
        E_min_global = E_trial;
        need_to_resize_global = 1;
    }
    if (E_trial > E_max_global){
        E_max_global  = E_trial;
        need_to_resize_global = 1;
    }
    if (M_trial < M_min_global and  constants::rw_dims == 2){
        M_min_global = M_trial;
        M_min_local = M_trial;
        need_to_resize_global = 1;
    }
    if (M_trial > M_max_global and  constants::rw_dims == 2){
        M_max_global = M_trial;
        M_max_local  = M_trial;
        need_to_resize_global = 1;
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

void class_worker::next_WL_iteration(){
    lnf = fmax(1e-12, lnf*constants::reduce_factor_lnf);
    histogram.setZero();
    saturation.clear();
    random_walk.clear();
    counter::walks++;
}

void class_worker::prev_WL_iteration(){
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

void class_worker::rewind_to_lowest_walk(){
    int min_walks;
    MPI_Allreduce(&counter::walks, &min_walks, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    counter::walks = min_walks;
    lnf = pow(constants::reduce_factor_lnf, min_walks);
    finish_line = 0;
    counter::MCS            = (int) (1.0/lnf);
    random_walk.clear();
    saturation.clear();
}

void class_worker::rewind_to_zero(){
    counter::walks = 0;
    lnf = 1;
    int save_vol_merges  = counter::merges;
    start_counters();
    counter::merges  = save_vol_merges;
    finish_line = 0;
    dos.setZero();
    histogram.setZero();
    random_walk.clear();
    saturation.clear();
    rewind_timers();
    flag_one_over_t = 0;
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
    if (counter::merges == constants::max_merges){
        rate_increment = max(1, (int)std::sqrt(fmax(dos_width,dos_height)));
    }else{
        rate_increment = 1;
    }
}

void class_worker::add_hist_volume(){
    timer::add_hist_volume = 0;
    if (flag_one_over_t == 0) {
        t_check_convergence.tic();
        math::subtract_min_nonzero_one(histogram);
        saturation.push_back(histogram.sum());
        t_check_convergence.toc();
    }
}

void class_worker::check_saturation(){
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



