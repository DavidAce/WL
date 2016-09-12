//
// Created by david on 2016-07-24.
//
#include "simulation.h"
using namespace std;
#define debug_sweep                     0
#define debug_convergence               0
#define debug_global_limits             0
#define debug_saturation                0

void do_simulations(class_worker &worker){
    for (int i = 0; i < constants::simulation_reps; i++){
        worker.iteration = i;
        wanglandau(worker);
        worker.rewind_to_zero();
    }
}

void wanglandau(class_worker &worker){
    int finish_line = 0;
    outdata out(worker.world_ID,worker.iteration);
    while(finish_line == 0){
        sweep               (worker)              ;
        mpi::swap           (worker)              ;
        check_global_limits (worker)              ;
        check_convergence   (worker, finish_line) ;
        print_status        (worker)              ;
        divide_range        (worker)              ;
        backup_data         (worker,out)          ;
    }
    out.write_data_worker (worker) ;
    mpi::merge            (worker) ;
    out.write_data_master (worker) ;
}

void sweep(class_worker &worker){
    if (debug_sweep){
        if (worker.world_ID == 0) {
            cout << endl << "Starting MCS " << counter::MCS << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    worker.t_sweep.tic();
    for (int i = 0; i < constants::N ; i++){
        worker.make_MC_trial();
        worker.acceptance_criterion();
        if(worker.accept){
            worker.accept_MC_trial();

        }else{
            worker.reject_MC_trial();
        }
    }
    counter::MCS++;
    worker.t_sweep.toc();
}

void check_convergence(class_worker &worker, int &finish_line){
    switch(worker.flag_one_over_t){
        case 0:
            if (timer::add_hist_volume > constants::rate_add_hist_volume) {
                timer::add_hist_volume = 0;
                if (debug_convergence){debug_print(worker,"Add hist volume ");}
                worker.t_check_convergence.tic();
                add_hist_volume(worker);
                worker.t_check_convergence.toc();

            }else{
                timer::add_hist_volume++;
            }
            if (timer::check_saturation >= constants::rate_check_saturation) {
                timer::check_saturation = 0;
                if (debug_convergence){debug_print(worker,"Check_saturation ");}
                worker.t_check_convergence.tic();
                check_saturation(worker);
                worker.t_check_convergence.toc();
            }else{
                timer::check_saturation++;
            }
            if(worker.lnf < pow(constants::one_over_t_factor/counter::MCS, constants::one_over_t_exponent)){
                worker.lnf =  pow(constants::one_over_t_factor/counter::MCS, constants::one_over_t_exponent);
                worker.flag_one_over_t = 1;         //Change to 1/t algorithm
            }
            break;
        case 1:
            if (debug_convergence){debug_print(worker,"Check one over t ");}
            worker.t_check_convergence.tic();
            check_one_over_t(worker);
            worker.t_check_convergence.toc();
            break;
        default:
            cout << "Error: check_convergence has wrong flag" << endl;
            exit(2);
    }
    timer::check_finish_line++;
    if (timer::check_finish_line > constants::rate_check_finish_line) {
        timer::check_finish_line = 0;
        if (debug_convergence){debug_print(worker,"Check Finish line ");}
        MPI_Allreduce(&worker.finish_line, &finish_line, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    }
}

void check_global_limits(class_worker &worker){
    timer::check_limits++;
    if (timer::check_limits >= constants::rate_check_limits) {
        timer::check_limits = 0;
        MPI_Allreduce(MPI_IN_PLACE, &worker.need_to_resize_global,  1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (worker.need_to_resize_global == 1) {
            if (debug_global_limits){debug_print(worker,"Check Global Limits ");}
            worker.resize_global_range();
            worker.divide_global_range_energy();
            worker.resize_local_bins();
            worker.prev_WL_iteration();
            worker.in_window = worker.check_in_window(worker.E);
            if (worker.in_window){
                worker.find_current_state();
            }
            worker.need_to_resize_global = 0;
        }
    }
}

void divide_range(class_worker &worker){
    timer::split_windows++;
    if (timer::split_windows >= constants::rate_split_windows) {
        timer::split_windows = 0;
        int min_walks, need_to_resize;
        MPI_Allreduce(&counter::walks, &min_walks, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&worker.need_to_resize_global, &need_to_resize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (min_walks > constants::min_walks && need_to_resize == 0 &&   counter::merges< constants::max_merges ) {
            mpi::merge(worker);
            mpi::broadcast (worker) ;
            mpi::divide_global_range_dos_volume(worker);
            worker.rewind_to_lowest_walk();
        }
    }
}

void add_hist_volume(class_worker &worker) {
    //Subtract the smallest positive number (subtracts away regions impossible to reach as well)
    math::subtract_min_nonzero_one(worker.histogram);
    worker.saturation.push_back(worker.histogram.sum());


}

void check_saturation(class_worker &worker) {
    int i, j;
    //counter::saturation tells how many elements are in worker.saturation
    int idx = (int) (constants::check_saturation_from * (worker.saturation.size()-1));
    double Sx = 0, Sxy = 0, mX = 0, mY = 0;
    j = 0;
    //Compute means of the last 10%:
    for (i = idx; i < worker.saturation.size(); i++) {
        mX += i;//X(i);
        mY += worker.saturation[i];
        j++;
    }
    mX /= fmax(j,1);
    mY /= fmax(j,1);
    for (i = idx; i < worker.saturation.size(); i++) {
        Sx += pow(i - mX, 2);
        Sxy += (worker.saturation[i] - mY) * (i - mX);
    }
    worker.slope = Sxy / fmax(Sx,1);
    if (debug_saturation && counter::merges == 1 ) {
        cout << setprecision(0);
        cout << "ID " << worker.world_ID << ": "
             << worker.histogram.sum() << " "
             << idx  << " "
             << worker.saturation.size() << " "
             << worker.histogram.rows() <<  " x " << worker.histogram.cols() << " "
             << worker.slope << "  "
             << worker.saturation
             << endl;
        cout << flush;
        std::this_thread::sleep_for(std::chrono::microseconds(10));

    }
    math::subtract_min_nonzero(worker.histogram);
    if (worker.slope < 0) {
        worker.next_WL_iteration();
        if (worker.lnf < constants::minimum_lnf) {
            worker.finish_line = 1;
        }
    }
}

void check_one_over_t (class_worker &worker){
    worker.lnf =  pow(constants::one_over_t_factor/counter::MCS, constants::one_over_t_exponent);
    if (worker.lnf < constants::minimum_lnf){
        worker.finish_line = 1;
    }
}

void backup_data(class_worker &worker, outdata &out){
    if(timer::backup > constants::rate_backup_data) {
       timer::backup = 0;
       int all_in_window;
       int need_to_resize;
        MPI_Allreduce(&worker.in_window     , &all_in_window,  1, MPI_INT, MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(&worker.need_to_resize_global, &need_to_resize, 1, MPI_INT, MPI_MAX,MPI_COMM_WORLD);
       if (need_to_resize == 0 && all_in_window == 1){
           mpi::merge(worker);
           out.write_data_master(worker);
       }
    }else{
        timer::backup++;
    }
}

void debug_print(class_worker & worker, string input){
    if (worker.world_ID == 0) {
        cout << input;
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(10));
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void print_status(class_worker &worker) {
    timer::print++;
    if (timer::print >= constants::rate_print_status){
        timer::print = 0;
        timer::total_toc = std::chrono::high_resolution_clock::now();
        timer::print_toc = std::chrono::high_resolution_clock::now();
        timer::elapsed_time_print = timer::print_toc - timer::print_tic;
        timer::elapsed_time_total = timer::total_toc - timer::total_tic;
        for (int i = 0; i < worker.world_size; i++){
            if(worker.world_ID == i){
                if (worker.world_ID == 0){cout << endl;}
                cout    << "ID: "        << left << setw(3) << i
                        << " Walk: "  << left << setw(3) << counter::walks
                        << " f: "     << left << setw(16)<< fixed << setprecision(12) << exp(worker.lnf)
                        << " Bins: [" << left << setw(4) << worker.dos.rows() << " " << worker.dos.cols() << "]"
                        << " E: "     << left << setw(9) << setprecision(2)   << worker.E
                        << " M: "     << left << setw(9) << setprecision(2)   << worker.M
                        << " dE: "    << left << setw(7) << setprecision(2)   << worker.E_max_local - worker.E_min_local
                        << " : ["     << left << setw(7) << setprecision(1)   << worker.E_bins(0) << " " << left << setw(7) << setprecision(1) << worker.E_bins(worker.E_bins.size()-1) << "]"
                        << " Sw: "    << left << setw(5) << counter::swap_accepts
                        << " iw: "    << worker.in_window
                        << " NR: "    << worker.need_to_resize_global
                        << " 1/t: "   << worker.flag_one_over_t
                        << " Fin: "   << worker.finish_line
                        << " slope "  << left << setw(10) << worker.slope
                        << " MCS: "   << left << setw(10) << counter::MCS
                        << " Time: "  << fixed << setprecision(3) << timer::elapsed_time_print.count() << " s ";
                if(profiling_sweep){
                    cout << " t_sweep = " << worker.t_sweep;
                    worker.t_sweep.reset();
                }
                if(profiling_swap){
                    cout << " t_swap = " << worker.t_swap;
                    worker.t_swap.reset();
                }
                if(profiling_check_global_limits){
                    cout << " t_glob = " << worker.t_check_global_limits;
                    worker.t_check_global_limits.reset();
                }
                if(profiling_check_convergence){
                    cout << " t_conv = " << worker.t_check_convergence;
                    worker.t_check_convergence.reset();
                }
                if(profiling_make_MC_trial){
                    cout << " t_mkMC = " << worker.t_make_MC_trial;
                    worker.t_make_MC_trial.reset();
                }
                if(profiling_acceptance_criterion){
                    cout << " t_accr = " << worker.t_acceptance_criterion;
                    worker.t_acceptance_criterion.reset();
                }
                
                

                cout << endl;
            }
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));
            MPI_Barrier(MPI_COMM_WORLD);
        }
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(10));
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0){
            cout    << "-----"
                    << " MaxWalks: "   << fixed << setprecision(0) << ceil(log(constants::minimum_lnf)/log(constants::reduce_factor_lnf))
                    << " Iteration: "   << fixed << setprecision(0) << worker.iteration
                    << " Total Time: " << fixed << setprecision(3) << timer::elapsed_time_total.count() << " s "
                    << "  -----"
                    << endl;
        }
        timer::print_tic = std::chrono::high_resolution_clock::now();
    }

}
