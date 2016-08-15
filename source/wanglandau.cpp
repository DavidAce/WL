//
// Created by david on 2016-07-24.
//
#include <iomanip>
#include "constants.h"
#include "counters_timers.h"
#include "class_data.h"
#include "wanglandau.h"
#include "MPI_algorithms.h"
#include "class_profiling.h"
#include <thread>
#include <chrono>
#define debug_sweep                     0
#define debug_convergence               0
#define debug_global_limits             0
#define profiling_sweep                 1
#define profiling_swap                  1
#define profiling_check_global_limits   0
#define profiling_check_convergence     1
using namespace std;

void WangLandau(class_worker &worker){
    int finish_line = 0;
    class_profiling t_sweep(profiling_sweep),
                    t_swap(profiling_swap),
                    t_check_global_limits(profiling_check_global_limits),
                    t_check_convergence(profiling_check_convergence);

    outdata out(worker.world_ID);
    while(finish_line == 0){
        sweep(worker, t_sweep);
        mpi::swap(worker, t_swap);
        check_global_limits(worker, t_check_global_limits);
        check_convergence(worker, finish_line, t_check_convergence);
        print_status(worker, t_sweep, t_swap, t_check_global_limits, t_check_convergence);
        divide_range(worker);
        backup_data(worker,out);
    }
    out.write_data_worker(worker);
    mpi::merge(worker);
    out.write_data_master(worker);


}

void sweep(class_worker &worker, class_profiling &t_sweep){
    if (debug_sweep){
        if (worker.world_ID == 0) {
            cout << endl << "Starting MCS " << counter::MCS << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    t_sweep.tic();
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
    t_sweep.toc();
}

void check_convergence(class_worker &worker, int &finish_line, class_profiling &t_check_convergence){
    switch(worker.flag_one_over_t){
        case 0:
            if (timer::add_hist_volume > constants::rate_add_hist_volume) {
                timer::add_hist_volume = 0;
                if (debug_convergence){debug_print(worker,"Add hist volume ");}
                t_check_convergence.tic();
                add_hist_volume(worker);
                t_check_convergence.toc();

            }else{
                timer::add_hist_volume++;
            }
            if (timer::check_saturation >= constants::rate_check_saturation) {
                timer::check_saturation = 0;
                if (debug_convergence){debug_print(worker,"Check_saturation ");}
                t_check_convergence.tic();
                check_saturation(worker);
                t_check_convergence.toc();
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
            t_check_convergence.tic();
            check_one_over_t(worker);
            t_check_convergence.toc();
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

void check_global_limits(class_worker &worker, class_profiling & t_check_global_limits){
    timer::check_limits++;
    if (timer::check_limits >= constants::rate_check_limits) {
        timer::check_limits = 0;
        MPI_Allreduce(&worker.need_to_resize_global, &worker.need_to_resize_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (worker.need_to_resize_global == 1) {
            if (debug_global_limits){debug_print(worker,"Check Global Limits ");}
            worker.resize_global_range();
            worker.divide_global_range();
            worker.resize_local_bins();
            worker.prev_WL_iteration();
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
        if (min_walks > constants::min_walks && need_to_resize == 0 &&   counter::merges * counter::merges * constants::max_merges <  min_walks ) {
            mpi::merge(worker);
            mpi::divide_global_range_dos_volume(worker);
            worker.rewind_to_lowest_walk();
        }
    }
}

void add_hist_volume(class_worker &worker) {
    //Find minimum histogram entry larger than 2.
    math::subtract_min_nonzero(worker.histogram);
    //Resize saturation vector
    if (worker.saturation.size() <= counter::saturation) {
        worker.saturation.conservativeResize(2 * std::max(1, counter::saturation));
        worker.saturation.bottomRows(counter::saturation).fill(0);
    }
    worker.saturation(counter::saturation) = worker.histogram.sum();
    counter::saturation++;
}

void check_saturation(class_worker &worker) {
    int i, j;
    int idx = (int) (constants::check_saturation_from * counter::saturation);
    double Sx = 0, Sxy = 0, mX = 0, mY = 0;
    j = 0;
    //Compute means of the last 10%:
    for (i = idx; i < counter::saturation; i++) {
        mX += i;//X(i);
        mY += worker.saturation(i);
        j++;
    }
    mX /= j;
    mY /= j;
    for (i = idx; i < counter::saturation; i++) {
        Sx += pow(i - mX, 2);
        Sxy += (worker.saturation(i) - mY) * (i - mX);
    }
    double slope = Sxy / Sx;

    if (slope < 0) {
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

void print_status(class_worker &worker, class_profiling &t_sweep,
                                        class_profiling &t_swap,
                                        class_profiling &t_check_global_limits,
                                        class_profiling &t_check_convergence) {
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
                        << " dE: "    << left << setw(7) << setprecision(2)   << worker.E_max_local - worker.E_min_local
                        << " E : ["   << left << setw(7) << setprecision(1)   << worker.E_bins(0) << " " << left << setw(7) << setprecision(1) << worker.E_bins(worker.E_bins.size()-1) << "]"
                        << " Sw: "    << left << setw(5) << counter::swap_accepts
                        << " iw: "    << worker.in_window
                        << " 1/t: "   << worker.flag_one_over_t
                        << " Fin: "   << worker.finish_line
                        << " MCS: "   << left << setw(10) << counter::MCS
                        << " Time: "  << fixed << setprecision(3) << timer::elapsed_time_print.count() << " s ";
                if(profiling_sweep){
                    cout << " t_sweep = " << t_sweep;
                    t_sweep.reset();
                }
                if(profiling_swap){
                    cout << " t_swap = " << t_swap;
                    t_swap.reset();
                }
                if(profiling_check_global_limits){
                    cout << " t_glob = " << t_check_global_limits;
                    t_check_global_limits.reset();
                }
                if(profiling_check_convergence){
                    cout << " t_conv = " << t_check_convergence;
                    t_check_convergence.reset();
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
                    << " Total Time: " << fixed << setprecision(3) << timer::elapsed_time_total.count() << " s "
                    << "  -----"
                    << endl;
        }
        timer::print_tic = std::chrono::high_resolution_clock::now();
    }

}
