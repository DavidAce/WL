//
// Created by david on 2016-07-24.
//
#include "simulation.h"
#include "algorithm/nmspc_WL_parallelization.h"
#include "class_WL_teams.h"
#include "class_WL_worker.h"
#include "general/class_tic_toc.h"
#include "IO/class_WL_print_data.h"
#include "nmspc_WL_counters_timers.h"
#include "params/nmspc_WL_constants.h"
#include <chrono>
#include <iomanip>
#include <numeric>
#include <ratio>
#include <thread>
#define debug_wanglandau 0
#define debug_check_finish_line 0
#define debug_divide_range 0
#define debug_status 0

using namespace std;

void do_simulations(class_worker &worker) {
    // Begin by finding the global energy range in which to work.
    for(int i = 0; i < constants::simulation_reps; i++) {
        worker.iteration = i;
        worker.rewind_to_zero();
        wanglandau(worker);
    }
}

void do_sampling(class_worker &worker) {
    if(constants::collect_samples) {
        outdata out;
        out.create_set_folder("outdata/samples/" + to_string(worker.world_ID) + "/");
        int samples = 0;
        worker.debug_print_all<1>("Starting sampling.\n");
        while(samples < constants::samples_to_collect) {
            worker.sweep();
            if(timer::swap >= constants::rate_swap) { parallel::swap2(worker); }
            if(timer::print >= constants::rate_print_status) { print_status(worker, false); }
            if(timer::sampling >= constants::rate_sampling) {
                timer::sampling = 0;
                out.write_sample(worker);
                samples++;
            }
            timer::print++;
            timer::swap++;
            timer::sampling++;
        }
        worker.debug_print_all<1>("Finished sampling.\n");
    }
}

void wanglandau(class_worker &worker) {
    int          finish_line = 0;
    outdata      out;
    class_backup backup;
    out.create_iteration_folder_commander(worker, worker.iteration);
    worker.t_total.tic();
    worker.t_print.tic();
    while(finish_line == 0) {
        worker.sweep();
        if(timer::sync_team >= constants::rate_sync_team) { sync_teams(worker); }
        if(timer::setup_team >= constants::rate_setup_team) { setup_teams(worker); }
        if(timer::swap >= constants::rate_swap) { parallel::swap2(worker); }
        if(timer::add_hist_volume >= constants::rate_add_hist_volume) { worker.add_hist_volume(); }
        if(timer::check_saturation >= constants::rate_check_saturation) { worker.check_saturation(); }
        if(timer::divide_range >= constants::rate_divide_range) { divide_range(worker, backup); }
        if(timer::print >= constants::rate_print_status) { print_status(worker, false); }
        if(timer::check_finish_line >= constants::rate_check_finish_line) { check_finish_line(worker, backup, finish_line); }

        counter::MCS++;
        timer::add_hist_volume++;
        timer::check_finish_line++;
        timer::check_saturation++;
        timer::backup++;
        timer::print++;
        timer::swap++;
        timer::sync_team++;
        timer::setup_team++;
        timer::divide_range++;
    }
    worker.debug_print_all<1>("MADE IT OUT ALIVE!\n");
    print_status(worker, true);
    MPI_Barrier(MPI_COMM_WORLD);

    backup.restore_state(worker);
    MPI_Barrier(MPI_COMM_WORLD);

    out.write_data_team_leader(worker);
    MPI_Barrier(MPI_COMM_WORLD);

    parallel::merge(worker, false, true);
    MPI_Barrier(MPI_COMM_WORLD);

    out.write_data_commander(worker);
}

void check_finish_line(class_worker &worker, class_backup &backup, int &finish_line) {
    timer::check_finish_line = 0;
    if(debug_check_finish_line) { worker.debug_print("Check Finish line "); }
    if(worker.lnf < constants::minimum_lnf) {
        worker.finish_line = 1;
        backup.backup_state(worker);
    }
    MPI_Allreduce(&worker.finish_line, &finish_line, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
}

void divide_range(class_worker &worker, class_backup &backup) {
    // Three situations, and they can only happen IF nobody is helping out or is finished.
    // 1) We need to resize global range.
    // 2) We have done far too few walks to do a dos_volume resize. Then we do a dos_area instead, only if everybody is in window!!
    // 3) We have done enough walks to do a dos_volume resize
    timer::divide_range = 0;
    worker.t_divide_range.tic();
    int need_to_resize;
    MPI_Allreduce(&worker.need_to_resize_global, &need_to_resize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if(need_to_resize) {
        if(debug_divide_range) { worker.debug_print("Dividing global energy\n"); }
        parallel::resize_global_range(worker);
        parallel::divide_global_range_uniform(worker);
        parallel::synchronize_sets(worker);
        parallel::adjust_local_bins(worker);
        worker.need_to_resize_global = 0;
        worker.state_is_valid        = false;
        counter::merges              = 0;
        worker.set_rate_increment();
        worker.prev_WL_iteration();
        worker.t_divide_range.tic();
        print_status(worker, true);
        return;
    }
    if(counter::merges < constants::max_merges) {
        int all_in_window;
        int min_walks;
        MPI_Allreduce(&counter::walks, &min_walks, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&worker.state_in_window, &all_in_window, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if(all_in_window && min_walks >= counter::merges) {
            if(worker.world_ID == 0) { cout << "Dividing according to dos" << endl; }

            backup.restore_state(worker);
            parallel::merge(worker, true, false);
            counter::merges++;
            if(counter::merges < constants::max_merges) {
                parallel::divide_global_range_dos_area(worker);
                worker.set_rate_increment();
                worker.rewind_to_lowest_walk();
            } else {
                parallel::divide_global_range_dos_volume(worker);
                worker.set_rate_increment();
                worker.rewind_to_zero();
            }
            worker.state_is_valid = false;
            print_status(worker, true);
            worker.t_divide_range.toc();
            return;
        }
    }
}

void sync_teams(class_worker &worker) {
    timer::sync_team = 0;
    if(worker.team->is_active()) {
        //        worker.debug_print_all<1>("size: " + to_string(worker.random_walk.size())+ "\n");
        worker.t_sync_team.tic();
        worker.team->sync_teams();
        counter::MCS += constants::rate_sync_team * (worker.team->get_team_size() - 1);
        timer::add_hist_volume += constants::rate_sync_team * (worker.team->get_team_size() - 1);
        timer::check_saturation += constants::rate_sync_team * (worker.team->get_team_size() - 1);
        worker.t_sync_team.toc();
    }
}

void setup_teams(class_worker &worker) {
    timer::setup_team = 0;
    worker.t_setup_team.tic();
    worker.team->setup_teams();
    worker.t_setup_team.toc();
}

void print_status(class_worker &worker, bool force) {
    timer::print = 0;
    worker.t_total.toc();
    worker.t_print.toc();
    for(int i = 0; i < worker.world_size; i++) {
        if(worker.world_ID == i) {
            if(worker.world_ID == 0) { cout << endl; }
            cout << fixed << showpoint;
            cout << "ID: " << left << setw(3) << i << " Walk: " << left << setw(3) << counter::walks << " lnf: " << left << setw(16) << fixed
                 << setprecision(12) << worker.lnf << " Bins: [" << left << setw(4) << worker.dos.rows() << " " << worker.dos.cols() << "]";
            if(debug_status) {
                cout << " E: " << left << setw(9) << setprecision(2) << worker.E << " M: " << left << setw(9) << setprecision(2) << worker.M;
                cout << " E_tr: " << left << setw(9) << setprecision(2) << worker.E_trial << " M_tr: " << left << setw(9) << setprecision(2) << worker.M_trial;
            }
            cout << " dE: " << left << setw(7) << setprecision(2) << worker.E_max_local - worker.E_min_local << " : [" << left << setw(7) << setprecision(1)
                 << worker.E_bins(0) << " " << left << setw(7) << setprecision(1) << worker.E_bins(worker.E_bins.size() - 1) << "]"
                 << " Sw: " << left << setw(7) << counter::swap_accepts << " Team: " << left << setw(2) << worker.team->get_team_id() << " I: " << left
                 << setw(3) << worker.rate_increment << " iw: " << worker.state_in_window << " NR: " << worker.need_to_resize_global
                 << " 1/t: " << worker.flag_one_over_t << " Fin: " << worker.finish_line << " slope " << left << setw(10) << worker.slope;
            cout << " MCS: " << left << setw(10) << counter::MCS;
            worker.t_print.print_delta();
            worker.t_sweep.print_total_reset();
            worker.t_acceptance_criterion.print_total_reset();
            worker.t_swap.print_total_reset();
            worker.t_merge.print_total_reset();
            worker.t_setup_team.print_total_reset();
            worker.t_sync_team.print_total_reset();
            worker.t_divide_range.print_total_reset<double, std::milli>();
            std::cout << "ms";
            worker.t_check_convergence.print_total_reset<double, std::milli>();
            std::cout << "ms";
            cout << endl;
        }
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(10));
        MPI_Barrier(MPI_COMM_WORLD);
    }
    cout.flush();
    std::this_thread::sleep_for(std::chrono::microseconds(10));
    MPI_Barrier(MPI_COMM_WORLD);
    if(worker.world_ID == 0) {
        cout << "-----"
             << " MaxWalks: " << fixed << setprecision(0) << (int) ceil(log(constants::minimum_lnf) / log(constants::reduce_factor_lnf)) << " Merges: " << fixed
             << setprecision(0) << counter::merges << "(" << constants::max_merges << ")"
             << " Iteration: " << fixed << setprecision(0) << worker.iteration + 1 << "(" << constants::simulation_reps << ")";
        worker.t_total.print_total<double>();
        cout << " s";
        if(debug_status) {
            cout << " Edge dos: " << fixed << setprecision(3) << worker.dos.topLeftCorner(1, 1) << " " << worker.dos.topRightCorner(1, 1) << " ";
            //                                << worker.dos << endl << endl
            //                                << worker.histogram << endl;
        }
        cout << "  -----" << endl;
    }
    //        timer::print_tic = std::chrono::high_resolution_clock::now();
    worker.t_total.tic();
    worker.t_print.tic();
}
