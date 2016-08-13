//
// Created by david on 2016-07-24.
//
#include <iomanip>
#include "constants.h"
#include "counters_timers.h"
#include "class_data.h"
#include "wanglandau.h"
#include "MPI_algorithms.h"
#include <thread>
#include <chrono>
#define debug_all 0
#define debug_check_global_limits 1
using namespace std;

void WangLandau(class_worker &worker){
    int finish_line = 0;
    outdata out(worker.world_ID);
    while(finish_line == 0){
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_all) {
            cout << "Starting MCS " << counter::MCS;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));

        }

        sweep(worker);
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_all) {
            cout << " Sweep ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));

        }

        mpi::swap(worker);
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_all) {
            cout << "Swap ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));

        }

        check_global_limits(worker);
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_all) {
            cout << "check_limit ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));

        }

        check_convergence(worker, finish_line);
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_all) {
            cout << "check_conv ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));

        }

        print_status(worker);
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_all) {
            cout << "print ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));

        }

        MPI_Barrier(MPI_COMM_WORLD);
        //backup_data(worker,out);
        if (worker.world_ID == 0 && debug_all) {
            cout << "backup ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }

        MPI_Barrier(MPI_COMM_WORLD);
        divide_range(worker);
        if (worker.world_ID == 0 && debug_all) {
            cout << "divide ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }


        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_all) {
            cout << "OK "<< endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }

    }
    out.write_data_worker(worker);
    mpi::merge(worker);
    out.write_data_master(worker);


}

void sweep(class_worker &worker){
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
}

void check_convergence(class_worker &worker, int &finish_line){
    switch(worker.flag_one_over_t){
        case 0:
            add_hist_volume(worker);
            check_saturation(worker);
            break;
        case 1:
            check_one_over_t(worker);
            break;
        default:
            cout << "Error: check_convergence has wrong flag" << endl;
            exit(2);
    }
    if (timer::check_finish_line > constants::rate_check_finish_line) {
        timer::check_finish_line = 0;
        MPI_Allreduce(&worker.finish_line, &finish_line, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    }else{
        timer::check_finish_line++;
    }
}

void check_global_limits(class_worker &worker){
    if (timer::check_limits >= constants::rate_check_limits) {
        timer::check_limits = 0;
        MPI_Allreduce(&worker.need_to_resize_global, &worker.need_to_resize_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (worker.need_to_resize_global == 1) {
            worker.resize_global_range();
            worker.divide_global_range();
            worker.resize_local_bins();
            worker.prev_WL_iteration();
            worker.need_to_resize_global = 0;
        }
    }else{
        timer::check_limits++;
    }
}

void divide_range(class_worker &worker){
    if (timer::split_windows >= constants::rate_split_windows) {
        timer::split_windows = 0;
        int min_walks, need_to_resize;
        MPI_Allreduce(&counter::walks, &min_walks, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&worker.need_to_resize_global, &need_to_resize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (min_walks > 2 && need_to_resize == 0) {
            mpi::merge(worker);
            mpi::divide_global_range_dos_volume(worker);
            worker.rewind_to_lowest_walk();
        }
    }else{
        timer::split_windows++;
    }
}

void add_hist_volume(class_worker &worker) {
    if (timer::add_hist_volume > constants::rate_add_hist_volume) {
        timer::add_hist_volume = 0;
        //Find minimum histogram entry larger than 2.
        math::subtract_min_nonzero(worker.histogram);

        //Resize saturation vector
        if (worker.saturation.size() <= counter::saturation) {
            worker.saturation.conservativeResize(2 * std::max(1,counter::saturation));
            worker.saturation.bottomRows(counter::saturation).fill(0);
        }
        worker.saturation(counter::saturation) = worker.histogram.sum();
        counter::saturation++;

    }else{
        timer::add_hist_volume++;
    }
}

void check_saturation(class_worker &worker) {
    if (timer::check_saturation >= constants::rate_check_saturation) {
        timer::check_saturation = 0;
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
            if (worker.lnf < constants::minimum_lnf){
                worker.finish_line = 1;
            }
            if(worker.lnf < 1/counter::MCS){
            //if(counter::MCS > 1e5){
                worker.lnf = 1.0/counter::MCS;
                worker.flag_one_over_t = 1;         //Change to 1/t algorithm
            }
        }
    }else{
        timer::check_saturation++;
    }
}

void check_one_over_t (class_worker &worker){
    worker.lnf = 1.0/counter::MCS;
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

void print_status(class_worker &worker){
    if (timer::print >= constants::rate_print_status){
        timer::print = 0;
        for (int i = 0; i < worker.world_size; i++){
            if(worker.world_ID == i){
                if (worker.world_ID == 0){cout << endl;}
                cout    << "ID: "        << left << setw(3) << i
                        << " Walk: "  << left << setw(3) << counter::walks
                        << " f: "     << left << setw(16)<< fixed << setprecision(12) << exp(worker.lnf)
                        << " Bins: [" << left << setw(4) << worker.dos.rows() << " " << worker.dos.cols() << "]"
                        << " dE: "    << left << setw(7) << setprecision(2)   << worker.E_max_local - worker.E_min_local
                        << " E : ["   << left << setw(7) << setprecision(1)   << worker.E_min_local << " " << left << setw(7) << setprecision(1) << worker.E_max_local << "]"
                        << " Sw: "    << left << setw(5) << counter::swap_accepts
                        << " iw: "    << worker.in_window
                        << " 1/t: "   << worker.flag_one_over_t
                        << " Fin: "   << worker.finish_line
                        << " MCS: "   << left << setw(10) << counter::MCS
                        << endl;
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
                    <<  " MaxWalks: "<< fixed << setprecision(0) << ceil(log(constants::minimum_lnf)/log(constants::reduce_factor_lnf))
                    << "  -----"
                    << endl;
        }

    }else{
        timer::print++;
    }

}
