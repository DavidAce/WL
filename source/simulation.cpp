//
// Created by david on 2016-07-24.
//
#include "simulation.h"
#include "class_WL_read_data.h"


#define debug_sweep                     0
#define debug_convergence               0
#define debug_divide_range              0
#define debug_status                    0


using namespace std;

void do_simulations(class_worker &worker){
    //Begin by finding the global energy range in which to work.
//    worker.iteration = -1;
//    find_global_range(worker);
    for (int i = 0; i < constants::simulation_reps; i++){
        worker.iteration = i;
        worker.rewind_to_zero();
        wanglandau(worker);
    }
}


void wanglandau(class_worker &worker){
    int finish_line = 0;
    outdata out;
    class_backup backup;
    out.create_iteration_folder_master(worker.iteration, worker.world_ID);
    worker.t_total.tic();
    worker.t_print.tic();
    while(finish_line == 0){
        sweep(worker);
        if (timer::add_dos              >= constants::rate_add_dos          ){worker.add_dos         ()                       ;}
        if (timer::swap                 >= constants::rate_swap             ){mpi::swap              (worker)                 ;}
        if (timer::add_hist_volume      >= constants::rate_add_hist_volume  ){worker.add_hist_volume ()                       ;}
        if (timer::check_saturation     >= constants::rate_check_saturation ){worker.check_saturation()                       ;}
        if (timer::check_finish_line    >= constants::rate_check_finish_line){check_finish_line      (worker,out, finish_line);}
        if (timer::divide_range         >= constants::rate_divide_range     ){divide_range           (worker, backup,out)     ;}
//        if (timer::check_help           >= constants::rate_check_help       ){mpi::check_help        (worker)                 ;}
        if (timer::take_help            >= constants::rate_take_help        ){mpi::take_help         (worker)                 ;}
        if (timer::setup_help           >= constants::rate_setup_help       ){mpi::setup_help        (worker,backup)          ;}
        if (timer::print                >= constants::rate_print_status     ){print_status           (worker,false)           ;}

        counter::MCS++;
        timer::add_dos++;
        timer::add_hist_volume++;
        timer::check_finish_line++;
        timer::check_saturation++;
        timer::backup++;
        timer::print++;
        timer::swap++;
//        timer::check_help++;
        timer::take_help++;
        timer::setup_help++;
        timer::divide_range++;

    }
    print_status           (worker,true);
    backup.restore_state   (worker) ;
    out.write_data_worker  (worker) ;
    mpi::merge             (worker,false,true) ;
    out.write_data_master  (worker) ;
}

void sweep(class_worker &worker){
    if(debug_sweep){debug_print(worker,"Starting MCS " + std::to_string(counter::MCS));}
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
    if (worker.flag_one_over_t) {
        worker.lnf = 1.0 / counter::MCS;
    } else {
        worker.flag_one_over_t = worker.lnf < 1.0 / max(1, counter::MCS) ? 1 : 0;
    }
    worker.t_sweep.toc();
}

void check_finish_line(class_worker &worker, outdata &out, int &finish_line){
        timer::check_finish_line = 0;
        if (debug_convergence) { debug_print(worker, "Check Finish line "); }
        if (!worker.help.giving_help){
            if (worker.lnf < constants::minimum_lnf){
                worker.finish_line = 1;
            }
        }
        MPI_Allreduce(&worker.finish_line, &finish_line, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
}

void divide_range(class_worker &worker, class_backup &backup, outdata &out) {
    //Three situations, and they can only happen IF nobody is helping out or is finished.
    //1) We need to resize global range.
    //2) We have done far too few walks to do a dos_volume resize. Then we do a dos_area instead, only if everybody is in window!!
    //3) We have done enough walks to do a dos_volume resize
    timer::divide_range = 0;
    worker.t_divide_range.tic();
    int need_to_resize;
    MPI_Allreduce(&worker.need_to_resize_global, &need_to_resize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

//    if (min_walks > 6 && counter::vol_merges == 4){
//        out.write_data_worker(worker);
//        MPI_Barrier(MPI_COMM_WORLD);
//        exit(0);
//    }
    if (need_to_resize) {
        if (debug_divide_range) { debug_print(worker, "Dividing global energy\n"); }
        worker.resize_global_range();
        worker.divide_global_range_uniformly();
        worker.synchronize_sets();
        worker.adjust_local_bins();
        worker.need_to_resize_global = 0;
        worker.state_is_valid = false;
        counter::vol_merges = 0;
        worker.set_rate_increment();
        worker.prev_WL_iteration();
        worker.t_divide_range.tic();
        print_status(worker, true);

        return;
    }
    if (counter::vol_merges < constants::max_vol_merges) {
        int all_in_window;
        int min_walks;
        MPI_Allreduce(&counter::walks, &min_walks, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&worker.state_in_window, &all_in_window, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (all_in_window && min_walks >= counter::vol_merges) {
            //divide dos vol
            if (worker.world_ID == 0) { cout << "Dividing according to dos VOLUME" << endl; }
            backup.restore_state(worker);
            worker.help.reset();
            mpi::merge(worker, true, false);
            mpi::divide_global_range_dos_volume(worker);
            worker.rewind_to_lowest_walk();
            worker.set_rate_increment();
            counter::vol_merges++;
            if (counter::vol_merges == constants::max_vol_merges){
                worker.rewind_to_zero();
            }
            worker.state_is_valid = false;
            print_status(worker, true);
            worker.t_divide_range.toc();
            return;
        }
    }

}

void backup_to_file(class_worker &worker, outdata &out){
    if (!worker.help.giving_help) {
        if (timer::backup >= constants::rate_backup_data) {
            timer::backup = 0;
            int all_in_window;
            int need_to_resize;
            MPI_Allreduce(&worker.state_in_window, &all_in_window, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&worker.need_to_resize_global, &need_to_resize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if (need_to_resize == 0 && all_in_window == 1 && counter::vol_merges > 0) {
                mpi::merge(worker,false,false);
                out.write_data_master(worker);
            }
        }
    }
}

void print_status(class_worker &worker, bool force) {
    timer::print = 0;
    worker.t_total.toc();
    worker.t_print.toc();
    for (int i = 0; i < worker.world_size; i++){
        if(worker.world_ID == i){
            if (worker.world_ID == 0){cout << endl;}
            cout << fixed << showpoint;
            cout    << "ID: "        << left << setw(3) << i
                    << " Walk: "  << left << setw(3) << counter::walks
                    << " lnf: "     << left << setw(16)<< fixed << setprecision(12) << worker.lnf
                    << " Bins: [" << left << setw(4) << worker.dos.rows() << " " << worker.dos.cols() << "]";
                    if(debug_status){
               cout << " E: "     << left << setw(9) << setprecision(2)   << worker.E
                    << " M: "     << left << setw(9) << setprecision(2)   << worker.M;
               cout << " E_tr: "     << left << setw(9) << setprecision(2)   << worker.E_trial
                    << " M_tr: "     << left << setw(9) << setprecision(2)   << worker.M_trial;
                    }
               cout << " dE: "    << left << setw(7) << setprecision(2)   << worker.E_max_local - worker.E_min_local
                    << " : ["     << left << setw(7) << setprecision(1)   << worker.E_bins(0) << " " << left << setw(7) << setprecision(1) << worker.E_bins(worker.E_bins.size()-1) << "]"
                    << " Sw: "    << left << setw(7) << counter::swap_accepts
                    << " H: "     << right<< setw(3) <<  worker.help.helping_id
                    << " I: "     << left << setw(5) << worker.rate_increment
                    << " iw: "    << worker.state_in_window
                    << " NR: "    << worker.need_to_resize_global
                    << " 1/t: "   << worker.flag_one_over_t
                    << " Fin: "   << worker.finish_line
                    << " slope "  << left << setw(10) << worker.slope;
                    if(true) {
                cout<< " Edge dos: " << fixed << setprecision(3)
                    << left << setw(10) << worker.dos.topLeftCorner(1, 1).sum() << " "
                    << left << setw(10) << worker.dos.topRightCorner(1, 1).sum();
                    }

               cout << " MCS: "   << left << setw(10) << counter::MCS;
                worker.t_print               .print_delta();
                worker.t_sweep               .print_total_reset();
                worker.t_make_MC_trial       .print_total_reset();
                worker.t_acceptance_criterion.print_total_reset();
                worker.t_swap                .print_total_reset();
                worker.t_merge               .print_total_reset();
                worker.t_help_setup          .print_total_reset();
                worker.t_help                .print_total_reset();
                worker.t_divide_range        .print_total_reset<double,std::milli>();std::cout << "ms";
                worker.t_check_convergence   .print_total_reset<double,std::milli>();std::cout << "ms";
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
                << " MaxWalks: "    << fixed << setprecision(0) << (int) ceil(log(constants::minimum_lnf)/log(constants::reduce_factor_lnf))
                << " Merges: "      << fixed << setprecision(0) << counter::vol_merges  << "("<< constants::max_vol_merges  << ")"
                << " Iteration: "   << fixed << setprecision(0) << worker.iteration+1   << "("<< constants::simulation_reps << ")";
                worker.t_total.print_total<double>(); cout << " s";
                if(true){
                    cout    << " Edge dos: " << fixed << setprecision(3)
                            << worker.dos.topLeftCorner(1,1) << " "
                            << worker.dos.topRightCorner(1,1) << " " ;
//                                << worker.dos << endl << endl
//                                << worker.histogram << endl;
                }
                cout << "  -----"
                << endl;
    }
//        timer::print_tic = std::chrono::high_resolution_clock::now();
    worker.t_total.tic();
    worker.t_print.tic();
}
