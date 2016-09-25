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
    for (int i = 0; i < constants::simulation_reps; i++){
        worker.iteration = i;
        wanglandau(worker);
        worker.rewind_to_zero();
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
        mpi::help           (worker,backup)       ;
        mpi::swap           (worker)              ;
        check_convergence   (worker, out,finish_line) ;
        divide_range        (worker, backup)              ;
//        backup_to_file         (worker,out)          ;
        print_status        (worker)              ;
    }
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

    counter::MCS++;
    if (worker.flag_one_over_t) {
            worker.lnf = 1.0 / counter::MCS;
    } else {
        if (worker.lnf < 1.0 / max(1, counter::MCS)) {
            worker.flag_one_over_t = 1;
        }else{
            worker.flag_one_over_t = 0;
        }
    }

    worker.t_sweep.toc();
}

void check_convergence(class_worker &worker, outdata &out, int &finish_line){
    if(!worker.help.giving_help) {
        if (debug_convergence) { cout << "ID: " << worker.world_ID << " Add hist volume "<< endl; }
        worker.add_hist_volume();
        if (debug_convergence) { cout << "ID: " << worker.world_ID << " Check Saturation "<< endl; }
        worker.check_saturation();
    }
    if (timer::check_finish_line > constants::rate_check_finish_line) {
        timer::check_finish_line = 0;
        if (debug_convergence) { debug_print(worker, "Check Finish line "); }
        if (!worker.help.giving_help){
            if (worker.lnf < constants::minimum_lnf){
                worker.finish_line = 1;
                worker.help.available = 1;
            }
        }
        MPI_Allreduce(&worker.finish_line, &finish_line, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    }else{
        timer::check_finish_line++;
    }

}


void divide_range(class_worker &worker, class_backup &backup){
    //Three situations, and they can only happen IF nobody is helping out or is finished.

    //1) We need to resize global range.
    //2) We have done far too few walks to do a dos_volume resize. Then we do a dos_area instead, only if everybody is in window!!
    //3) We have done enough walks to do a dos_volume resize
    if (timer::divide_range >= constants::rate_divide_range) {
        timer::divide_range = 0;
        worker.t_divide_range.tic();
        int any_helping = worker.help.getting_help;
        int all_in_window, min_walks, need_to_resize;
        int in_window = worker.in_window;
        MPI_Allreduce(MPI_IN_PLACE, &any_helping, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (any_helping == 0){
            MPI_Allreduce(&worker.need_to_resize_global,&need_to_resize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&counter::walks, &min_walks, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&in_window, &all_in_window, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
//            cout << "NR: " << worker.need_to_resize_global
//                 << " MW: " << min_walks << "(" << counter::walks << ")"
//                 << " AIW:" << all_in_window <<  "(" << in_window << ")" << endl;
            if (need_to_resize){
                if(debug_divide_range){debug_print(worker,"Dividing global energy\n");}
                //divide global energy
                worker.resize_global_range();
                worker.divide_global_range_energy();
                worker.resize_local_bins();
                worker.prev_WL_iteration();
                worker.dos.fill(0);
                worker.need_to_resize_global = 0;
                worker.find_current_state();
                worker.P_increment = 1.0 / sqrt(worker.E_bins.size());

//                worker.P_increment = 1.0 / sqrt(math::count_num_elements(worker.dos));

            }else if(min_walks < constants::min_walks && counter::merges < constants::max_merges && all_in_window == 1){
                //divide dos area
                if(worker.world_ID == 0){cout << "Dividing dos area" << endl;}
                mpi::merge(worker,true,false);
                mpi::divide_global_range_dos_area(worker);
                worker.P_increment = 1.0 / sqrt(worker.E_bins.size());

//                worker.P_increment = 1.0 / sqrt(math::count_num_elements(worker.dos));

            }else if (counter::merges < constants::max_merges && all_in_window == 1){
                //divide dos vol
                if(worker.world_ID == 0){cout << "Dividing dos vol" << endl;}
                //If anybody had started to help they need to be restored
                backup.restore_state(worker);
                worker.help.reset();
                mpi::merge(worker,true,false);
                mpi::divide_global_range_dos_volume(worker);
                counter::merges++;
                worker.P_increment = 1.0 / sqrt(worker.E_bins.size());
//                worker.P_increment = 1.0 / sqrt(math::count_num_elements(worker.dos));
                worker.rewind_to_zero();

            }

        }
        worker.t_divide_range.toc();
    }else{
        timer::divide_range++;
    }

}


//void add_hist_volume(class_worker &worker) {
//    //Subtract the smallest positive number plus one
//    math::subtract_min_nonzero_one(worker.histogram);
//    worker.saturation.push_back(worker.histogram.sum());
//}

//void check_saturation(class_worker &worker) {
//    int i, j;
//    //counter::saturation tells how many elements are in worker.saturation
//    int idx_to   = (int) worker.saturation.size()-1;
//    int idx_from = (int) (constants::check_saturation_from * idx_to);
//    double Sx = 0, Sxy = 0;
//    double mX, mY;
//    //Compute means of the last 10%:
//    mY = std::accumulate(worker.saturation.begin()+idx_from, worker.saturation.end(), 0);
//    mX = (idx_to + idx_from)*(idx_to - idx_from + 1)/2;
//    mX /= fmax(idx_to-idx_from,1);
//    mY /= fmax(idx_to-idx_from,1);
//
//    for (i = idx_from; i <= idx_to; i++) {
//        Sx  += (i-mX)*(i-mX) ;//pow(i - mX, 2);
//        Sxy += (worker.saturation[i] - mY) * (i - mX);
//    }
//    worker.slope = Sxy / fmax(Sx,1);
//    if (worker.slope < 0) {
//        worker.next_WL_iteration();
//        if (worker.lnf < constants::minimum_lnf) {
//            worker.finish_line = 1;
//        }
//    }
//}

void backup_to_file(class_worker &worker, outdata &out){
    if (!worker.help.giving_help) {
        if (timer::backup > constants::rate_backup_data) {
            timer::backup = 0;
            int all_in_window;
            int need_to_resize;
            MPI_Allreduce(&worker.in_window, &all_in_window, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&worker.need_to_resize_global, &need_to_resize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if (need_to_resize == 0 && all_in_window == 1 && counter::merges > 0) {
                mpi::merge(worker,false,false);
                out.write_data_master(worker);
            }
        } else {
            timer::backup++;
        }
    }
}

void print_status(class_worker &worker) {
    if (timer::print >= constants::rate_print_status){
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
                        << " Sw: "    << left << setw(5) << counter::swap_accepts
                        << " H: "     << right<< setw(3) <<  worker.help.helping_id
                        << " P: "     << left << setw(5) << 1/worker.P_increment
                        << " iw: "    << worker.in_window
                        << " NR: "    << worker.need_to_resize_global
                        << " 1/t: "   << worker.flag_one_over_t
                        << " Fin: "   << worker.finish_line;
                        if(debug_status){
                   cout << " slope "  << left << setw(10) << worker.slope;
                        }
                   cout << " MCS: "   << left << setw(10) << counter::MCS;
                    worker.t_print               .print_delta();
                    worker.t_sweep               .print_total_reset();
                    worker.t_make_MC_trial       .print_total_reset();
                    worker.t_acceptance_criterion.print_total_reset();
                    worker.t_swap                .print_total_reset();
                    worker.t_help                .print_total_reset();
                    worker.t_divide_range        .print_total_reset<double,std::milli>();
                    worker.t_check_convergence   .print_total_reset<double,std::milli>();
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
                    << " Iteration: "   << fixed << setprecision(0) << worker.iteration;
                    worker.t_total.print_total<double>(); cout << " s";
                    if(debug_status){
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
    }else{
        timer::print++;
    }
}
