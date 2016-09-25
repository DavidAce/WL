//
// Created by david on 2016-08-11.
//

#include "nmspc_WL_parallel_algorithms.h"
#define debug_swap      0
#define debug_merge     0
#define debug_bcast     0
#define debug_divide    0
#define debug_take_help 0
#define debug_setup_help 0
using namespace std;

namespace mpi {
    void swap(class_worker &worker) {
        //Use MPI Tag in the 100-200 range
        if (timer::swap > constants::rate_swap) {
            timer::swap = 0;
            int abort;
            worker.t_swap.tic();
            MPI_Allreduce(&worker.need_to_resize_global, &abort, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if (abort) { worker.t_swap.toc();return; }
            counter::swaps++;
            int swap, copy;
            double dos_X, dos_Y;
            double E_X, E_Y, M_X, M_Y;
            int E_X_idx, E_Y_idx, M_X_idx, M_Y_idx;
            int in_window_up;
            double E_min_up, E_max_up;
            double P_swap;      //Swap probability
            bool myTurn = math::mod(worker.world_ID, 2) == math::mod(counter::swaps, 2);

            int up = math::mod(worker.world_ID + 1, worker.world_size);
            int dn = math::mod(worker.world_ID - 1, worker.world_size);
            if (debug_swap) {
                for (int w = 0; w < worker.world_size; w++) {
                    if (w == worker.world_ID) {
                        cout << "ID: " << w << " Starting Swap. Myturn = " << myTurn << endl;
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }

            if (worker.in_window != worker.check_in_window(worker.E)){
                cout << "Worker is not really in window! Swap failed!" << endl;
                exit(1);
            }
            //Send current E and M to neighbors up and down. Receive X from below, Y from above.
            MPI_Sendrecv(&worker.E, 1, MPI_DOUBLE, up, 100, &E_X, 1, MPI_DOUBLE, dn, 100, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.M, 1, MPI_DOUBLE, up, 101, &M_X, 1, MPI_DOUBLE, dn, 101, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.E, 1, MPI_DOUBLE, dn, 102, &E_Y, 1, MPI_DOUBLE, up, 102, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.M, 1, MPI_DOUBLE, dn, 103, &M_Y, 1, MPI_DOUBLE, up, 103, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
            //Check if the neighbors position is within my overlap region. If so, find the indices.
            MPI_Sendrecv(&worker.E_min_local, 1, MPI_DOUBLE, dn, 104, &E_min_up, 1, MPI_DOUBLE, up, 104, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.E_max_local, 1, MPI_DOUBLE, dn, 105, &E_max_up, 1, MPI_DOUBLE, up, 105, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
//            MPI_Sendrecv(&worker.E_max_local, 1, MPI_DOUBLE, up, 105, &E_max_dn, 1, MPI_DOUBLE, dn, 105, MPI_COMM_WORLD,
//                         MPI_STATUS_IGNORE);
//            MPI_Sendrecv(&worker.in_window, 1, MPI_INT, up, 1055, &in_window_dn, 1, MPI_INT, dn, 1055, MPI_COMM_WORLD,
//                         MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.in_window, 1, MPI_INT, dn, 1044, &in_window_up, 1, MPI_INT, up, 1044, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
            //Now both swappees need to know if it is ok to go ahead with a swap.
            if (debug_swap) {
                for (int w = 0; w < worker.world_size; w++) {
                    if (w == worker.world_ID) {
                        cout << "ID: " << w << " SendRecv Successful"
                             << " E_X = " << E_X
                             << " E_Y = " << E_Y
                             << " E_idx = " << worker.E_idx
                             << " E_idx_trial = " << worker.E_idx_trial
                             << endl;
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
            int go_ahead;
            if (myTurn) {
                //Go ahead if inside the upper window and vice versa
                go_ahead = worker.E >= E_min_up && E_Y <= worker.E_max_local;
                go_ahead = worker.E <= E_max_up && E_Y >= worker.E_min_local && go_ahead;
                copy = !go_ahead && !worker.in_window && ((worker.E < E_Y && worker.E < worker.E_min_local) ||
                                                          (worker.E > E_Y && worker.E > worker.E_max_local));
                copy = copy && !worker.help.giving_help;
                //Make sure the last worker doesn't swap!
                go_ahead = go_ahead && worker.world_ID != worker.world_size - 1;
                MPI_Send(&go_ahead, 1, MPI_INT, up, 106, MPI_COMM_WORLD);
                MPI_Send(&copy, 1, MPI_INT, up, 107, MPI_COMM_WORLD);
            } else {
                MPI_Recv(&go_ahead, 1, MPI_INT, dn, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&copy, 1, MPI_INT, dn, 107, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if (debug_swap) {
                for (int w = 0; w < worker.world_size; w++) {
                    if (w == worker.world_ID) {
                        cout << "ID: " << w << " Received goahead = " << go_ahead << " Copy = " << copy << endl;
//                        cout << worker << endl;
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
            if (myTurn) {
                if (go_ahead) {
                    MPI_Recv(&dos_X, 1, MPI_DOUBLE, up, 108, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&dos_Y, 1, MPI_DOUBLE, up, 109, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    E_Y_idx = math::binary_search(worker.E_bins, E_Y);
                    M_Y_idx = math::binary_search(worker.M_bins, M_Y);
                    P_swap = exp(worker.dos(worker.E_idx, worker.M_idx)
                                 - worker.dos(E_Y_idx, M_Y_idx)
                                 + dos_Y
                                 - dos_X);
                    if (rn::uniform_double_1() < fmin(1, P_swap)) {
                        swap = 1;
                    } else {
                        swap = 0;
                    }
                } else {
                    swap = 0;
                }
                MPI_Send(&swap, 1, MPI_INT, up, 110, MPI_COMM_WORLD);
            } else {
                E_X_idx = math::binary_search(worker.E_bins, E_X);
                M_X_idx = math::binary_search(worker.M_bins, M_X);
                MPI_Send(&worker.dos(E_X_idx, M_X_idx), 1, MPI_DOUBLE, dn, 108, MPI_COMM_WORLD);
                MPI_Send(&worker.dos(worker.E_idx, worker.M_idx), 1, MPI_DOUBLE, dn, 109, MPI_COMM_WORLD);
                MPI_Recv(&swap, 1, MPI_INT, dn, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (debug_swap) {
                for (int w = 0; w < worker.world_size; w++) {
                    if (w == worker.world_ID) {
                        cout << "ID: " << w << " Swap = " << swap  << endl;
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
            if(debug_swap){debug_print(worker,"Starting swap\n");}
            //Now do the swapping if you got lucky
            if (myTurn) {
                if (swap == 1) {
                    MPI_Sendrecv_replace(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, up,
                                         111, up, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    worker.E = E_Y;
                    worker.M = M_Y;
                    worker.find_current_state();
                } else if (copy == 1) {
                    MPI_Recv(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, up, 112,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    worker.E = E_Y;
                    worker.M = M_Y;
                    worker.find_current_state();
                }
            } else {
                if (swap == 1) {
                    MPI_Sendrecv_replace(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, dn,
                                         111, dn, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    worker.E = E_X;
                    worker.M = M_X;
                    worker.find_current_state();
                } else if (copy == 1) {
                    MPI_Send(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, dn, 112,
                             MPI_COMM_WORLD);
                }
            }
            counter::swap_accepts += swap;
            if (debug_swap) {
                for (int w = 0; w < worker.world_size; w++) {
                    if (w == worker.world_ID) {
                        cout << "ID: " << w << " Swap OK " << endl;
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
            worker.t_swap.toc();
        } else {
            timer::swap++;
        }
    }

    void merge(class_worker &worker, bool broadcast, bool trim) {
        if(debug_merge){debug_print(worker,"\nMerging. ");}

        ArrayXXd dos_total,  dos_recv;
        ArrayXd E_total, M_total, E_recv, M_recv;
        int E_sizes[worker.world_size];
        int M_sizes[worker.world_size];
        int E_size = (int) worker.E_bins.size();
        int M_size = (int) worker.M_bins.size();

        double diff;

        MPI_Gather(&E_size, 1, MPI_INT, E_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&M_size, 1, MPI_INT, M_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(debug_merge){debug_print(worker," Gathering dos ");}

        if (worker.world_ID == 0) {
            dos_total = worker.dos;
            E_total = worker.E_bins;
            M_total = worker.M_bins;
        }
//        int from_total, from_up, rows_total, rows_up;
//        int E_merge_idx, E_merge_idx_up;
//        int M_merge_idx, M_merge_idx_up;
        int E_shared_low   , E_shared_high;
        int E_shared_low_up, E_shared_high_up;
        vector<ArrayXd> dos_merge;
        vector<double> E_merge;
        for (int w = 1; w < worker.world_size; w++) {
            if (worker.world_ID == 0) {
                if (debug_merge) {
                    cout << "Receiving from " << w << endl;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
//                cout << dos_total << endl<< endl;

                dos_recv.resize(E_sizes[w], M_sizes[w]);
                E_recv.resize(E_sizes[w]);
                M_recv.resize(M_sizes[w]);
                MPI_Recv(dos_recv.data(), E_sizes[w] * M_sizes[w], MPI_DOUBLE, w, 60, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                MPI_Recv(E_recv.data(), E_sizes[w], MPI_DOUBLE, w, 61, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(M_recv.data(), M_sizes[w], MPI_DOUBLE, w, 62, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


                //Find average height of shared sections
                E_shared_low    = math::binary_search(E_total, fmax(E_recv.minCoeff(), E_total.minCoeff()));
                E_shared_high   = math::binary_search(E_total, fmin(E_recv.maxCoeff(), E_total.maxCoeff()));
                E_shared_low_up = math::binary_search(E_recv , fmax(E_recv.minCoeff(), E_total.minCoeff()));
                E_shared_high_up= math::binary_search(E_recv , fmin(E_recv.maxCoeff(), E_total.maxCoeff()));

//                cout << " E:low     "<<E_shared_low       << endl;
//                cout << " E:high    "<< E_shared_high      << endl;
//                cout << " E:low up  "<< E_shared_low_up    << endl;
//                cout << " E:high up "<< E_shared_high_up   << endl;
//                cout << " diff      "<< math::nanzeromean(dos_total.middleRows(E_shared_low   , E_shared_high    - E_shared_low))<< endl;
//                cout << " diff up   "<< math::nanzeromean(dos_recv .middleRows(E_shared_low_up, E_shared_high_up - E_shared_low_up))<< endl;
                diff = math::nanzeromean(dos_total.middleRows(E_shared_low   , E_shared_high    - E_shared_low))
                      -math::nanzeromean(dos_recv .middleRows(E_shared_low_up, E_shared_high_up - E_shared_low_up));




//                //Find coordinates on dos_total
//                E_merge_idx = math::find_matching_slope(dos_total, dos_recv, E_total, E_recv, M_total, M_recv);
//                M_merge_idx = math::nanmaxCoeff_idx(dos_total.row(E_merge_idx));
//
//                dos_total.row(E_merge_idx).maxCoeff(&M_merge_idx);
//                //Find coordinates on received dos;
//                E_merge_idx_up = math::binary_search(E_recv, E_total(E_merge_idx));
//                M_merge_idx_up = math::nanmaxCoeff_idx(dos_recv.row(E_merge_idx_up));
//                dos_recv.row(E_merge_idx_up).maxCoeff(&M_merge_idx_up);
//                //Find difference between heights at these points
//                diff = dos_total(E_merge_idx, M_merge_idx) - dos_recv(E_merge_idx_up, M_merge_idx_up);
                //Add that difference to the next
                if (std::isnan(diff)) {
                    cout << "Tried to concatenate " << w - 1 << " and " << w << ". Diff between two DOS is NaN, exiting"
                         << endl;
                    cout << "dos_tot" << endl << dos_total << endl;
                    cout << "dos_rec" << endl << dos_recv << endl;
                    exit(1);
                }
                math::add_to_nonzero_nonnan(dos_recv, diff);
                //Now now all doses should be the same height
                if (debug_merge) {
                    cout << "Concatenating " << w - 1 << " and " << w << endl;
                     cout << "dos_tot" << endl << dos_total << endl;
                    cout << "dos_rec" << endl << dos_recv << endl;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }


                double weight;
                int j = 0;
                int rows_total      = E_shared_high     -  E_shared_low;
                int rows_total_up   = E_shared_high_up  -  E_shared_low_up;
                if (rows_total  != rows_total_up){
//                    cout << "Rows mismatch!!" << endl;
//                    cout << "Rows_total     = " << rows_total << endl;
//                    cout << "Rows_total_up  = " << rows_total_up << endl;
//                    cout << dos_total << endl << endl;
//                    cout << dos_recv << endl << endl << endl;
//                    for (int i = 0; i < E_total.size(); i++){
//                        if (i == E_shared_low){cout << "[";}
//                        cout << E_total(i)  <<" ";
//                        if (i == E_shared_high){cout << "]";}
//                    }
//                    cout << endl ;
//                    for (int i = 0; i < E_sizes[w]; i++){
//                        if (i == E_shared_low_up){cout << "[";}
//                        cout << E_recv(i)  <<" ";
//                        if (i == E_shared_high_up){cout << "]";}
//                    }
//                    cout << endl ;


                }
//                cout << "Starting pushback 1" << endl;
//                cout << "Rows: " << rows_total << endl;
                for (int i = 0; i < E_total.size() ; i++){
                    if (i < E_shared_low){
                        dos_merge.push_back(dos_total.row(i));
                        E_merge.push_back(E_total(i));
                    }else  if(i >= E_shared_low && i <= E_shared_high){
                        weight = (double) j / rows_total;
//                        cout << setprecision(5)<< "weight = " << weight << " j = " << j << endl;
                        int E_idx_up = math::binary_search(E_recv, E_total(i));
                        dos_merge.push_back((1-weight)*dos_total.row(i) + weight*dos_recv.row(E_idx_up) );
                        E_merge.push_back(E_total(i));
                        j++;
                    }
                }
//                cout << "Starting pushback 2" << endl;

                for (int i = 0; i < E_recv.size();i ++){
                    if (i > E_shared_high_up){
                        dos_merge.push_back(dos_recv.row(i) );
                        E_merge.push_back(E_recv(i));
                    }
                }
//                cout << "Starting Resize" << endl;
//                cout << "dos_merge.size() = " << dos_merge.size() << endl;
//                for (int i = 0; i < dos_merge.size(); i++){
//                    cout << dos_merge[i].transpose() << endl;
//
//                }
//                for (int i = 0; i < E_merge.size(); i++){
//                    cout << E_merge[i] << " ";
//
//                }
//                cout << endl ;

                dos_total.resize(dos_merge.size(), M_sizes[w]);
                E_total.resize(E_merge.size());
                for (unsigned int i = 0; i < dos_merge.size(); i++){
                    dos_total.row(i)    = dos_merge[i];
                    E_total(i)          = E_merge[i];
                }
//                cout << "Clearing" << endl;
                dos_merge.clear();
                E_merge.clear();
//                from_total  = 0;
//                from_up     = E_shared_low_up;
//                rows_total  = E_shared_high     - E_shared_low;
//                rows_up     = E_shared_high_up  - E_shared_low_up;
//                dos_temp.resize(dos_total.middleRows(from_total, rows_total).rows() +
//                                dos_recv.middleRows(from_up, rows_up).rows(), M_size);
//
//                //Compute the starting points and number of rows for concatenation
//                from_total = 0;
//                from_up = E_merge_idx_up + 1;
//
//                rows_total = E_merge_idx + 1;
//                rows_up = E_sizes[w] - E_merge_idx_up - 1;
//
//                E_temp.resize(E_total.segment(from_total, rows_total).size() + E_recv.segment(from_up, rows_up).size());
//                dos_temp.resize(dos_total.middleRows(from_total, rows_total).rows() +
//                                dos_recv.middleRows(from_up, rows_up).rows(), M_size);
//
//                dos_temp << dos_total.middleRows(from_total, rows_total),
//                        dos_recv.middleRows(from_up, rows_up);
//                E_temp << E_total.segment(from_total, rows_total),
//                        E_recv.segment(from_up, rows_up);
//                dos_total = dos_temp;
//                E_total = E_temp;

                if (debug_merge) {
                    cout << "OK ";
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
            } else if (w == worker.world_ID) {
                MPI_Send(worker.dos.data(), E_size * M_size, MPI_DOUBLE, 0, 60, MPI_COMM_WORLD);
                MPI_Send(worker.E_bins.data(), E_size, MPI_DOUBLE, 0, 61, MPI_COMM_WORLD);
                MPI_Send(worker.M_bins.data(), M_size, MPI_DOUBLE, 0, 62, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (worker.world_ID == 0) {
            //Let the master worker have the results
            worker.dos_total    = dos_total;
            worker.E_bins_total = E_total;
            worker.M_bins_total = M_total;
        }
        if (trim){
            if (worker.world_ID == 0) {
                //Trim away nan rows
                math::subtract_min_nonzero_nan(worker.dos_total);
                math::remove_nan_rows(worker.dos_total, worker.E_bins_total);
            }
        }
        if (broadcast){
            broadcast_merger(worker);
        }
    }


    void merge2(class_worker &worker, bool broadcast, bool trim) {
        if(debug_merge){debug_print(worker,"\nMerging. ");}

        ArrayXXd dos_total, dos_temp, dos_recv;
        ArrayXd E_total, M_total, E_temp, M_temp, E_recv, M_recv;
        int E_sizes[worker.world_size];
        int M_sizes[worker.world_size];
        int E_size = (int) worker.E_bins.size();
        int M_size = (int) worker.M_bins.size();

        double diff;

        MPI_Gather(&E_size, 1, MPI_INT, E_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&M_size, 1, MPI_INT, M_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(debug_merge){debug_print(worker," Gathering dos ");}

        if (worker.world_ID == 0) {
            dos_total = worker.dos;
            E_total = worker.E_bins;
            M_total = worker.M_bins;
        }
        int from_total, from_up, rows_total, rows_up;
        int E_merge_idx, E_merge_idx_up;
        int M_merge_idx, M_merge_idx_up;
        for (int w = 1; w < worker.world_size; w++) {
            if (worker.world_ID == 0) {
                if (debug_merge) {
                    cout << "Receiving from " << w << endl;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
                dos_recv.resize(E_sizes[w], M_sizes[w]);
                E_recv.resize(E_sizes[w]);
                M_recv.resize(M_sizes[w]);
                MPI_Recv(dos_recv.data(), E_sizes[w] * M_sizes[w], MPI_DOUBLE, w, 60, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                MPI_Recv(E_recv.data(), E_sizes[w], MPI_DOUBLE, w, 61, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(M_recv.data(), M_sizes[w], MPI_DOUBLE, w, 62, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //Find coordinates on dos_total
                E_merge_idx = math::find_matching_slope(dos_total, dos_recv, E_total, E_recv, M_total, M_recv);
                M_merge_idx = math::nanmaxCoeff_idx(dos_total.row(E_merge_idx));

//                dos_total.row(E_merge_idx).maxCoeff(&M_merge_idx);
                //Find coordinates on received dos;
                E_merge_idx_up = math::binary_search(E_recv, E_total(E_merge_idx));
                M_merge_idx_up = math::nanmaxCoeff_idx(dos_recv.row(E_merge_idx_up));
//                dos_recv.row(E_merge_idx_up).maxCoeff(&M_merge_idx_up);
                //Find difference between heights at these points
                diff = dos_total(E_merge_idx, M_merge_idx) - dos_recv(E_merge_idx_up, M_merge_idx_up);
                //Add that difference to the next
                if (std::isnan(diff)) {
                    cout << "Tried to concatenate " << w - 1 << " and " << w << ". Diff between two DOS is NaN, exiting"
                         << endl;
                    cout << "Merge_idx    = [" << E_merge_idx << ", " << M_merge_idx << "]" << endl;
                    cout << "Merge_idx_up = [" << E_merge_idx_up << ", " << M_merge_idx_up << "]" << endl;

                    cout << "E_bins_total   " << E_total(E_merge_idx) << " [" << E_merge_idx << "]" << "    = "
                         << E_total.transpose() << endl;
                    cout << "E_bins_up      " << E_recv(E_merge_idx_up) << " [" << E_merge_idx_up << "]" << "    = "
                         << E_recv.transpose() << endl;
                    cout << "dos_tot" << endl << dos_total << endl;
                    cout << "dos_rec" << endl << dos_recv << endl;
                    MPI_Finalize();
                    exit(1);
                }
                math::add_to_nonzero_nonnan(dos_recv, diff);
                //Now now all doses should be the same height
                if (debug_merge) {
                    cout << "Concatenating " << w - 1 << " and " << w << endl;
                    cout << "Merge_idx    = [" << E_merge_idx << ", " << M_merge_idx << "]" << endl;
                    cout << "Merge_idx_up = [" << E_merge_idx_up << ", " << M_merge_idx_up << "]" << endl;

                    cout << "E_bins_total   " << E_total(E_merge_idx) << " [" << E_merge_idx << "]" << "    = "
                         << E_total.transpose() << endl;
                    cout << "E_bins_up      " << E_recv(E_merge_idx_up) << " [" << E_merge_idx_up << "]" << "    = "
                         << E_recv.transpose() << endl;
                    cout << "dos_tot" << endl << dos_total << endl;
                    cout << "dos_rec" << endl << dos_recv << endl;

                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }

                //Compute the starting points and number of rows for concatenation
                from_total = 0;
                from_up = E_merge_idx_up + 1;

                rows_total = E_merge_idx + 1;
                rows_up = E_sizes[w] - E_merge_idx_up - 1;

                E_temp.resize(E_total.segment(from_total, rows_total).size() + E_recv.segment(from_up, rows_up).size());
                dos_temp.resize(dos_total.middleRows(from_total, rows_total).rows() +
                                dos_recv.middleRows(from_up, rows_up).rows(), M_size);

                dos_temp << dos_total.middleRows(from_total, rows_total),
                        dos_recv.middleRows(from_up, rows_up);
                E_temp << E_total.segment(from_total, rows_total),
                        E_recv.segment(from_up, rows_up);
                dos_total = dos_temp;
                E_total = E_temp;


                if (debug_merge) {
                    cout << "OK ";
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
            } else if (w == worker.world_ID) {
                MPI_Send(worker.dos.data(), E_size * M_size, MPI_DOUBLE, 0, 60, MPI_COMM_WORLD);
                MPI_Send(worker.E_bins.data(), E_size, MPI_DOUBLE, 0, 61, MPI_COMM_WORLD);
                MPI_Send(worker.M_bins.data(), M_size, MPI_DOUBLE, 0, 62, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (worker.world_ID == 0) {
            //Let the master worker have the results
            worker.dos_total    = dos_total;
            worker.E_bins_total = E_total;
            worker.M_bins_total = M_total;
        }
        if (trim){
            if (worker.world_ID == 0) {
                //Trim away nan rows
                math::subtract_min_nonzero_nan(worker.dos_total);
                math::remove_nan_rows(worker.dos_total, worker.E_bins_total);
            }
        }
        if (broadcast){
            broadcast_merger(worker);
        }
    }

    void broadcast_merger(class_worker &worker) {
        if(debug_bcast){debug_print(worker,"\nBroadcasting raw to workers ");}
        mpi::bcast_dynamic(worker.dos_total, MPI_DOUBLE, 0);
        mpi::bcast_dynamic(worker.E_bins_total, MPI_DOUBLE, 0);
        mpi::bcast_dynamic(worker.M_bins_total, MPI_DOUBLE, 0);
        if(debug_bcast){debug_print(worker,"Broadcast OK! \n ");}
    }

    void divide_global_range_dos_area(class_worker &worker) {
        //Update limits
        if(debug_divide){debug_print(worker,"\n Dividing Area. ");}
        double global_range = math::area(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
        double local_range = global_range / worker.world_size;
        double x = constants::overlap_factor_dos_area / (1 - constants::overlap_factor_dos_area / 2);
        double overlap_range = local_range * x;//2.0*(worker.world_size - 2.0 + x)/worker.world_size;
        //The overlap_range is the total range in a domain that will have overlap, either up or down.
        //Find the boundaries of the DOS domain that gives every worker the  same DOS volume to work on
        int E_min_local_idx, E_max_local_idx;
        if      (worker.world_ID == 0)                      { overlap_range = overlap_range / 2; }
        else if (worker.world_ID == worker.world_size - 1)  { overlap_range = overlap_range / 2; }
        else                                                { overlap_range = overlap_range / 4; }
        E_min_local_idx = math::area_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                           worker.world_ID * local_range - overlap_range);
        E_max_local_idx = math::area_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                           (worker.world_ID + 1) * local_range + overlap_range);

        while (E_max_local_idx - E_min_local_idx < constants::bins) {
            E_min_local_idx--;
            E_max_local_idx++;
            E_min_local_idx = max(E_min_local_idx, 0);
            E_max_local_idx = min(E_max_local_idx, (int) worker.E_bins_total.size() - 1);
        }

        worker.E_min_local = worker.E_bins_total(E_min_local_idx);
        worker.E_max_local = worker.E_bins_total(E_max_local_idx);
        worker.M_min_local = worker.M_min_global;
        worker.M_max_local = worker.M_max_global;
        if(debug_divide){debug_print(worker," Inheriting from dos_total. ");}

        //Inherit the corresponding part of the dos
        int from = E_min_local_idx;
        int rows = E_max_local_idx - E_min_local_idx + 1;

        if (from + rows > worker.E_bins_total.size()) {
            cout << "TOO MANY ROWS   |  from + rows = " << from + rows << " E_bins_total.size() = "
                 << worker.E_bins_total.size() << endl;
            cout << worker << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::seconds(1));
            MPI_Finalize();
            exit(15);
        }


        worker.dos          = worker.dos_total.middleRows(from, rows);
        worker.E_bins       = worker.E_bins_total.segment(from, rows);
        worker.histogram    = ArrayXXi::Zero(worker.dos.rows(), worker.dos.cols());
        worker.find_current_state();

        if (worker.model.discrete_model) {
            worker.E_set.clear();
            for (int i = 0; i < worker.E_bins.size(); i++) {
                worker.E_set.insert(worker.E_bins(i));
            }
        }

        if (debug_divide) {
            for (int w = 0; w < worker.world_size; w++) {
                if (w == worker.world_ID) {
                    cout << setprecision(2);
                    cout << "ID: " << w << " Bounds : " << worker.E_min_local << " " << worker.E_max_local << endl;
                    cout << worker.E_bins.transpose() << endl;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        if(debug_divide){debug_print(worker," Merged OK. \n ");}

    }

    void divide_global_range_dos_volume(class_worker &worker) {

        if(debug_divide){debug_print(worker," \nDividing dos volume ");}

        double global_range = math::volume(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
        double local_range = global_range / worker.world_size;
        double x = constants::overlap_factor_dos_vol / (1 - constants::overlap_factor_dos_vol / 2);
        double overlap_range = local_range * x;//2.0*(worker.world_size - 2.0 + x)/worker.world_size;
        //The overlap_range is the total range in a domain that will have overlap, either up or down.
        //Find the boundaries of the DOS domain that gives every worker the  same DOS volume to work on
        int E_min_local_idx, E_max_local_idx;

        if (worker.world_ID == 0) { overlap_range = overlap_range / 2; }
        else if (worker.world_ID == worker.world_size - 1) { overlap_range = overlap_range / 2; }
        else { overlap_range = overlap_range / 4; }
        E_min_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                           worker.world_ID * local_range - overlap_range);
        E_max_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                           (worker.world_ID + 1) * local_range + overlap_range);

        while (E_max_local_idx - E_min_local_idx < constants::bins) {
            E_min_local_idx--;
            E_max_local_idx++;
            E_min_local_idx = max(E_min_local_idx, 0);
            E_max_local_idx = min(E_max_local_idx, (int) worker.E_bins_total.size() - 1);
        }

        worker.E_min_local = worker.E_bins_total(E_min_local_idx);
        worker.E_max_local = worker.E_bins_total(E_max_local_idx);
        worker.M_min_local = worker.M_min_global;
        worker.M_max_local = worker.M_max_global;
        if(debug_divide){debug_print(worker," Inheriting from dos_total ");}

        //Inherit the corresponding part of the dos
        int from = E_min_local_idx;
        int rows = E_max_local_idx - E_min_local_idx + 1;

        if (from + rows > worker.E_bins_total.size()) {
            cout << "TOO MANY ROWS   |  from + rows = " << from + rows << " E_bins_total.size() = "
                 << worker.E_bins_total.size() << endl;
            cout << worker << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::seconds(1));
            MPI_Finalize();
            exit(15);
        }


        worker.dos          = worker.dos_total.middleRows(from, rows);
        worker.E_bins       = worker.E_bins_total.segment(from, rows);
        worker.histogram    = ArrayXXi::Zero(worker.dos.rows(), worker.dos.cols());
        worker.in_window    = worker.check_in_window(worker.E);
        worker.find_current_state();

        if (worker.model.discrete_model) {
            worker.E_set.clear();
            for (int i = 0; i < worker.E_bins.size(); i++) {
                worker.E_set.insert(worker.E_bins(i));
            }
        }

        if (debug_divide) {
            for (int w = 0; w < worker.world_size; w++) {
                if (w == worker.world_ID) {
                    cout << setprecision(2);
                    cout << "ID: " << w << " Bounds : " << worker.E_min_local << " " << worker.E_max_local << endl;
                    cout << worker.E_bins.transpose() << endl;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        if(debug_divide){debug_print(worker," Merge OK.\n ");}

    }

    void take_help(class_worker &worker) {
        //Offload help
        for (int w = 0; w < worker.world_size; w++) {
            if (worker.world_ID == worker.help.whos_helping_who(w)) {
                //You should receive help
                MPI_Recv(worker.help.histogram_recv.data(), (int) worker.help.histogram_recv.size(), MPI_INT, w, w, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (worker.flag_one_over_t == 0){
                    worker.histogram        += worker.help.histogram_recv;
                }
                worker.dos += worker.help.histogram_recv.cast<double>() * worker.lnf;
                counter::MCS                += constants::rate_take_help;
                timer::add_hist_volume      += constants::rate_take_help - 1;
                timer::check_saturation     += constants::rate_take_help - 1;
                worker.add_hist_volume();
                worker.check_saturation();

            } else if (worker.world_ID == w && worker.help.whos_helping_who(w) >= 0) {
                //You should send help
                MPI_Send(worker.histogram.data(), (int) worker.histogram.size(), MPI_INT, worker.help.whos_helping_who(w), w, MPI_COMM_WORLD);
                worker.histogram.fill(0);
            }
        }

        //Now get an updated dos:
        //If there has been no change since last time, simply share the dos and lnf and return;
        for (int w = 0; w < worker.world_size; w++) {
            if (worker.world_ID == worker.help.whos_helping_who(w)) {
                //You should send
                MPI_Send(worker.dos.data(), (int) worker.dos.size(), MPI_DOUBLE, w, w, MPI_COMM_WORLD);
                MPI_Send(&worker.lnf, 1, MPI_DOUBLE, w, w+worker.world_size, MPI_COMM_WORLD);
                MPI_Send(&counter::MCS, 1, MPI_INT, w, w+worker.world_size, MPI_COMM_WORLD);
            } else if (worker.world_ID == w && worker.help.whos_helping_who(w) >= 0) {
                //You should receive
                MPI_Recv(worker.dos.data(), (int) worker.dos.size(), MPI_DOUBLE, worker.help.whos_helping_who(w), w,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&worker.lnf, 1, MPI_DOUBLE, worker.help.whos_helping_who(w), w+worker.world_size, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&counter::MCS, 1, MPI_INT, worker.help.whos_helping_who(w), w+worker.world_size, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                worker.flag_one_over_t = worker.lnf < 1.0 / max(1, counter::MCS) ? 1 : 0;

            }
        }


    }

    void setup_help(class_worker &worker, class_backup &backup) {
        //Find out who's finished
        int any_finished = worker.finish_line;
        MPI_Allreduce(MPI_IN_PLACE, &any_finished, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        //if nobody is finished, return.
        if (any_finished == 0) { return; }
        //If some are finished, start by offloading any help that may be had already
        if(debug_setup_help){debug_print(worker,"Setting up help. ");}

        //Now all the help.available workers need to be set up again. Two scenarios, either they've (1) helped before,
        //or they've (2) never helped.
        //If (1), just take_help() and return.
        //If (2), setup from scratch.

        worker.help.available = worker.finish_line;
        worker.help.giving_help = false;
        worker.help.getting_help = false;
        worker.help.helping_id = -1;
        ArrayXi all_available(worker.world_size);
        MPI_Allgather(&worker.help.available, 1, MPI_INT, all_available.data(), 1, MPI_INT, MPI_COMM_WORLD);

        //Now find out who to help next. compute max amount of helpers per worker.
        int max_helpers = (int) ceil((double) all_available.sum() / (worker.world_size - all_available.sum()));
        ArrayXi given_help(worker.world_size);
        ArrayXi whos_helping_who_old = worker.help.whos_helping_who; //For comparison!
        given_help.fill(0);
        worker.help.whos_helping_who.fill(-1);
        for (int w = 0; w < worker.world_size; w++) {
            if (w == worker.world_ID) {
                if (worker.finish_line) {
                    for (int i = 0; i < worker.world_size; i++) {
                        if (all_available(i) == 0 && given_help(i) < max_helpers) {
                            given_help(i)++;
                            worker.help.available = 0;
                            worker.help.giving_help = true;
                            worker.help.helping_id = i;
                            break;
                        }
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(given_help.data(), worker.world_size, MPI_INT, w, MPI_COMM_WORLD);
        }
        worker.help.getting_help = given_help(worker.world_ID) > 0;
        MPI_Allgather(&worker.help.helping_id, 1, MPI_INT, worker.help.whos_helping_who.data(), 1, MPI_INT,
                      MPI_COMM_WORLD);
        //Now every available guy knows who to help, and the helpees know who to send their info to.

        //If there has been no change since last time, simply return;
        if (worker.help.whos_helping_who.cwiseEqual(whos_helping_who_old).all()) {
            if(debug_setup_help){debug_print(worker,"Same help settings\n");}
//            if(debug_setup_help){debug_print(worker,worker.help.whos_helping_who.transpose().eval());}
            return;
        }

        //Otherwise send out the details for help to begin.
        for (int w = 0; w < worker.world_size; w++) {
            if (worker.world_ID == worker.help.whos_helping_who(w)) {
                //You should send
                mpi::send_dynamic(worker.dos, MPI_DOUBLE, w);
                mpi::send_dynamic(worker.E_bins, MPI_DOUBLE, w);
                mpi::send_dynamic(worker.M_bins, MPI_DOUBLE, w);
                MPI_Send(&worker.lnf, 1, MPI_DOUBLE, w, 4, MPI_COMM_WORLD);
                MPI_Send(&worker.E, 1, MPI_DOUBLE, w, 5, MPI_COMM_WORLD);
                MPI_Send(&worker.M, 1, MPI_DOUBLE, w, 6, MPI_COMM_WORLD);
                MPI_Send(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, w,7, MPI_COMM_WORLD);
                MPI_Send(&counter::MCS, 1, MPI_INT, w, 8, MPI_COMM_WORLD);

                worker.help.histogram_recv.resizeLike(worker.histogram);
                worker.help.histogram_recv.fill(0);

            } else if (worker.world_ID == w && worker.help.whos_helping_who(w) >= 0) {
                //You should receive
                //Backup first
                backup.backup_state(worker);
                mpi::recv_dynamic(worker.dos, MPI_DOUBLE, worker.help.whos_helping_who(w));
                mpi::recv_dynamic(worker.E_bins, MPI_DOUBLE, worker.help.whos_helping_who(w));
                mpi::recv_dynamic(worker.M_bins, MPI_DOUBLE, worker.help.whos_helping_who(w));
                MPI_Recv(&worker.lnf, 1, MPI_DOUBLE, worker.help.whos_helping_who(w), 4, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
                MPI_Recv(&worker.E, 1, MPI_DOUBLE, worker.help.whos_helping_who(w), 5, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
                MPI_Recv(&worker.M, 1, MPI_DOUBLE, worker.help.whos_helping_who(w), 6, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
                MPI_Recv(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, worker.help.whos_helping_who(w),7, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&counter::MCS, 1, MPI_INT, worker.help.whos_helping_who(w), 8, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
                worker.flag_one_over_t = worker.lnf < 1.0 / max(1, counter::MCS) ? 1 : 0;
                //Set all values needed to sweep
                worker.E_min_local = worker.E_bins.minCoeff();
                worker.E_max_local = worker.E_bins.maxCoeff();
                worker.M_min_local = worker.M_bins.minCoeff();
                worker.M_max_local = worker.M_bins.maxCoeff();
                worker.in_window = worker.check_in_window(worker.E);
                worker.find_current_state();
//                worker.P_increment = 1.0 / sqrt(math::count_num_elements(worker.dos));
                worker.P_increment = 1.0 / sqrt(worker.E_bins.size());

                worker.histogram.resizeLike(worker.dos);
                worker.histogram.fill(0);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            //Now if you didn't get to help anybody out you might as well restore your own dos.
        }
        if (!worker.help.giving_help) {
            backup.restore_state(worker);
        }
        if(debug_setup_help){debug_print(worker,"New help settings\n");}
        if(debug_setup_help){debug_print(worker,worker.help.whos_helping_who.transpose().eval());}
//        cout << endl <<"ID: " << worker.world_ID << " In_window = " << worker.in_window << " [" << worker.E_min_local << " " << worker.E << " " << worker.E_max_local <<"]" <<endl;
    }

    void help(class_worker &worker, class_backup &backup) {
        if (timer::take_help > constants::rate_take_help) {
            timer::take_help = 0;
            if (worker.help.giving_help || worker.help.getting_help) {
                worker.t_help.tic();
                take_help(worker);
                worker.t_help.toc();
            }

        } else {
            timer::take_help++;
        }
        if (timer::setup_help > constants::rate_setup_help) {
            timer::setup_help = 0;
//            if (counter::merges >= constants::max_merges) {
                setup_help(worker, backup);
//            }
        } else {
            timer::setup_help++;
        }
    }

}















//
//
//
//
//
//
//
//    void help_out(class_worker &worker, class_backup &backup){
//        //With some time interval
//        timer::help_out++;
//        if (timer::help_out > constants::rate_help_out) {
//            timer::help_out = 0;
//            //Find out who's finished
//            int any_finished = worker.finish_line;
//            MPI_Allreduce(MPI_IN_PLACE, &any_finished,1, MPI_INT, MPI_MAX,MPI_COMM_WORLD );
//            //if nobody is finished, return.
//            if (any_finished == 0){return;}
//            if (debug_help_out){debug_print(worker,"Help has arrived!\n");}
//            //If some are finished, start by offloading any help that may be had already
//
//            for (int w = 0; w < worker.world_size; w++) {
//                if(worker.world_ID == worker.help.whos_helping_who(w)){
//                    //You should receive help
////                cout << "ID: " << worker.world_ID << " Receiving help from ID: "<< w <<endl;
//                    ArrayXXi histogram_temp(worker.histogram.rows(), worker.histogram.cols());
//                    mpi::recv_dynamic(histogram_temp, MPI_INT, w);
//                    if (histogram_temp.size() != worker.histogram.size() ){
//                        cout << "Size mismatch!!!!" << endl;
//                        exit(1);
//                    }
//                    worker.histogram += histogram_temp;
//
//                    for (int j = 0; j < worker.histogram.cols(); j++){
//                        for (int i = 0; i < worker.histogram.rows(); i++){
//                            worker.dos(i,j) += histogram_temp(i,j)*worker.lnf;
//                        }
//                    }
//                    counter::MCS += constants::rate_help_out;
//
//                }else if(worker.world_ID == w && worker.help.helping_out && worker.help.whos_helping_who(w) >= 0){
//                    //You should send help
////                cout << "ID: " << worker.world_ID << " Sending help to ID: "<< worker.help.whos_helping_who(w) << endl;
//                    math::subtract_min_nonzero_one(worker.histogram);
//                    mpi::send_dynamic(worker.histogram, MPI_INT, worker.help.whos_helping_who(w));
//                    //Relieve this hero of his duties
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//            }
//            if (debug_help_out){debug_print(worker,"Offloaded help \n");}
//
//
//            //Now all the help.available workers need to be set up again. Two scenarios, either they've (1) helped before,
//            //or they've (2) never helped.
//            //If (1), just take over the new updated dos.
//            //If (2), setup from scratch. (default)
//            worker.help.available = worker.finish_line;
//            worker.help.helping_out = false;
//            worker.help.helping_id = -1;
//            ArrayXi all_available(worker.world_size);
//            MPI_Allgather(&worker.help.available,1,MPI_INT, all_available.data(),1, MPI_INT,MPI_COMM_WORLD );
//
//            if (debug_help_out){debug_print(worker," New helpers available!\n");}
//
//            //Now find out who to help next. compute max amount of helpers per worker.
//            int max_helpers = (int)ceil((double)all_available.sum()/(worker.world_size - all_available.sum()));
//            ArrayXi given_help(worker.world_size);
//            ArrayXi whos_helping_who_old = worker.help.whos_helping_who;
//            given_help.fill(0);
//
//            worker.help.whos_helping_who.fill(-1);
//            for (int w = 0; w < worker.world_size; w++){
//                if (w == worker.world_ID){
//                    if (worker.finish_line){
//                        for(int i = 0; i < worker.world_size; i++){
//                            if (all_available(i) == 0 && given_help(i) < max_helpers){
//                                given_help(i)++;
//                                worker.help.available   = 0;
//                                worker.help.helping_out = true;
//                                worker.help.helping_id  = i;
//                                backup.backup_state(worker);
//                                break;
//                            }
//                        }
//                    }
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//                MPI_Bcast(given_help.data(),worker.world_size,MPI_INT, w, MPI_COMM_WORLD );
//            }
//            MPI_Allgather(&worker.help.helping_id, 1, MPI_INT, worker.help.whos_helping_who.data(),1, MPI_INT,MPI_COMM_WORLD);
//            //Now every available guy knows who to help, and the helpees know who to send their info to.
//
//            //If there has been no change since last time, simply share the dos and lnf and return;
//            if (worker.help.whos_helping_who.cwiseEqual(whos_helping_who_old).all()){
//                for (int w = 0; w < worker.world_size; w++) {
//                    if (worker.world_ID == worker.help.whos_helping_who(w)) {
//                        //You should send
////                    mpi::send_dynamic(worker.dos, MPI_DOUBLE, w);
//                        MPI_Send(worker.dos.data(), (int)worker.dos.size(), MPI_DOUBLE, w, 0, MPI_COMM_WORLD);
//                        MPI_Send(&worker.lnf, 1, MPI_DOUBLE, w, 1, MPI_COMM_WORLD);
//                    } else if (worker.world_ID == w && worker.help.whos_helping_who(w) >= 0) {
//                        //You should receive
//                        MPI_Recv(worker.dos.data(), (int)worker.dos.size(), MPI_DOUBLE, worker.help.whos_helping_who(w), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
////                    mpi::recv_dynamic(worker.dos, MPI_DOUBLE, worker.help.whos_helping_who(w));
//                        MPI_Recv(&worker.lnf, 1, MPI_DOUBLE, worker.help.whos_helping_who(w), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//                        //Set all values needed to sweep
//                        worker.histogram.fill(0);
////                cout << "ID: " << worker.world_ID << " Received Size "<< worker.histogram.size() << " from " << worker.help.whos_helping_who(w) << endl;
//                    }
//                }
//                return;
//            }
//            //Start sending out.
//            for (int w = 0; w < worker.world_size; w++) {
//                if(worker.world_ID == worker.help.whos_helping_who(w)){
//                    //You should send
//                    mpi::send_dynamic(worker.dos   , MPI_DOUBLE,w);
//                    mpi::send_dynamic(worker.E_bins, MPI_DOUBLE,w);
//                    mpi::send_dynamic(worker.M_bins, MPI_DOUBLE,w);
//                    MPI_Send(&worker.lnf, 1, MPI_DOUBLE, w, 4,MPI_COMM_WORLD);
////                cout << "ID: " << worker.world_ID << " Sent Size "<< worker.histogram.size() << " to " << worker.help.whos_helping_who(w) << endl;
//
//                }else if(worker.world_ID == w && worker.help.whos_helping_who(w) >= 0){
//                    //You should receive
//                    mpi::recv_dynamic(worker.dos   ,MPI_DOUBLE,worker.help.whos_helping_who(w));
//                    mpi::recv_dynamic(worker.E_bins,MPI_DOUBLE,worker.help.whos_helping_who(w));
//                    mpi::recv_dynamic(worker.M_bins,MPI_DOUBLE,worker.help.whos_helping_who(w));
//                    MPI_Recv(&worker.lnf, 1 , MPI_DOUBLE, worker.help.whos_helping_who(w), 4,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//                    //Set all values needed to sweep
//                    worker.E_min_local = worker.E_bins.minCoeff();
//                    worker.E_max_local = worker.E_bins.maxCoeff();
//                    worker.M_min_local = worker.M_bins.minCoeff();
//                    worker.M_max_local = worker.M_bins.maxCoeff();
//                    worker.in_window   = worker.check_in_window(worker.E);
//                    worker.P_increment = 1.0/sqrt(math::count_num_elements(worker.dos));
//                    worker.histogram.resizeLike(worker.dos);
//                    worker.histogram.fill(0);
////                cout << "ID: " << worker.world_ID << " Received Size "<< worker.histogram.size() << " from " << worker.help.whos_helping_who(w) << endl;
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//                //Now if you didn't get to help anybody out you might as well restore your own dos.
//                if (!worker.help.helping_out) {
//                    backup.restore_state(worker);
//                }
//            }
//        }else{
//            timer::help_out++;
//        }
//    }









//    void divide_global_range_dos_volume(class_worker &worker) {
//        //Update limits
//        if (worker.world_ID == 0 && debug_divide) {
//            cout << "Dividing. ";
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        }
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        math::subtract_min_nonnan(worker.dos_total);
//        if (worker.world_ID == 0){
//            cout << "dos size " << worker.dos_total.rows()    << " x " << worker.dos_total.cols() << endl;
//            cout << "E   size " << worker.E_bins_total.rows() << " x " << worker.E_bins_total.cols() << endl;
//        }
//        if (worker.world_ID == 0 && debug_divide) {
//            cout << "Computing Volume ";
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        }
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        double global_volume = math::volume(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
//        double local_volume = global_volume / worker.world_size;
//        double x = constants::overlap_factor_dos_vol / (1 - constants::overlap_factor_dos_vol / 2)  ;
//        //Add a little bit if there are too many workers (Add nothing if world_size == 2, and up to local_volume if world_size == inf)
//        double overlap_range = local_volume * x ;//2.0*(worker.world_size - 2.0 + x)/worker.world_size;
////        cout << local_volume / global_volume << overlap_range << " " << local_volume*x << endl << endl;
//        //Find the boundaries of the DOS domain that gives every worker the  same DOS volume to work on
//        int E_min_local_idx, E_max_local_idx;
////        if (worker.world_ID == 0){
////            cout << "dos volume " << global_volume << endl;
////        }
//        int min_width = 10;
//        if (worker.world_ID == 0) {
//            E_min_local_idx = 0;
//            E_max_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
//                                               (worker.world_ID + 1) * local_volume + overlap_range / 2);
//            while (E_max_local_idx - E_min_local_idx < min_width) { E_max_local_idx++; }
//
//        } else if (worker.world_ID == worker.world_size - 1) {
//            E_min_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
//                                               worker.world_ID * local_volume - overlap_range / 2);
//            E_max_local_idx = (int) worker.E_bins_total.size() - 1;
//            while (E_max_local_idx - E_min_local_idx < min_width) { E_min_local_idx--; }
//
//        } else {
//            E_min_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
//                                               worker.world_ID * local_volume - overlap_range / 4);
//            E_max_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
//                                               (worker.world_ID + 1) * local_volume + overlap_range / 4);
//            while (E_max_local_idx - E_min_local_idx < min_width) {
//                E_min_local_idx--;
//                E_max_local_idx++;
//            }
//        }
//        E_min_local_idx = max(E_min_local_idx, 0);
//        E_max_local_idx = min(E_max_local_idx, (int)worker.E_bins_total.size()-1);
//
//        worker.E_min_local = worker.E_bins_total(E_min_local_idx);
//        worker.E_max_local = worker.E_bins_total(E_max_local_idx);
//        worker.M_min_local = worker.M_min_global;
//        worker.M_max_local = worker.M_max_global;
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (worker.world_ID == 0 && debug_divide) {
//            cout << "...OK. Inherit from dos_total ";
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//
//        //Inherit the corresponding part of the dos
//        int from = E_min_local_idx;
//        int rows = E_max_local_idx - E_min_local_idx + 1;
//
//        if (from + rows > worker.E_bins_total.size()) {
//            cout << "TOO MANY ROWS   |  from + rows = " << from+rows << " E_bins_total.size() = " << worker.E_bins_total.size() << endl;
//            cout << worker << endl;
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::seconds(1));
//            MPI_Finalize();
//            exit(15);
//        }
//
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        if (E_max_local_idx <= E_min_local_idx) {
//            cout << "Local range backwards! " << endl;
//            cout << worker << endl;
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::seconds(1));
//            MPI_Finalize();
//            exit(16);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        worker.dos = worker.dos_total.middleRows(from, rows);
//        worker.E_bins = worker.E_bins_total.segment(from, rows);
//        worker.in_window = worker.check_in_window(worker.E);
//        worker.histogram = ArrayXXi::Zero(worker.dos.rows(), worker.dos.cols());
//
//        if (worker.in_window) {
//            worker.E_idx = math::binary_search(worker.E_bins, worker.E);
//            worker.M_idx = math::binary_search(worker.M_bins, worker.M);
//        }
//        if (worker.model.discrete_model) {
//            worker.E_set.clear();
//            for (int i = 0; i < worker.E_bins.size(); i++) {
//                worker.E_set.insert(worker.E_bins(i));
//            }
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        if (worker.world_ID == 0 && debug_divide) {
//            cout << "...OK " << endl;
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::microseconds(1000));
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (debug_divide) {
//            for (int w = 0; w < worker.world_size; w++) {
//                if (w == worker.world_ID) {
//                    cout << setprecision(2);
//                    cout << "ID: " << w << " Bounds : " << worker.E_min_local << " " << worker.E_max_local << endl;
//                    cout << worker.E_bins.transpose() << endl;
//                    cout.flush();
//                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//            }
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::microseconds(1000));
//        }
//        worker.dos_total.resize(1, 1);
//        worker.E_bins_total.resize(1);
//        worker.M_bins_total.resize(1);
//
//
//    }

