//
// Created by david on 2016-08-11.
//

#include "nmspc_WL_parallel_algorithms.h"
#define debug_swap      0
#define debug_merge     0
#define debug_bcast     0
#define debug_divide    0
using namespace std;

namespace mpi {
    void swap(class_worker &worker) {
        //Use MPI Tag in the 100-200 range
        if (timer::swap > constants::rate_swap) {
            timer::swap = 0;
            int abort;
            MPI_Allreduce(&worker.need_to_resize_global, &abort, 1, MPI_INT , MPI_MAX, MPI_COMM_WORLD);
            if (abort){return;}
            worker.t_swap.tic();
            counter::swaps++;
            int swap, copy;
            double dos_X, dos_Y;
            double E_X, E_Y, M_X, M_Y;
            int    E_X_idx, E_Y_idx, M_X_idx, M_Y_idx;
            int    in_window_up, in_window_dn;
            double E_min_up, E_max_up;
            double P_swap;      //Swap probability
            bool myTurn = math::mod(worker.world_ID, 2) == math::mod(counter::swaps, 2);

            int up = math::mod(worker.world_ID + 1, worker.world_size);
            int dn = math::mod(worker.world_ID - 1, worker.world_size);
            if (debug_swap){
                for (int w = 0; w < worker.world_size; w++){
                    if(w == worker.world_ID){
                        cout << "ID: "<< w << " Starting Swap. Myturn = " << myTurn << endl;
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
            worker.in_window = worker.check_in_window(worker.E);
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
            if (debug_swap){
                for (int w = 0; w < worker.world_size; w++){
                    if(w == worker.world_ID){
                        cout << "ID: "<< w<< " SendRecv Successful"
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
                //Make sure the last worker doesn't swap!
                go_ahead = go_ahead && worker.world_ID != worker.world_size - 1;
                MPI_Send(&go_ahead, 1, MPI_INT, up, 106, MPI_COMM_WORLD);
                MPI_Send(&copy, 1, MPI_INT, up, 107, MPI_COMM_WORLD);
            } else {
                MPI_Recv(&go_ahead, 1, MPI_INT, dn, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&copy, 1, MPI_INT, dn, 107, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if (debug_swap){
                for (int w = 0; w < worker.world_size; w++){
                    if(w == worker.world_ID){
                        cout << "ID: "<< w<< " Received goahead = " << go_ahead << endl;
                        cout << worker << endl;
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
            if (debug_swap){
                for (int w = 0; w < worker.world_size; w++){
                    if(w == worker.world_ID){
                        cout << "ID: "<< w<< "Starting swap. Swap = "  << swap << endl;
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
            //Now do the swapping if you got lucky
            if (myTurn) {
                if (go_ahead == 1) {
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
                if (go_ahead == 1) {
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
            worker.t_swap.toc();
            if (debug_swap){
                for (int w = 0; w < worker.world_size; w++){
                    if(w == worker.world_ID){
                        cout << "ID: "<< w<< "Swap OK "  << endl;
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        } else {
            timer::swap++;
        }
    }

    void merge(class_worker &worker) {
        if (worker.world_ID == 0 && debug_merge) {
            cout << "Merging" << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }

//        //set zero_values to nan
//        for (int w = 0; w < worker.world_size ; w++){
//            if(w == worker.world_ID){
//                cout << worker.dos << endl << endl;
//                cout << worker.dos << endl << endl << endl;
//
//                cout.flush();
//                std::this_thread::sleep_for(std::chrono::microseconds(1000));
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
        worker.dos = math::Zero_to_NaN(worker.dos);


        ArrayXXd dos_total, dos_temp, dos_recv;
        ArrayXd E_total, M_total, E_temp, M_temp, E_recv, M_recv;
        int E_sizes[worker.world_size];
        int M_sizes[worker.world_size];
        int E_size = (int) worker.E_bins.size();
        int M_size = (int) worker.M_bins.size();

        double diff;

        MPI_Gather(&E_size, 1, MPI_INT, E_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&M_size, 1, MPI_INT, M_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (worker.world_ID == 0 && debug_merge) {
            cout << "Gathering dos" << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        MPI_Barrier(MPI_COMM_WORLD);
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
                if (std::isnan(diff)){
                    cout << "Tried to concatenate "<< w-1 << " and " << w << ". Diff between two DOS is NaN, exiting" << endl;
                    cout << "Merge_idx    = [" <<  E_merge_idx    << ", " << M_merge_idx    << "]" << endl;
                    cout << "Merge_idx_up = [" <<  E_merge_idx_up << ", " << M_merge_idx_up << "]" << endl;

                    cout << "E_bins_total   " << E_total(E_merge_idx) << " [" << E_merge_idx << "]" << "    = "
                         << E_total.transpose() << endl;
                    cout << "E_bins_up      " << E_recv(E_merge_idx_up) << " [" << E_merge_idx_up << "]" << "    = "
                         << E_recv.transpose() << endl;
                    cout << "dos_tot" <<endl << dos_total << endl;
                    cout << "dos_rec" <<endl << dos_recv << endl;
                    MPI_Finalize();
                    exit(1);
                }
                math::add_to_nonzero_nonnan(dos_recv, diff);
                //Now now all doses should be the same height
                if (debug_merge) {
                    cout << "Concatenating " << w-1 << " and "<< w << endl;
                    cout << "Merge_idx    = [" <<  E_merge_idx    << ", " << M_merge_idx    << "]" << endl;
                    cout << "Merge_idx_up = [" <<  E_merge_idx_up << ", " << M_merge_idx_up << "]" << endl;

                    cout << "E_bins_total   " << E_total(E_merge_idx) << " [" << E_merge_idx << "]" << "    = "
                         << E_total.transpose() << endl;
                    cout << "E_bins_up      " << E_recv(E_merge_idx_up) << " [" << E_merge_idx_up << "]" << "    = "
                         << E_recv.transpose() << endl;
                    cout << "dos_tot" <<endl << dos_total << endl;
                    cout << "dos_rec" <<endl << dos_recv << endl;

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
            math::subtract_min_nonnan(dos_total);
            //Trim away nan rows for each worker
            math::remove_nan_rows(dos_total, E_total);

            worker.find_current_state();
            worker.dos_total = dos_total;
            worker.E_bins_total = E_total;
            worker.M_bins_total = M_total;
            if (debug_merge) {
                cout << "Finished Merging" << endl;
                cout << "E_bins_total = " << E_total.transpose() << endl<<endl;
                cout << "dos_total = " << endl << dos_total << endl;
                cout.flush();
                std::this_thread::sleep_for(std::chrono::microseconds(1000));
            }
        }
    }

    void broadcast(class_worker &worker) {
        if (worker.world_ID == 0 && debug_bcast) {
            cout << "Sending back to workers " << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        mpi::bcast_dynamic(worker.dos_total, MPI_DOUBLE, 0, worker.world_ID );
        mpi::bcast_dynamic(worker.E_bins_total, MPI_DOUBLE, 0, worker.world_ID );
        mpi::bcast_dynamic(worker.M_bins_total, MPI_DOUBLE, 0, worker.world_ID );

        counter::merges++;
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_bcast) {
            cout << "Merge OK " << endl;
            cout << "E_bins_total = " << worker.E_bins_total.transpose() << endl;
            cout << "Diff         = "
                 << (worker.E_bins_total.tail(worker.E_bins_total.size() - 1) - worker.E_bins_total.head(worker.E_bins_total.size() - 1)).transpose() << endl;
            cout << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10000));
        }
    }

    void divide_global_range_dos_volume(class_worker &worker) {
        //Update limits
        if (worker.world_ID == 0 && debug_divide) {
            cout << "Dividing. ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10000));
        }

        MPI_Barrier(MPI_COMM_WORLD);
        math::subtract_min_nonnan(worker.dos_total);
        if (worker.world_ID == 0){
            cout << "dos size " << worker.dos_total.rows()    << " x " << worker.dos_total.cols() << endl;
            cout << "E   size " << worker.E_bins_total.rows() << " x " << worker.E_bins_total.cols() << endl;
        }
        if (worker.world_ID == 0 && debug_divide) {
            cout << "Computing Volume ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10000));
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double global_volume = math::volume(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
        double local_volume = global_volume / worker.world_size;
        double x = constants::overlap_factor_dos_vol / (1 - constants::overlap_factor_dos_vol / 2)  ;
        //Add a little bit if there are too many workers (Add nothing if world_size == 2, and up to local_volume if world_size == inf)
        double overlap_range = local_volume * x ;//2.0*(worker.world_size - 2.0 + x)/worker.world_size;
//        cout << local_volume / global_volume << overlap_range << " " << local_volume*x << endl << endl;
        //Find the boundaries of the DOS domain that gives every worker the  same DOS volume to work on
        int E_min_local_idx, E_max_local_idx;
//        if (worker.world_ID == 0){
//            cout << "dos volume " << global_volume << endl;
//        }
        int min_width = 10;
        if (worker.world_ID == 0) {
            E_min_local_idx = 0;
            E_max_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                               (worker.world_ID + 1) * local_volume + overlap_range / 2);
            while (E_max_local_idx - E_min_local_idx < min_width) { E_max_local_idx++; }

        } else if (worker.world_ID == worker.world_size - 1) {
            E_min_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                               worker.world_ID * local_volume - overlap_range / 2);
            E_max_local_idx = (int) worker.E_bins_total.size() - 1;
            while (E_max_local_idx - E_min_local_idx < min_width) { E_min_local_idx--; }

        } else {
            E_min_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                               worker.world_ID * local_volume - overlap_range / 4);
            E_max_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                               (worker.world_ID + 1) * local_volume + overlap_range / 4);
            while (E_max_local_idx - E_min_local_idx < min_width) {
                E_min_local_idx--;
                E_max_local_idx++;
            }
        }
        E_min_local_idx = max(E_min_local_idx, 0);
        E_max_local_idx = min(E_max_local_idx, (int)worker.E_bins_total.size()-1);

        worker.E_min_local = worker.E_bins_total(E_min_local_idx);
        worker.E_max_local = worker.E_bins_total(E_max_local_idx);
        worker.M_min_local = worker.M_min_global;
        worker.M_max_local = worker.M_max_global;
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_divide) {
            cout << "...OK. Inherit from dos_total ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10000));
        }
        MPI_Barrier(MPI_COMM_WORLD);


        //Inherit the corresponding part of the dos
        int from = E_min_local_idx;
        int rows = E_max_local_idx - E_min_local_idx + 1;

        if (from + rows > worker.E_bins_total.size()) {
            cout << "TOO MANY ROWS   |  from + rows = " << from+rows << " E_bins_total.size() = " << worker.E_bins_total.size() << endl;
            cout << worker << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::seconds(1));
            MPI_Finalize();
            exit(15);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (E_max_local_idx <= E_min_local_idx) {
            cout << "Local range backwards! " << endl;
            cout << worker << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::seconds(1));
            MPI_Finalize();
            exit(16);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        worker.dos = worker.dos_total.middleRows(from, rows);
        worker.E_bins = worker.E_bins_total.segment(from, rows);
        worker.in_window = worker.check_in_window(worker.E);
        worker.histogram = ArrayXXi::Zero(worker.dos.rows(), worker.dos.cols());

        if (worker.in_window) {
            worker.E_idx = math::binary_search(worker.E_bins, worker.E);
            worker.M_idx = math::binary_search(worker.M_bins, worker.M);
        }
        if (worker.model.discrete_model) {
            worker.E_set.clear();
            for (int i = 0; i < worker.E_bins.size(); i++) {
                worker.E_set.insert(worker.E_bins(i));
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if (worker.world_ID == 0 && debug_divide) {
            cout << "...OK " << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        MPI_Barrier(MPI_COMM_WORLD);
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
        worker.dos_total.resize(1, 1);
        worker.E_bins_total.resize(1);
        worker.M_bins_total.resize(1);


    }

}