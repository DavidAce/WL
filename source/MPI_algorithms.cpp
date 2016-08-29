//
// Created by david on 2016-08-11.
//

#include "MPI_algorithms.h"

using namespace std;

namespace mpi {
    void swap(class_worker &worker) {
        //Use MPI Tag in the 100-200 range
        if (timer::swap > constants::rate_swap) {
            timer::swap = 0;
            worker.t_swap.tic();
            counter::swaps++;
            int swap, copy;
            double dos_X, dos_Y;
            double E_X, E_Y, M_X, M_Y;
            int E_X_idx, E_Y_idx, M_X_idx, M_Y_idx;
            double E_min_up, E_max_dn;
            double P_swap;      //Swap probability
            bool myTurn = math::mod(worker.world_ID, 2) == math::mod(counter::swaps, 2);

            int up = math::mod(worker.world_ID + 1, worker.world_size);
            int dn = math::mod(worker.world_ID - 1, worker.world_size);

            //Send current E and M to neighbors up and down. Receive X from below, Y from above.
            MPI_Sendrecv(&worker.E, 1, MPI_DOUBLE, up, 100, &E_X, 1, MPI_DOUBLE, dn, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.M, 1, MPI_DOUBLE, up, 101, &M_X, 1, MPI_DOUBLE, dn, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.E, 1, MPI_DOUBLE, dn, 102, &E_Y, 1, MPI_DOUBLE, up, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.M, 1, MPI_DOUBLE, dn, 103, &M_Y, 1, MPI_DOUBLE, up, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //Check if the neighbors position is within my overlap region. If so, find the indices.
            MPI_Sendrecv(&worker.E_min_local, 1, MPI_DOUBLE, dn, 104, &E_min_up, 1, MPI_DOUBLE, up, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.E_max_local, 1, MPI_DOUBLE, up, 105, &E_max_dn, 1, MPI_DOUBLE, dn, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //Now both swappees need to know if it is ok to go ahead with a swap.
            int go_ahead;
            if (myTurn) {
                go_ahead = worker.E >= E_min_up && E_Y <= worker.E_max_local;
                copy = !go_ahead && !worker.in_window && ((worker.E < E_Y && worker.E < worker.E_min_local) ||
                                                          (worker.E > E_Y && worker.E > worker.E_max_local));
                go_ahead = go_ahead && worker.world_ID != worker.world_size - 1;
                MPI_Send(&go_ahead, 1, MPI_INT, up, 106, MPI_COMM_WORLD);
                MPI_Send(&copy, 1, MPI_INT, up, 107, MPI_COMM_WORLD);
            } else {
                MPI_Recv(&go_ahead, 1, MPI_INT, dn, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&copy, 1, MPI_INT, dn, 107, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            //Now we should only be left with workers going ahead with a swap.
            if (myTurn) {
                if (go_ahead) {
                    MPI_Recv(&dos_X, 1, MPI_DOUBLE, up, 108, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&dos_Y, 1, MPI_DOUBLE, up, 109, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    E_Y_idx = math::binary_search(worker.E_bins.data(), E_Y, worker.E_bins.size());
                    M_Y_idx = math::binary_search(worker.M_bins.data(), M_Y, worker.M_bins.size());
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
                E_X_idx = math::binary_search(worker.E_bins.data(), E_X, worker.E_bins.size());
                M_X_idx = math::binary_search(worker.M_bins.data(), M_X, worker.M_bins.size());
                MPI_Send(&worker.dos(E_X_idx, M_X_idx), 1, MPI_DOUBLE, dn, 108, MPI_COMM_WORLD);
                MPI_Send(&worker.dos(worker.E_idx, worker.M_idx), 1, MPI_DOUBLE, dn, 109, MPI_COMM_WORLD);
                MPI_Recv(&swap, 1, MPI_INT, dn, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }


            //Now do the swapping if you got lucky
            if (myTurn) {
                if (go_ahead == 1) {
                    MPI_Sendrecv_replace(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, up, 111, up, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    worker.E = E_Y;
                    worker.M = M_Y;
                    worker.find_current_state();
                } else if (copy == 1) {
                    MPI_Recv(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, up, 112, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    worker.E = E_Y;
                    worker.M = M_Y;
                    worker.find_current_state();
                }
            }
            else {
                if (go_ahead == 1) {
                    MPI_Sendrecv_replace(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, dn, 111, dn, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    worker.E = E_X;
                    worker.M = M_X;
                    worker.find_current_state();
                } else if (copy == 1) {
                    MPI_Send(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, dn, 112, MPI_COMM_WORLD);
                }
            }
            counter::swap_accepts += swap;
            worker.t_swap.toc();

        } else {
            timer::swap++;
        }
    }

    void merge(class_worker &worker){
        if (worker.world_ID == 0 && debug_merge) {
            cout << "Merging" << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MatrixXd dos_total, dos_temp;
        VectorXd E_total, M_total, E_temp, M_temp;
        MatrixXd dos_recv;
        VectorXd E_recv;
        VectorXd M_recv;
        int E_sizes[worker.world_size];
        int M_sizes[worker.world_size];
        int E_size = (int)worker.E_bins.size();
        int M_size = (int)worker.M_bins.size();

        double diff;

        MPI_Gather(&E_size, 1, MPI_INT, E_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&M_size, 1, MPI_INT, M_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (worker.world_ID == 0 && debug_merge) {
            cout << "Receiving from " ;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0) {
            dos_total = worker.dos;
            E_total = worker.E_bins;
            M_total = worker.M_bins;
        }
        int from_total, from_up, rows_total,rows_up;
        int E_merge_idx, E_merge_idx_up;
        int M_merge_idx, M_merge_idx_up;
        for(int w = 1; w < worker.world_size; w++){
            if (worker.world_ID == 0){
                if (debug_merge) {
                    cout << w << " ";
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
                dos_recv.resize(E_sizes[w],M_sizes[w]);
                E_recv.resize(E_sizes[w]);
                M_recv.resize(M_sizes[w]);
                MPI_Recv(dos_recv.data(), E_sizes[w]*M_sizes[w], MPI_DOUBLE, w, 60, MPI_COMM_WORLD,MPI_STATUS_IGNORE );
                MPI_Recv(E_recv.data(), E_sizes[w], MPI_DOUBLE, w, 61, MPI_COMM_WORLD,MPI_STATUS_IGNORE );
                MPI_Recv(M_recv.data(), M_sizes[w], MPI_DOUBLE, w, 62, MPI_COMM_WORLD,MPI_STATUS_IGNORE );

                //Find coordinates on dos_total
                E_merge_idx        = math::find_matching_slope(dos_total, dos_recv, E_total, E_recv, M_total, M_recv);
                dos_total.row(E_merge_idx).maxCoeff(&M_merge_idx);
                //Find coordinates on received dos;
                E_merge_idx_up     = math::binary_search(E_recv.data(), E_total(E_merge_idx), E_sizes[w]);
                dos_recv.row(E_merge_idx_up).maxCoeff(&M_merge_idx_up);
                //Find difference between heights at these points
                diff = dos_total(E_merge_idx, M_merge_idx) - dos_recv(E_merge_idx_up, M_merge_idx_up) ;
                //Add that difference to the next
                math::add_to_nonzero(dos_recv, diff);
                //Now now all doses should be the same height
                if (debug_merge) {
                    cout << "Concatenating ID: "<< w << endl;
                    cout << "E_bins_up      " << E_recv(E_merge_idx_up) << " [" << E_merge_idx_up << "]" << "    = " << E_recv.transpose() << endl;
                    cout << "E_bins_total   " << E_total(E_merge_idx)      << " [" << E_merge_idx    << "]" << "    = " << E_total.transpose() << endl;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }

                //Compute the starting points and number of rows for concatenation
                from_total = 0;
                from_up    = E_merge_idx_up + 1;

                rows_total = E_merge_idx + 1;
                rows_up    = E_sizes[w] - E_merge_idx_up - 1;

                E_temp.resize(E_total.segment(from_total,rows_total).size() + E_recv.segment(from_up, rows_up).size());
                dos_temp.resize(dos_total.middleRows(from_total, rows_total).size() + dos_recv.middleRows(from_up,rows_up).size()  , M_size);

                dos_temp << dos_total.middleRows(from_total, rows_total),
                        dos_recv.middleRows(from_up,rows_up);
                E_temp   << E_total.segment(from_total,rows_total),
                        E_recv.segment(from_up, rows_up);
                dos_total = dos_temp;
                E_total   = E_temp  ;

                if (debug_merge) {
                    cout << "OK ";
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
            }else if(w == worker.world_ID){
                MPI_Send(worker.dos.data(),    E_size*M_size, MPI_DOUBLE, 0, 60, MPI_COMM_WORLD);
                MPI_Send(worker.E_bins.data(), E_size, MPI_DOUBLE, 0, 61, MPI_COMM_WORLD);
                MPI_Send(worker.M_bins.data(), M_size, MPI_DOUBLE, 0, 62, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (worker.world_ID == 0 && debug_merge) {
            cout << "Sending back to workers " << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }
        if(worker.world_ID == 0){
            worker.dos_total    = dos_total;
            worker.E_bins_total = E_total;
            worker.M_bins_total = M_total;
            E_size              = (int)E_total.size();
            M_size              = (int)M_total.size();
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&E_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&M_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        worker.dos_total.conservativeResize(E_size,M_size);
        worker.E_bins_total.conservativeResize(E_size);
        worker.M_bins_total.conservativeResize(M_size);
        MPI_Bcast(worker.dos_total.data()   , E_size*M_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(worker.E_bins_total.data(), E_size       , MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(worker.M_bins_total.data(), M_size       , MPI_DOUBLE, 0, MPI_COMM_WORLD);
        counter::merges++;
        MPI_Barrier(MPI_COMM_WORLD);

        if (worker.world_ID == 0 && debug_merge) {
            cout << "Merge OK " << endl;
            cout << "E_bins_total = " << worker.E_bins_total.transpose() << endl;
            cout << "Diff         = " << (worker.E_bins_total.tail(E_size - 1) - worker.E_bins_total.head(E_size-1)).transpose() << endl;
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
        math::subtract_min_nonzero(worker.dos_total);
        if (worker.world_ID == 0 && debug_divide) {
            cout << "Computing Volume ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10000));
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double global_volume = math::volume(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
        double local_volume  = global_volume / worker.world_size;
        double x = local_volume * constants::overlap_factor / (1 - constants::overlap_factor / 2);
        cout << x <<" " << local_volume << endl;
        //Find the boundaries of the DOS domain that gives every worker the  same DOS volume to work on
        int E_min_local_idx, E_max_local_idx;

        if (worker.world_ID == 0) {
            E_min_local_idx  = 0;
            E_max_local_idx     = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total, (worker.world_ID + 1)*local_volume + x/2);
            while(E_max_local_idx - E_min_local_idx < 3){E_max_local_idx++;}

        }else if (worker.world_ID == worker.world_size - 1){
            E_min_local_idx     = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total, worker.world_ID      *local_volume - x/2);
            E_max_local_idx     = (int) worker.E_bins_total.size()-1;
            while(E_max_local_idx - E_min_local_idx < 3){E_min_local_idx--;}

        }else{
            E_min_local_idx     = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total, worker.world_ID      *local_volume - x/4);
            E_max_local_idx     = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total, (worker.world_ID + 1)*local_volume + x/4);
            while(E_max_local_idx - E_min_local_idx < 3){E_min_local_idx--; E_max_local_idx++;}
        }

        worker.E_min_local = worker.E_bins_total(E_min_local_idx) ;
        worker.E_max_local = worker.E_bins_total(E_max_local_idx) ;
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

        if (from+rows > worker.E_bins_total.size()){
            cout << "TOO MANY ROWS "<< endl;
            cout << worker << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::seconds(1));
            exit(15);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (E_max_local_idx <= E_min_local_idx){
            cout << "Local range backwards! "<< endl;
            cout << worker << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::seconds(1));
            exit(16);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        worker.dos         = worker.dos_total.middleRows(from, rows);
        worker.E_bins      = worker.E_bins_total.segment(from, rows);
        worker.in_window   = worker.check_in_window(worker.E);
        if (worker.in_window){
            worker.E_idx = math::binary_search(worker.E_bins.data(), worker.E, worker.E_bins.size());
            worker.M_idx = math::binary_search(worker.M_bins.data(), worker.M, worker.M_bins.size());
        }
        worker.histogram.conservativeResize(worker.dos.rows(), worker.dos.cols());
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
            for (int w = 0; w < worker.world_size; w++){
                if (w == worker.world_ID){
                    cout << setprecision(2);
                    cout << "ID: "<< w << " Bounds : " << worker.E_min_local<< " " << worker.E_max_local <<  endl;
                    cout << worker.E_bins.transpose() << endl;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        worker.dos_total.resize(1,1);
        worker.E_bins_total.resize(1);
        worker.M_bins_total.resize(1);


    }

    void divide_global_range_dos_volume2(class_worker &worker) {
        //Update limits
        if (worker.world_ID == 0 && debug_divide) {
            cout << "Dividing. ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10000));
        }

        MPI_Barrier(MPI_COMM_WORLD);
        math::subtract_min_nonzero(worker.dos_total);
        if (worker.world_ID == 0 && debug_divide) {
            cout << "Computing Volume ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10000));
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double global_volume    = math::volume(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
        double local_volume     = global_volume / worker.world_size;
        int E_idx_low           = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total, local_volume * (worker.world_ID));
        int E_idx_high          = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total, local_volume * (worker.world_ID + 1));

        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_divide) {
            cout << "...OK. Finding Ranges ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10000));
        }
        MPI_Barrier(MPI_COMM_WORLD);

        double local_range = worker.E_bins_total(E_idx_high) - worker.E_bins_total(E_idx_low);
        double overlap_range_up = constants::overlap_factor*local_range / (1 - 2.0*constants::overlap_factor) / 2;
        double overlap_range_dn;
        int dn = worker.world_ID > 0 ? worker.world_ID - 1 : worker.world_size - 1;
        int up = worker.world_ID < worker.world_size - 1 ? worker.world_ID + 1 : 0;
        MPI_Sendrecv(&overlap_range_up, 1, MPI_DOUBLE, up, 350, &overlap_range_dn, 1, MPI_DOUBLE, dn, 350, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        int E_min_local_idx, E_max_local_idx;
        if (worker.world_ID == 0) {
            E_min_local_idx = 0;
            E_max_local_idx = math::binary_search(worker.E_bins_total.data(),  worker.E_bins_total(E_idx_high) + overlap_range_up, worker.E_bins_total.size());
            while(E_max_local_idx - E_min_local_idx < 3){E_max_local_idx++;}
            worker.E_min_local = worker.E_bins_total(E_min_local_idx);
            worker.E_max_local = worker.E_bins_total(E_max_local_idx);
            if (worker.E_min_local > worker.E_min_global){
                cout << endl << "Lowest local E is larger than global E!" << endl;
                cout << "E_min_local = "  << worker.E_min_local << " E_min_global = " << worker.E_min_global << endl;
                cout.flush();
                std::this_thread::sleep_for(std::chrono::microseconds(10000));
//                exit(40);
            }


        } else if (worker.world_ID == worker.world_size - 1) {
            E_min_local_idx = math::binary_search(worker.E_bins_total.data(), worker.E_bins_total(E_idx_low) - overlap_range_dn, worker.E_bins_total.size());
            E_max_local_idx = (int)worker.E_bins_total.size() - 1;
            while(E_max_local_idx - E_min_local_idx < 3){E_min_local_idx--;}
            worker.E_min_local = worker.E_bins_total(E_min_local_idx);
            worker.E_max_local = worker.E_bins_total(E_max_local_idx);
            if (worker.E_max_local < worker.E_max_global){
                cout << endl << "Highest local E is lower than global E!" << endl;
                cout << "E_max_local = "  << worker.E_max_local << " E_max_global = " << worker.E_max_global << endl;
                cout.flush();
                std::this_thread::sleep_for(std::chrono::microseconds(10000));
//                exit(41);
            }

        } else {
            E_min_local_idx = math::binary_search(worker.E_bins_total.data(), worker.E_bins_total(E_idx_low) - overlap_range_dn, worker.E_bins_total.size());
            E_max_local_idx = math::binary_search(worker.E_bins_total.data(), worker.E_bins_total(E_idx_high) + overlap_range_up, worker.E_bins_total.size());
            while(E_max_local_idx - E_min_local_idx < 3){E_min_local_idx--; E_max_local_idx++;}
            worker.E_min_local = worker.E_bins_total(E_min_local_idx) ;
            worker.E_max_local = worker.E_bins_total(E_max_local_idx) ;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0 && debug_divide) {
            cout << "...OK. Inherit from dos_total ";
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(10000));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        worker.M_min_local = worker.M_min_global;
        worker.M_max_local = worker.M_max_global;

        //Inherit the corresponding part of the dos
        int from = E_min_local_idx;
        int rows = E_max_local_idx - E_min_local_idx + 1;

        if (from+rows > worker.E_bins_total.size()){
            cout << "TOO MANY ROWS "<< endl;
            cout << worker << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::seconds(1));
            exit(15);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (E_max_local_idx <= E_min_local_idx){
            cout << "Local range backwards! "<< endl;
            cout << worker << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::seconds(1));
            exit(16);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        worker.dos         = worker.dos_total.middleRows(from, rows);
        worker.E_bins      = worker.E_bins_total.segment(from, rows);
        worker.in_window   = worker.check_in_window(worker.E);
        if (worker.in_window){
            worker.E_idx = math::binary_search(worker.E_bins.data(), worker.E, worker.E_bins.size());
            worker.M_idx = math::binary_search(worker.M_bins.data(), worker.M, worker.M_bins.size());
        }
        worker.histogram.conservativeResize(worker.dos.rows(), worker.dos.cols());
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
            for (int w = 0; w < worker.world_size; w++){
                if (w == worker.world_ID){
                    cout << setprecision(2);
                    cout << "ID: "<< w << " Bounds : " << worker.E_min_local<< " " << worker.E_max_local <<  endl;
                    cout << worker.E_bins.transpose() << endl;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        worker.dos_total.resize(1,1);
        worker.E_bins_total.resize(1);
        worker.M_bins_total.resize(1);
        

    }

}


//void merge3(class_worker &worker) {
//    //Use MPI-tag in the 200-300 range
//    //Get id of neighbours up and down.
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "Merging " << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10));
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    int dn = worker.world_ID > 0 ? worker.world_ID - 1 : worker.world_size - 1;
//    int up = worker.world_ID < worker.world_size - 1 ? worker.world_ID + 1 : 0;
//    int E_size_up, E_size = (int) worker.dos.rows();
//    int M_size_up, M_size = (int) worker.dos.cols();
//
//    //Find out the size of the data of the neighbor above
//    MPI_Sendrecv(&E_size, 1, MPI_INT, dn, 200, &E_size_up, 1, MPI_INT, up, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    MPI_Sendrecv(&M_size, 1, MPI_INT, dn, 201, &M_size_up, 1, MPI_INT, up, 201, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//    MatrixXd dos_up(E_size_up, M_size_up);
//    VectorXd E_bins_up(E_size_up), M_bins_up(M_size_up);
//
//    //Receive it as well
//    MPI_Sendrecv(worker.dos.data(), E_size * M_size, MPI_DOUBLE, dn, 202, dos_up.data(), E_size_up * M_size_up,
//                 MPI_DOUBLE, up, 202, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    MPI_Sendrecv(worker.E_bins.data(), E_size, MPI_DOUBLE, dn, 203, E_bins_up.data(), E_size_up, MPI_DOUBLE, up,
//                 203, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    MPI_Sendrecv(worker.M_bins.data(), M_size, MPI_DOUBLE, dn, 204, M_bins_up.data(), M_size_up, MPI_DOUBLE, up,
//                 204, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    int error = 0;
//    if (worker.E_bins.maxCoeff() < E_bins_up.minCoeff() && up != 0 && debug_merge) {
//        cout << endl << "ID: " << worker.world_ID << " Error in merge_windows(). Windows do not overlap." << endl;
//        cout << "Max E = " << worker.E_bins.maxCoeff() << " Min E above = " << E_bins_up.maxCoeff() << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10));
//        error = 1;
//    }
//    if (worker.E_bins.maxCoeff() > E_bins_up.maxCoeff() && up != 0 && debug_merge) {
//        cout << "ID: " << worker.world_ID << " Error in merge_windows(). Overlap too large." << endl;
//        cout << worker << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10));
//        error = 1;
//    }
//    MPI_Allreduce(&error, &error, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//    if (error) {
//        for (int w = 0; w < worker.world_size - 1; w++) {
//            if (worker.world_ID == w) {
//                cout << "E_bins_up: " << E_bins_up.transpose() << endl;
//                cout << worker << endl;
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//        exit(4);
//    }
//    Vector3d u1, u2; //Vectors connecting adjacent 3 orthogonal points on DOS
//    Vector3d v_up, v_rt, v_dn, v_lf; //Vectors connecting adjacent 3 orthogonal points on DOS
//    VectorXd sum(E_size);
//    sum.fill(0);
//    int col = 1, num;
//    int E_merge_idx;    //Index of merging point
//    int M_merge_idx;    //Index of merging point
//    double E_merge;
//
//    int x, y; //Coordinates closest to i, j on dos of neighbor above.
//    for (int i = 1; i < E_size - 1; i++) {
//        num = 0;
//        for (int j = 1; j < M_size - 1; j++) {
//            if (worker.E_bins(i) <= E_bins_up.minCoeff()) { continue; }
//            if (worker.E_bins(i) >= E_bins_up.maxCoeff()) { continue; }
//            if (worker.M_bins(j) >= M_bins_up.maxCoeff()) { continue; }
//            if (worker.M_bins(j) <= M_bins_up.minCoeff()) { continue; }
//            x = math::upper_bound(E_bins_up.data(), worker.E_bins(i), E_size_up);
//            y = math::upper_bound(M_bins_up.data(), worker.M_bins(i), M_size_up);
//            if (x <= 0 || x >= E_size_up - 1) { continue; }
//            if (y <= 0 || y >= M_size_up - 1) { continue; }
//
//            v_up << worker.E_bins(i - 1) - worker.E_bins(i), 0, worker.dos(i - 1, j) - worker.dos(i, j);
//            v_dn << worker.E_bins(i + 1) - worker.E_bins(i), 0, worker.dos(i + 1, j) - worker.dos(i, j);
//            v_rt << 0, worker.M_bins(j + 1) - worker.M_bins(j), worker.dos(i, j + 1) - worker.dos(i, j);
//            v_lf << 0, worker.M_bins(j - 1) - worker.M_bins(j), worker.dos(i, j - 1) - worker.dos(i, j);
//            u1 = (v_up.cross(v_rt) + v_rt.cross(v_dn) + v_dn.cross(v_lf) + v_lf.cross(v_dn)).normalized();
//
//            v_up << E_bins_up(x - 1) - E_bins_up(x), 0, dos_up(x - 1, y) - dos_up(x, y);
//            v_dn << E_bins_up(x + 1) - E_bins_up(x), 0, dos_up(x + 1, y) - dos_up(x, y);
//            v_rt << 0, M_bins_up(y + 1) - M_bins_up(y), dos_up(x, y + 1) - dos_up(x, y);
//            v_lf << 0, M_bins_up(y - 1) - M_bins_up(y), dos_up(x, y - 1) - dos_up(x, y);
//
//            u2 = (v_up.cross(v_rt) + v_rt.cross(v_dn) + v_dn.cross(v_lf) + v_lf.cross(v_dn)).normalized();
//            sum(col) += u1.dot(u2);
//            num++;
//        }
//        sum(col) /= max(1, num);
//        col++;
//    }
//    sum.maxCoeff(&E_merge_idx);
//    E_merge_idx = sum.sum() == 0 ? E_size - 1 : E_merge_idx;
//    if (worker.world_ID == worker.world_size-1){
//        E_merge_idx = E_size_up - 1;
//    }
//    worker.dos.row(E_merge_idx).maxCoeff(&M_merge_idx);
//    E_merge     = worker.E_bins(E_merge_idx);
//
//    if (debug_merge) {
//        cout << "ID: " << worker.world_ID
//             << " E_size = " << E_size
//             << " E_size_up = " << E_size_up
//             << " M_size = " << M_size
//             << " M_size_up " << M_size_up
//             << " E_merge_idx "<< E_merge_idx
//             << " E_merge "     << E_merge
//             << endl;
//    }
//
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "Found merging points... " << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    //Now send all to world_ID == 0, to do the merging
//
//    worker.dos_total         = worker.dos;
//    worker.E_bins_total      = worker.E_bins;
//    worker.M_bins_total      = worker.M_bins;
//    E_size                  = (int) worker.dos.rows();
//    M_size                  = (int) worker.dos.cols();
//    int E_merge_idx_recv, M_merge_idx_recv;
//    int i,j;
//    int E_size_total, E_idx_start;
//    double diff;
//    MatrixXd dos_recv;
//    VectorXd E_bins_recv;
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "Received from ";
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int w = 1; w < worker.world_size; w++) {
//        if (worker.world_ID == 0){
//            MPI_Recv(&E_size_up, 1, MPI_INT, w, 207, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            MPI_Recv(&M_size_up, 1, MPI_INT, w, 208, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//            dos_recv.resize(E_size_up, M_size_up);
//            E_bins_recv.resize(E_size_up);
//            MPI_Recv(dos_recv.data()   , (int)dos_recv.size()    , MPI_DOUBLE, w, 209, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//            MPI_Recv(E_bins_recv.data(), (int)E_bins_recv.size() , MPI_DOUBLE, w, 210, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//
//            if (worker.E_bins_total.maxCoeff() <= E_bins_recv.minCoeff()){
//                cout << "ID: " << worker.world_ID << " NO OVERLAP!!" << endl;
//                cout << worker.E_bins_total.transpose() << endl;
//                exit(31);
//            }
////                E_size_total = (int) worker.E_bins_total.size() + E_size_up;
//            E_merge_idx_recv = math::binary_search(E_bins_recv.data(), E_merge, E_size_up);
//            E_merge_idx      = math::binary_search(worker.E_bins_total.data(), E_merge, worker.E_bins_total.size());
//            worker.dos_total.row(E_merge_idx).maxCoeff(&M_merge_idx);
//            dos_recv.row(E_merge_idx_recv).maxCoeff(&M_merge_idx_recv);
//            diff = worker.dos_total(E_merge_idx, M_merge_idx) - dos_recv(E_merge_idx_recv, M_merge_idx_recv);
//            math::add_to_nonzero(dos_recv, diff);
//            i = E_merge_idx; j = E_merge_idx_recv;
//            while (j < E_size_up){
//
//                if (worker.E_bins_total.size() < i + 1){
//                    worker.dos_total.conservativeResize(i+1, Eigen::NoChange);
//                    worker.E_bins_total.conservativeResize(i+1);
//                }
//
//                worker.E_bins_total(i) = E_bins_recv(j);
//                worker.dos_total.row(i) = dos_recv.row(i);
//                i++;
//                j++;
//            }
//
//
//
//            MPI_Recv(&E_merge, 1 , MPI_DOUBLE, w, 211, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//            if (debug_merge) {
//                cout << w << " ";
//                cout.flush();
//                std::this_thread::sleep_for(std::chrono::microseconds(10000));
//            }
//
//        }else if(worker.world_ID == w){
//            MPI_Send(&E_size     , 1, MPI_INT, 0, 207, MPI_COMM_WORLD);
//            MPI_Send(&M_size     , 1, MPI_INT, 0, 208, MPI_COMM_WORLD);
//            MPI_Send(worker.dos.data(), E_size*M_size, MPI_DOUBLE, 0, 209, MPI_COMM_WORLD);
//            MPI_Send(worker.E_bins.data()  , E_size  , MPI_DOUBLE, 0, 210, MPI_COMM_WORLD);
//            MPI_Send(&E_merge, 1  , MPI_DOUBLE, 0, 211, MPI_COMM_WORLD);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "OK!" << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(worker.world_ID == 0){
//        E_size = (int)worker.dos_total.rows();
//        M_size = (int)worker.dos_total.cols();
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Bcast(&E_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&M_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    worker.dos_total.resize(E_size,M_size);
//    worker.E_bins_total.resize(E_size);
//    worker.M_bins_total.resize(M_size);
//    MPI_Bcast(worker.dos_total.data()   , E_size*M_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    MPI_Bcast(worker.E_bins_total.data(), E_size       , MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    MPI_Bcast(worker.M_bins_total.data(), M_size       , MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    counter::merges++;
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "Merge OK " << endl;
//        cout << "E_bins_total = " << worker.E_bins_total.transpose() << endl;
//        cout << "Diff         = " << (worker.E_bins_total.tail(E_size - 1) - worker.E_bins_total.head(E_size-1)).transpose() << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//}
//
//void merge2 (class_worker &worker) {
//    //Use MPI-tag in the 200-300 range
//    //Get id of neighbours up and down.
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "Merging " << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10));
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    int dn = worker.world_ID > 0 ? worker.world_ID - 1 : worker.world_size - 1;
//    int up = worker.world_ID < worker.world_size - 1 ? worker.world_ID + 1 : 0;
//    int E_size_up, E_size = (int) worker.dos.rows();
//    int M_size_up, M_size = (int) worker.dos.cols();
//
//    //Find out the size of the data of the neighbor above
//    MPI_Sendrecv(&E_size, 1, MPI_INT, dn, 200, &E_size_up, 1, MPI_INT, up, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    MPI_Sendrecv(&M_size, 1, MPI_INT, dn, 201, &M_size_up, 1, MPI_INT, up, 201, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    if (debug_merge) {
//        cout << "ID: " << worker.world_ID
//             << " E_size = " << E_size
//             << " E_size_up = " << E_size_up
//             << " M_size = " << M_size
//             << " M_size_up " << M_size_up << endl;
//    }
//
////    MPI_Barrier(MPI_COMM_WORLD);
//
//    MatrixXd dos_up(E_size_up, M_size_up);
//    VectorXd E_bins_up(E_size_up);
//    VectorXd M_bins_up(M_size_up);
//    //Receive it as well
//    MPI_Sendrecv(worker.dos.data(), E_size * M_size, MPI_DOUBLE, dn, 202, dos_up.data(), E_size_up * M_size_up,
//                 MPI_DOUBLE, up, 202, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    MPI_Sendrecv(worker.E_bins.data(), E_size, MPI_DOUBLE, dn, 203, E_bins_up.data(), E_size_up, MPI_DOUBLE, up, 203,
//                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    MPI_Sendrecv(worker.M_bins.data(), M_size, MPI_DOUBLE, dn, 204, M_bins_up.data(), M_size_up, MPI_DOUBLE, up, 204,
//                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    int error = 0;
//    if (worker.E_bins.maxCoeff() < E_bins_up.minCoeff() && up != 0 && debug_merge) {
//        cout << endl << "ID: " << worker.world_ID << " Error in merge_windows(). Windows do not overlap." << endl;
//        cout << "Max E = " <<worker.E_bins.maxCoeff() << " Min E above = " << E_bins_up.maxCoeff() << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10));
//        error = 1;
//    }
//    if (worker.E_bins.maxCoeff() > E_bins_up.maxCoeff() && up != 0 && debug_merge) {
//        cout << "ID: " << worker.world_ID << " Error in merge_windows(). Overlap too large." << endl;
//        cout << worker << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10));
//        error = 1;
//    }
//    MPI_Allreduce(&error, &error, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//    if (error) {
//        for (int w = 0; w < worker.world_size - 1; w++) {
//            if (worker.world_ID == w) {
//                cout << "E_bins_up: " << E_bins_up.transpose() << endl;
//                cout << worker << endl;
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//        exit(4);
//    }
//    Vector3d u1, u2; //Vectors connecting adjacent 3 orthogonal points on DOS
//    Vector3d v_up, v_rt, v_dn, v_lf; //Vectors connecting adjacent 3 orthogonal points on DOS
//    VectorXd sum(E_size);
//    sum.fill(0);
//    int col = 1, num;
//    int E_merge_idx;    //Index of merging point
//    int M_merge_idx;    //Index of merging point
//    int E_merge_idx_up; //Index of merging point above
//    int M_merge_idx_up; //Index of merging point above
//    int E_merge_idx_dn = 0;
//
//    int x, y; //Coordinates closest to i, j on dos of neighbor above.
//    for (int w = 0; w < worker.world_size; w++) {
//        if (worker.world_ID == w) {
//            for (int i = 1; i < E_size - 1; i++) {
//                num = 0;
//                for (int j = 1; j < M_size - 1; j++) {
//                    if (worker.E_bins(i) <= E_bins_up.minCoeff()) { continue; }
//                    if (worker.E_bins(i) >= E_bins_up.maxCoeff()) { continue; }
//                    if (worker.M_bins(j) >= M_bins_up.maxCoeff()) { continue; }
//                    if (worker.M_bins(j) <= M_bins_up.minCoeff()) { continue; }
//                    x = math::upper_bound(E_bins_up.data(), worker.E_bins(i), E_size_up);
//                    y = math::upper_bound(M_bins_up.data(), worker.M_bins(i), M_size_up);
//                    if (x <= 0 || x >= E_size_up - 1) { continue; }
//                    if (y <= 0 || y >= M_size_up - 1) { continue; }
//
//                    v_up << worker.E_bins(i - 1) - worker.E_bins(i) , 0                                      , worker.dos(i - 1, j) - worker.dos(i, j);
//                    v_dn << worker.E_bins(i + 1) - worker.E_bins(i) , 0                                      , worker.dos(i + 1, j) - worker.dos(i, j);
//                    v_rt << 0                                       , worker.M_bins(j + 1) - worker.M_bins(j), worker.dos(i, j + 1) - worker.dos(i, j);
//                    v_lf << 0                                       , worker.M_bins(j - 1) - worker.M_bins(j), worker.dos(i, j - 1) - worker.dos(i, j);
//                    u1 = (v_up.cross(v_rt) + v_rt.cross(v_dn) + v_dn.cross(v_lf) + v_lf.cross(v_dn)).normalized();
//
//                    v_up << E_bins_up(x - 1) - E_bins_up(x) , 0                              , dos_up(x - 1, y) - dos_up(x, y);
//                    v_dn << E_bins_up(x + 1) - E_bins_up(x) , 0                              , dos_up(x + 1, y) - dos_up(x, y);
//                    v_rt << 0                               , M_bins_up(y + 1) - M_bins_up(y), dos_up(x, y + 1) - dos_up(x, y);
//                    v_lf << 0                               , M_bins_up(y - 1) - M_bins_up(y), dos_up(x, y - 1) - dos_up(x, y);
//
//                    u2 = (v_up.cross(v_rt) + v_rt.cross(v_dn) + v_dn.cross(v_lf) + v_lf.cross(v_dn)).normalized();
//                    sum(col) += u1.dot(u2);
//                    num++;
//                }
//                sum(col) /= max(1, num);
//                col++;
//            }
//            sum.maxCoeff(&E_merge_idx);
//            E_merge_idx = sum.sum() == 0 ? E_size - 1 : max(E_merge_idx_dn, E_merge_idx);
//            worker.dos.row(E_merge_idx).maxCoeff(&M_merge_idx);
//            E_merge_idx_up = math::binary_search(E_bins_up.data(), worker.E_bins(E_merge_idx), E_size_up);
//            M_merge_idx_up = math::binary_search(M_bins_up.data(), worker.M_bins(M_merge_idx), M_size_up);
////                if (E_bins_up(E_merge_idx_up) > worker.E_bins(E_merge_idx)){
////                    E_merge_idx_up--;
////                }
//            if (E_merge_idx_up < 0){
//                cout << "ID: " << worker.world_ID << " E_merge_idx_up not found: " << E_merge_idx_up << endl;
//                cout.flush();
//                std::this_thread::sleep_for(std::chrono::microseconds(10000));
//                exit(20);
//            }
//            MPI_Send(&E_merge_idx_up, 1, MPI_INT, up, 266, MPI_COMM_WORLD);
//        }else if (worker.world_ID == w + 1){
//            MPI_Recv(&E_merge_idx_dn, 1, MPI_INT, dn, 266,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//    if (worker.world_ID == worker.world_size-1){
//        E_merge_idx = E_size_up - 1;
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (E_merge_idx_dn <= 0 && debug_merge) {
//        cout << "ID: " << worker.world_ID << " E_merge_idx_dn probably too small: " << E_merge_idx_dn << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    if (E_merge_idx_dn >= E_merge_idx && worker.world_ID != 0 && worker.world_ID < worker.world_size -1 && debug_merge) {
//        cout << "ID: " << worker.world_ID << " E_merge_idx_dn larger than E_merge_idx: " << E_merge_idx_dn << " < " << E_merge_idx << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "Found merging points... " << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    double diff;
//
//    for (int w = 0; w < worker.world_size - 1; w++) {
//        if (worker.world_ID == w) {
//            //Adjust the dos above. The index has to correspond to the same energy level!!
//            diff = worker.dos(E_merge_idx, M_merge_idx) - dos_up(E_merge_idx_up, M_merge_idx_up);
//
//            for (int j = 0; j < M_size_up; j++) {
//                for (int i = 0; i < E_size_up; i++) {
//                    if (dos_up(i, j) == 0) { continue; }
//                    else {
//                        dos_up(i, j) += diff;
//                        dos_up(i, j) = dos_up(i, j) <= 0 ? 0 : dos_up(i, j);
//                    }
//                }
//            }
//            MPI_Send(dos_up.data(),E_size_up*M_size_up, MPI_DOUBLE, up,206,MPI_COMM_WORLD);
//        }
//        else if ( worker.world_ID == w+1){
//            MPI_Recv(worker.dos.data(),E_size*M_size, MPI_DOUBLE, dn,206,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "Subtracted diff... " << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    //Now all workers have their dos aligned with the others. Now just send all to master?
//    //Start by trimming off excess
//
//
//    int from = E_merge_idx_dn;
//    int rows = E_merge_idx - E_merge_idx_dn + 1;
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (from + rows > E_size  && debug_merge){
//        cout << "ID: " << worker.world_ID<< " Rows out of bounds! "<< "from+rows = " << from+rows << " E_size = " << E_size <<  endl;
//        cout << "E_merge_idx = " << E_merge_idx << " E_merge_idx_dn" << E_merge_idx_dn << endl;
//        cout << "from = " << from << " rows = " << rows << endl;
//
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        exit(20);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    if (rows == 0  && debug_merge){
//        cout << "ID: " << worker.world_ID << "Warning! No rows! " << rows << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        rows = min(1,rows);
//
////        exit(20);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    if ((from < 0 || from > E_size - 1)  && debug_merge){
//        cout << "ID: " << worker.world_ID << "from_idx out of bounds! :" << from << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        exit(20);
//    }
//
//    worker.dos_total         = worker.dos.middleRows(from, rows);
//    worker.E_bins_total      = worker.E_bins.segment(from, rows);
//    worker.M_bins_total      = worker.M_bins;
//    E_size                  = (int) worker.dos_total.rows();
//    M_size                  = (int) worker.dos_total.cols();
//    MPI_Barrier(MPI_COMM_WORLD);
//
////        for (int w = 0 ; w < worker.world_size; w++){
////            if (worker.world_ID == w){
////                cout << endl << worker.dos_total << endl;
////                cout.flush();
////                std::this_thread::sleep_for(std::chrono::microseconds(10000));
////            }
////            MPI_Barrier(MPI_COMM_WORLD);
////        }
//
//
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "Initialized dos_total... " << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    //Now send all to world_ID == 0, to do the merging
//    MatrixXd dos_recv;
//    VectorXd E_bins_recv;
//    int E_new_size;
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "Received from ";
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    for (int w = 1; w < worker.world_size; w++) {
//        if (worker.world_ID == 0){
//            MPI_Recv(&E_size_up, 1, MPI_INT, w, 207, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            MPI_Recv(&M_size_up, 1, MPI_INT, w, 208, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            dos_recv.resize(E_size_up, M_size_up);
//            E_bins_recv.resize(E_size_up);
//            MPI_Recv(dos_recv.data()   , (int)dos_recv.size()    , MPI_DOUBLE, w, 209, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//            MPI_Recv(E_bins_recv.data(), (int)E_bins_recv.size() , MPI_DOUBLE, w, 210, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//            E_new_size = (int)worker.dos_total.rows() + (int)dos_recv.rows();
//            worker.dos_total.conservativeResize(E_new_size, M_size_up);
//            worker.E_bins_total.conservativeResize(E_new_size , 1);
//            worker.dos_total.bottomRows(E_size_up)        = dos_recv;
//            worker.E_bins_total.tail(E_size_up)           = E_bins_recv;
//            if (debug_merge) {
//                cout << w << " ";
//                cout.flush();
//                std::this_thread::sleep_for(std::chrono::microseconds(10000));
//            }
//
//        }else if(worker.world_ID == w){
//            MPI_Send(&E_size     , 1, MPI_INT, 0, 207, MPI_COMM_WORLD);
//            MPI_Send(&M_size     , 1, MPI_INT, 0, 208, MPI_COMM_WORLD);
//            MPI_Send(worker.dos_total.data(), E_size*M_size, MPI_DOUBLE, 0, 209, MPI_COMM_WORLD);
//            MPI_Send(worker.E_bins_total.data()  , E_size  , MPI_DOUBLE, 0, 210, MPI_COMM_WORLD);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "OK!" << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(worker.world_ID == 0){
//        E_size = (int)worker.dos_total.rows();
//        M_size = (int)worker.dos_total.cols();
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Bcast(&E_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&M_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    worker.dos_total.resize(E_size,M_size);
//    worker.E_bins_total.resize(E_size);
//    worker.M_bins_total.resize(M_size);
//    MPI_Bcast(worker.dos_total.data()   , E_size*M_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    MPI_Bcast(worker.E_bins_total.data(), E_size       , MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    MPI_Bcast(worker.M_bins_total.data(), M_size       , MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    counter::merges++;
//    if (worker.world_ID == 0 && debug_merge) {
//        cout << "Merge OK " << endl;
//        cout << "E_bins_total = " << worker.E_bins_total.transpose() << endl;
//        cout << "Diff         = " << (worker.E_bins_total.tail(E_size - 1) - worker.E_bins_total.head(E_size-1)).transpose() << endl;
//        cout.flush();
//        std::this_thread::sleep_for(std::chrono::microseconds(10000));
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//}

//    void divide_global_range_dos_volume(class_worker &worker) {
//        //Update limits
//        if (worker.world_ID == 0 && debug_divide) {
//            cout << "Dividing. ";
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        }
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        math::subtract_min_nonzero(worker.dos_total);
//        if (worker.world_ID == 0 && debug_divide) {
//            cout << "Computing Volume ";
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        }
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        double global_volume    = math::volume(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
//        double local_volume     = global_volume / worker.world_size;
//        int E_idx_low           = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total, local_volume * (worker.world_ID));
//        int E_idx_high          = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total, local_volume * (worker.world_ID + 1));
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (worker.world_ID == 0 && debug_divide) {
//            cout << "...OK. Finding Ranges ";
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        double local_range = worker.E_bins_total(E_idx_high) - worker.E_bins_total(E_idx_low);
//
////        double overlap_range = local_range * constants::overlap_factor;
//        double overlap_range = 1.0/(1.0/constants::overlap_factor + 0.5) * local_range / 2;
//        int E_min_local_idx, E_max_local_idx;
//        if (worker.world_ID == 0) {
//            E_min_local_idx = 0;
//            E_max_local_idx = math::binary_search(worker.E_bins_total.data(),  worker.E_bins_total(E_idx_high) + overlap_range, worker.E_bins_total.size());
//            worker.E_min_local = worker.E_bins_total(E_min_local_idx);
//            worker.E_max_local = worker.E_bins_total(E_max_local_idx);
//
//        } else if (worker.world_ID == worker.world_size - 1) {
//            E_min_local_idx = math::binary_search(worker.E_bins_total.data(), worker.E_bins_total(E_idx_low) - overlap_range, worker.E_bins_total.size());
//            E_max_local_idx = (int)worker.E_bins_total.size() - 1;
//            worker.E_min_local = worker.E_bins_total(E_min_local_idx);
//            worker.E_max_local = worker.E_bins_total(E_max_local_idx);
//
//
//        } else {
//            E_min_local_idx = math::binary_search(worker.E_bins_total.data(), worker.E_bins_total(E_idx_low) - overlap_range, worker.E_bins_total.size());
//            E_max_local_idx = math::binary_search(worker.E_bins_total.data(), worker.E_bins_total(E_idx_high) + overlap_range, worker.E_bins_total.size());
//            worker.E_min_local = worker.E_bins_total(E_min_local_idx) ;
//            worker.E_max_local = worker.E_bins_total(E_max_local_idx) ;
//             }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (worker.world_ID == 0 && debug_divide) {
//            cout << "...OK. Inherit from dos_total ";
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        worker.M_min_local = worker.M_min_global;
//        worker.M_max_local = worker.M_max_global;
//
//        //Inherit the corresponding part of the dos
//        int from = E_min_local_idx;
//        int rows = E_max_local_idx - E_min_local_idx;
//
//        if (from+rows >= worker.E_bins_total.size()){
//            cout << "TOO MANY ROWS "<< endl;
//            cout << worker << endl;
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::seconds(1));
//            exit(15);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        if (E_max_local_idx <= E_min_local_idx){
//            cout << "Local range backwards! "<< endl;
//            cout << worker << endl;
//            cout.flush();
//            std::this_thread::sleep_for(std::chrono::seconds(1));
//            exit(16);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        worker.dos         = worker.dos_total.middleRows(from, rows);
//        worker.E_bins      = worker.E_bins_total.segment(from, rows);
//        worker.in_window   = worker.check_in_window(worker.E);
//        if (worker.in_window){
//            worker.E_idx = math::binary_search(worker.E_bins.data(), worker.E, worker.E_bins.size());
//            worker.M_idx = math::binary_search(worker.M_bins.data(), worker.M, worker.M_bins.size());
//        }
//        worker.histogram.conservativeResize(worker.dos.rows(), worker.dos.cols());
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
//            std::this_thread::sleep_for(std::chrono::microseconds(10000));
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
