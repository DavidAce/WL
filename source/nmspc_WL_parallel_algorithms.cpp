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
        timer::swap = 0;
        worker.t_swap.tic();

        int abort;
        MPI_Allreduce(&worker.need_to_resize_global, &abort, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        if (abort || counter::vol_merges < 1) {
            worker.t_swap.toc();
            return;
        }
        MPI_Allreduce(&worker.state_in_window, &abort, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

        if (abort || counter::vol_merges < 1) {
            worker.t_swap.toc();
            return;
        }
        int swap;
        double dos_X, dos_Y;
        double E_X, E_Y, M_X, M_Y;
        int E_X_idx, E_Y_idx, M_X_idx, M_Y_idx;
        double   E_min_up, E_max_up;
        double P_swap;      //Swap probability
        bool myTurn = math::mod(worker.world_ID, 2) == math::mod(counter::swaps, 2);

        int up = math::mod(worker.world_ID + 1, worker.world_size);
        int dn = math::mod(worker.world_ID - 1, worker.world_size);
//        if (debug_swap) {
//            for (int w = 0; w < worker.world_size; w++) {
//                if (w == worker.world_ID) {
//                    cout << "ID: " << w << " Starting Swap. Myturn = " << myTurn << endl;
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//            }
//        }
        MPI_Request reqs[4];
        MPI_Status stats[4];
        if (myTurn) {
            //Receive relevant info
            MPI_Sendrecv(&worker.E, 1, MPI_DOUBLE,up,100, &E_Y, 1,MPI_DOUBLE,up, 100,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.M, 1, MPI_DOUBLE,up,101, &M_Y, 1,MPI_DOUBLE,up, 101,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            E_Y_idx = math::binary_search_nearest(worker.E_bins, E_Y);
            M_Y_idx = math::binary_search_nearest(worker.M_bins, M_Y);
            MPI_Irecv(&E_min_up ,1, MPI_DOUBLE, up, 102, MPI_COMM_WORLD, &reqs[0]);
            MPI_Irecv(&E_max_up ,1, MPI_DOUBLE, up, 103, MPI_COMM_WORLD, &reqs[1]);
            MPI_Irecv(&dos_Y,    1, MPI_DOUBLE, up, 104, MPI_COMM_WORLD, &reqs[2]);
            MPI_Irecv(&dos_X,    1, MPI_DOUBLE, up, 105, MPI_COMM_WORLD, &reqs[3]);
            MPI_Waitall(4,reqs,stats);


            //Swap if inside the upper window and vice versa
            swap = 1;
            swap = swap && worker.E >= E_min_up && E_Y <= worker.E_max_local;
            swap = swap && worker.E <= E_max_up && E_Y >= worker.E_min_local;
            swap = swap && worker.check_in_window(worker.E);
            swap = swap && E_Y <= E_max_up && E_Y >= E_min_up;
            //Swap with probability P_swap
            if (swap) {
                P_swap = exp(worker.dos(worker.E_idx, worker.M_idx)
                             - worker.dos(E_Y_idx, M_Y_idx)
                             + dos_Y
                             - dos_X);
            } else { P_swap = 0; }
            swap = swap && rn::uniform_double_1() < fmin(1, P_swap);
            //Make sure the last worker doesn't swap!
            swap = swap && worker.world_ID != worker.world_size - 1;
            MPI_Send(&swap, 1, MPI_INT, up, 106, MPI_COMM_WORLD);

        } else {
            MPI_Sendrecv(&worker.E,1,MPI_DOUBLE,dn,100, &E_X, 1,MPI_DOUBLE,dn, 100,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.M,1,MPI_DOUBLE,dn,101, &M_X, 1,MPI_DOUBLE,dn, 101,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            E_X_idx = math::binary_search_nearest(worker.E_bins, E_X);
            M_X_idx = math::binary_search_nearest(worker.M_bins, M_X);
            MPI_Isend(&worker.E_min_local,                     1, MPI_DOUBLE, dn, 102, MPI_COMM_WORLD, &reqs[0]);
            MPI_Isend(&worker.E_max_local,                     1, MPI_DOUBLE, dn, 103, MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(&worker.dos(worker.E_idx, worker.M_idx), 1, MPI_DOUBLE, dn, 104, MPI_COMM_WORLD, &reqs[2]);
            MPI_Isend(&worker.dos(E_X_idx     , M_X_idx     ), 1, MPI_DOUBLE, dn, 105, MPI_COMM_WORLD, &reqs[3]);
            MPI_Waitall(4,reqs,stats);

            MPI_Recv(&swap, 1, MPI_INT, dn, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
//
//        if (debug_swap) {
//            for (int w = 0; w < worker.world_size; w++) {
//                if (w == worker.world_ID) {
//                    cout << "ID: " << w << " Swap = " << swap << endl;
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//            }
//        }
//
//
//        if (debug_swap) { debug_print(worker, "Starting swap\n"); }

        //Now do the swapping if you got lucky

        if (myTurn) {
            if (swap == 1) {
                MPI_Sendrecv_replace(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, up,
                                     111, up, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                worker.E = E_Y;
                worker.M = M_Y;
                worker.state_is_valid = false;
            }
        } else {
            if (swap == 1) {
                MPI_Sendrecv_replace(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, dn,
                                     111, dn, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                worker.E = E_X;
                worker.M = M_X;
                worker.state_is_valid = false;
            }
        }

        if (debug_swap) { debug_print(worker, "Swap OK\n"); }

        counter::swap_accepts += swap;
        counter::swaps++;
        worker.t_swap.toc();
    }

    void merge(class_worker &worker, bool broadcast, bool setNaN) {
        if (debug_merge) { debug_print(worker, "\nMerging. "); }
        ArrayXXd dos_total;
        ArrayXd E_total;
        ArrayXd M_total;
        //Start by trimming


        if (worker.world_ID == 0) {
            dos_total = worker.dos;
            E_total = worker.E_bins;
            M_total = worker.M_bins;
        }


        for (int w = 1; w < worker.world_size; w++) {
            if (worker.world_ID == 0) {
                if (debug_merge) {
                    cout << "Receiving from " << w << endl;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
                ArrayXXd dos_recv;
                ArrayXd E_recv, M_recv;
                recv_dynamic(dos_recv, MPI_DOUBLE, w);
                recv_dynamic(E_recv, MPI_DOUBLE, w);
                recv_dynamic(M_recv, MPI_DOUBLE, w);
                //Find average height of shared sections
                int E_shared_low_idx = math::binary_search_nearest(E_total, fmax(E_recv.minCoeff(), E_total.minCoeff()));
                int E_shared_high_idx = math::binary_search_nearest(E_total, fmin(E_recv.maxCoeff(), E_total.maxCoeff()));
                int E_shared_low_up_idx = math::binary_search_nearest(E_recv, fmax(E_recv.minCoeff(), E_total.minCoeff()));
                int E_shared_high_up_idx = math::binary_search_nearest(E_recv, fmin(E_recv.maxCoeff(), E_total.maxCoeff()));

                double diff = math::dos_distance(dos_total.middleRows(E_shared_low_idx, E_shared_high_idx - E_shared_low_idx),
                                                 dos_recv.middleRows(E_shared_low_up_idx, E_shared_high_up_idx - E_shared_low_up_idx));

                //Add that difference to the next
                if (std::isnan(diff)) {
                    cout << setprecision(2) << std::fixed << std::showpoint;
                    cout << "Tried to concatenate " << w - 1 << " and " << w << ". Diff between two DOS is NaN, exiting" << endl;
                    cout << "dos_tot" << endl << dos_total << endl;
                    cout << "dos_rec" << endl << dos_recv << endl;
                    exit(1);
                }

                math::add_to_nonzero_nonnan(dos_recv, diff);
                //Now now all doses should be the same height

                if (debug_merge) {
                    int zero_rows = 0;
                    int rows_total = E_shared_high_idx - E_shared_low_idx;
                    int rows_total_up = E_shared_high_up_idx - E_shared_low_up_idx;
                    for (int i = 0; i < dos_recv.rows(); i++) {
                        if ((dos_recv.row(i) == 0).all()) {
                            zero_rows += 1;
                        }
                    }
                    cout << setprecision(2) << fixed;
                    if (rows_total != rows_total_up) {
                        cout << "Rows mismatch!!" << endl;
                    }
                    cout << "Concatenating " << w - 1 << " and " << w << endl;
                    cout << "Zero rows recv = " << zero_rows << endl;
                    cout << "E_total = " << E_total.size()
                         << " dos_total =  " << dos_total.rows()
                         << " Shared Rows_total = " << rows_total
                         << " E_idx_low = " << E_shared_low_idx
                         << " E_idx_high = " << E_shared_high_idx
                         << endl;
                    cout << "E_recv  = " << E_recv.size()
                         << " dos_recv  =  " << dos_recv.rows()
                         << " Shared Rows_up    = " << rows_total_up
                         << " E_idx_low = " << E_shared_low_up_idx
                         << " E_idx_high = " << E_shared_high_up_idx
                         << endl;
//                    cout << "dos_total:" << endl << dos_total << endl << endl;
//                    cout << "dos_recv:" << endl << dos_recv<< endl << endl;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }

                double weight;
                double E_span = fmax(E_total(E_shared_high_idx) - E_total(E_shared_low_idx), E_recv(E_shared_high_up_idx) - E_recv(E_shared_low_up_idx));
                vector<ArrayXd> dos_merge;
                vector<double> E_merge;
                if (debug_merge) { cout << "Merging " << w - 1 << " and " << w << endl; }
                for (int i = 0; i < E_total.size(); i++) {
                    if (i < E_shared_low_idx) {
                        //Take E_total(i)
                        dos_merge.push_back(dos_total.row(i));
                        E_merge.push_back(E_total(i));
                        if (debug_merge) {
                            cout << " Inserted " << "i: " << i << endl
                                 << " E_total(i): " << E_total(i) << endl
                                 << " E_total   : " << E_total.transpose() << endl
                                 << " E_recv    : " << E_recv.transpose() << endl
                                 << " E_merge   : " << E_merge << endl
                                 << endl;
                        }
                        continue;
                    } else if (i <= E_shared_high_idx) {
                        int j = math::binary_search_exact(E_recv, E_total(i));
                        if (j == -1) {
                            printf(" Could not find E_total(%d) = %f. Closest match: E_recv(%d) = %f \n", i, E_total(i), j, E_recv(j));
                            exit(1);
                        } else {
                            //Merge taking Average
                            weight = (E_total(i) - E_total(E_shared_low_idx)) / E_span;
                            dos_merge.push_back((1 - weight) * dos_total.row(i) + weight * dos_recv.row(j));
                            E_merge.push_back(E_total(i));
                            if (debug_merge) {
                                cout << " Merged  i : " << i << " and j: " << j << endl
                                     << " E_total(i): " << E_total(i) << endl
                                     << " E_total   : " << E_total.transpose() << endl
                                     << " E_recv    : " << E_recv.transpose() << endl
                                     << " E_merge   : " << E_merge << endl
                                     << endl;
                            }

                        }
                    }
                }

                for (int j = 0; j < E_recv.size(); j++) {
                    //Check that it hasn't been used already
                    if (j <= E_shared_high_up_idx) { continue; }
                    //Take E_recv(j)
                    dos_merge.push_back(dos_recv.row(j));
                    E_merge.push_back(E_recv(j));
                    if (debug_merge) {
                        cout << " Inserted " << "j: " << j << endl
                             << " E_recv    : " << E_recv.transpose() << endl
                             << " E_merge   : " << E_merge << endl
                             << endl;
                    }
                }

                dos_total.resize(dos_merge.size(), M_recv.size());
                E_total.resize(E_merge.size());
                for (unsigned int i = 0; i < dos_merge.size(); i++) {
                    dos_total.row(i) = dos_merge[i];
                    E_total(i) = E_merge[i];
                }
                dos_merge.clear();
                E_merge.clear();
                if (debug_merge) {
                    cout << "OK ";
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(1000));
                }
            } else if (w == worker.world_ID) {
                send_dynamic(worker.dos, MPI_DOUBLE, 0);
                send_dynamic(worker.E_bins, MPI_DOUBLE, 0);
                send_dynamic(worker.M_bins, MPI_DOUBLE, 0);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        if (setNaN) {
            if (worker.world_ID == 0) {
                math::subtract_min_nonzero_nan(dos_total);
            }
        }
        if (worker.world_ID == 0) {
            //Let the master worker have the results
            worker.dos_total = dos_total;
            worker.E_bins_total = E_total;
            worker.M_bins_total = M_total;
        }
        if (broadcast) {
            broadcast_merger(worker);
        }
    }

    void broadcast_merger(class_worker &worker) {
        if (debug_bcast) { debug_print(worker, "\nBroadcasting raw to workers "); }
        mpi::bcast_dynamic(worker.dos_total, MPI_DOUBLE, 0);
        mpi::bcast_dynamic(worker.E_bins_total, MPI_DOUBLE, 0);
        mpi::bcast_dynamic(worker.M_bins_total, MPI_DOUBLE, 0);
        if (debug_bcast) { debug_print(worker, "Broadcast OK! \n "); }
    }

    void divide_global_range_dos_area(class_worker &worker) {
        //Update limits
        if (debug_divide) { debug_print(worker, "\n Dividing Area. "); }
        double global_range = math::area(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
        double local_range = global_range / worker.world_size;
        double x = constants::overlap_factor_dos_area / (1 - constants::overlap_factor_dos_area / 2);
        double overlap_range = local_range * x;//2.0*(worker.world_size - 2.0 + x)/worker.world_size;
        //The overlap_range is the total range in a domain that will have overlap, either up or down.
        //Find the boundaries of the DOS domain that gives every worker the  same DOS volume to work on
        int E_min_local_idx, E_max_local_idx;
        if (worker.world_ID == 0) { overlap_range = overlap_range / 2; }
        else if (worker.world_ID == worker.world_size - 1) { overlap_range = overlap_range / 2; }
        else { overlap_range = overlap_range / 4; }
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
        if (debug_divide) { debug_print(worker, " Inheriting from dos_total. "); }

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


        worker.dos                  = worker.dos_total.middleRows(from, rows);
        worker.E_bins               = worker.E_bins_total.segment(from, rows);
        worker.histogram            = ArrayXXi::Zero(worker.dos.rows(), worker.dos.cols());
        worker.random_walk.clear();
        worker.saturation.clear();
        worker.state_is_valid = false;

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
        if (debug_divide) { debug_print(worker, " Merged OK. \n "); }

    }

    void divide_global_range_dos_volume(class_worker &worker) {

        if (debug_divide) { debug_print(worker, " \nDividing dos volume "); }

        double global_range = math::volume(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
        double local_range = global_range / worker.world_size;
        double x = constants::overlap_factor_dos_vol / (1 - constants::overlap_factor_dos_vol / 2);
        double overlap_range = local_range * x;
        //The overlap_range is the total range in a domain that will have overlap, either up or down.
        //Find the boundaries of the DOS domain that gives every worker the  same DOS volume to work on
        int E_min_local_idx, E_max_local_idx;
        if       (worker.world_ID == 0)                    { overlap_range = overlap_range / 2; }
        else if (worker.world_ID == worker.world_size - 1) { overlap_range = overlap_range / 2; }
        else                                               { overlap_range = overlap_range / 4; }
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
        if (debug_divide) { debug_print(worker, " Inheriting from dos_total "); }

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


        worker.dos                  = worker.dos_total.middleRows(from, rows);
        worker.E_bins               = worker.E_bins_total.segment(from, rows);
        worker.histogram            = ArrayXXi::Zero(worker.dos.rows(), worker.dos.cols());
        worker.random_walk.clear();
        worker.saturation.clear();
        worker.state_is_valid = false;

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
        if (debug_divide) { debug_print(worker, " Merge OK.\n "); }

    }


    void take_help(class_worker &worker) {
        //Check if anybody was able to do a walk
        if (worker.help.active) {
            timer::take_help = 0;
            timer::add_dos   = 0;
            worker.t_help.tic();
            int size = (int) worker.random_walk.size();
            ArrayXi sizes(worker.help.help_size);
            ArrayXi displ(worker.help.help_size);
            MPI_Allgather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, worker.help.MPI_COMM_HELP);
            if ((sizes != sizes(0)).any() ) {
                if(debug_take_help) {
                    if (worker.help.help_rank == 0) {
                        cout << "Random walk size mismatch! "
                             << sizes.transpose()
                             << endl;
                    }
                }
                worker.random_walk.clear();
                worker.t_help.toc();
                return;
            }

            vector<state> random_walk_recv((unsigned long) (sizes.sum()));
            MPI_Allgather(worker.random_walk.data(), size, MPI_2INT,  random_walk_recv.data(), size, MPI_2INT, worker.help.MPI_COMM_HELP);

            for (int i = 0; i < random_walk_recv.size(); i++) {
                worker.histogram(random_walk_recv[i].E_idx, random_walk_recv[i].M_idx) += 1;
                worker.dos(random_walk_recv[i].E_idx, random_walk_recv[i].M_idx) += worker.lnf;
            }
            counter::MCS            += constants::rate_take_help * (worker.help.help_size - 1);
            timer::add_hist_volume  += constants::rate_take_help * (worker.help.help_size - 1);
            timer::check_saturation += constants::rate_take_help * (worker.help.help_size - 1);
            worker.random_walk.clear();
            worker.t_help.toc();
        }
    }

    void setup_help(class_worker &worker, class_backup &backup) {
        timer::setup_help = 0;
        worker.t_help_setup.tic();

        //Find out who's finished
        int any_finished = worker.finish_line;
        MPI_Allreduce(MPI_IN_PLACE, &any_finished, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        //if nobody is finished, return.
        if (any_finished == 0) {
            worker.t_help_setup.toc();
            return;
        }
        if (debug_setup_help) { debug_print(worker, "Setting up help. "); }

        //Now all the help.available workers need to be set up again. Two scenarios, either they've (1) helped before,
        //or they've (2) never helped.
        //If (1), just take_help() and return.
        //If (2), setup from scratch.

        ArrayXi whos_helping_who_old = ArrayXi::Constant(worker.world_size, -1);  //For comparison!
        ArrayXi whos_helping_who_new = ArrayXi::Constant(worker.world_size, -1);  //For comparison!
        MPI_Allgather(&worker.help.helping_id, 1, MPI_INT, whos_helping_who_old.data(), 1, MPI_INT, MPI_COMM_WORLD);
        worker.help.reset();
        ArrayXi all_available(worker.world_size);
        MPI_Allgather(&worker.finish_line, 1, MPI_INT, all_available.data(), 1, MPI_INT, MPI_COMM_WORLD);

        //Now find out who to help next. compute max amount of helpers per worker.
        int max_helpers = (int) ceil((double) all_available.sum() / (worker.world_size - all_available.sum()));
        ArrayXi given_help(worker.world_size);
        given_help.fill(0);
        for (int w = 0; w < worker.world_size; w++) {
            if (w == worker.world_ID) {
                if (worker.finish_line) {
                    for (int i = 0; i < worker.world_size; i++) {
                        if (all_available(i) == 0 && given_help(i) < max_helpers) {
                            given_help(i)++;
                            worker.help.giving_help = true;
                            worker.help.helping_id  = i;
                            break;
                        }
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(given_help.data(), worker.world_size, MPI_INT, w, MPI_COMM_WORLD);
        }
        //If you got help, you are now rank (key) 0 in your own communicator (whose color is your ID).
        //Otherwise, if you are giving help you join color = helping_id, and pass key = color+1?

        worker.help.getting_help = given_help(worker.world_ID) > 0;
        int color, key;
        if (worker.help.getting_help) {
            color               = worker.world_ID;
            key                 = 0;
            worker.help.active  = true;
        } else if (worker.help.giving_help) {
            color               = worker.help.helping_id;
            key                 = color + 1;
            worker.help.active  = true;
        } else {
            color               = MPI_UNDEFINED;
            key                 = 0;
            worker.help.active  = false;
        }
        MPI_Comm_split(MPI_COMM_WORLD, color, key, &worker.help.MPI_COMM_HELP);
        MPI_Allgather(&worker.help.helping_id, 1, MPI_INT, whos_helping_who_new.data(), 1, MPI_INT, MPI_COMM_WORLD);


        //If there has been no change since last time, simply return;
        if (whos_helping_who_new.cwiseEqual(whos_helping_who_old).all()) {
            worker.t_help_setup.toc();
            return;
        }
        //If you are not getting and not giving help, just leave
        if (worker.help.MPI_COMM_HELP == MPI_COMM_NULL) {
            backup.restore_state(worker);
            worker.t_help_setup.toc();
            return;
        }
        cout <<"ID: " << worker.world_ID << ": New help!" << endl;
        MPI_Comm_rank(worker.help.MPI_COMM_HELP, &worker.help.help_rank);
        MPI_Comm_size(worker.help.MPI_COMM_HELP, &worker.help.help_size);

        //Now every available guy knows who to help (rank 0 in MPI_COMM_HELP) and the helpees know who to send their info to (Bcast).
        //Start by backing up and such
        if (worker.help.giving_help) { backup.backup_state(worker); }

        //Send out the details for help to begin.
        mpi::bcast_dynamic(worker.dos,       MPI_DOUBLE, 0, worker.help.MPI_COMM_HELP);
        mpi::bcast_dynamic(worker.histogram, MPI_INT   , 0, worker.help.MPI_COMM_HELP);
        mpi::bcast_dynamic(worker.E_bins,    MPI_DOUBLE, 0, worker.help.MPI_COMM_HELP);
        mpi::bcast_dynamic(worker.M_bins,    MPI_DOUBLE, 0, worker.help.MPI_COMM_HELP);
        MPI_Bcast(&worker.lnf,     1,        MPI_DOUBLE, 0, worker.help.MPI_COMM_HELP);
        MPI_Bcast(&counter::MCS,   1,        MPI_INT,    0, worker.help.MPI_COMM_HELP);
        MPI_Bcast(&counter::walks, 1,        MPI_INT,    0, worker.help.MPI_COMM_HELP);
        ArrayXi saturation_map = Map<ArrayXi>(worker.saturation.data(), (int)worker.saturation.size());
        mpi::bcast_dynamic(saturation_map,   MPI_INT   , 0, worker.help.MPI_COMM_HELP);
        if (worker.help.giving_help) {
            worker.E_min_local      = worker.E_bins.minCoeff();
            worker.E_max_local      = worker.E_bins.maxCoeff();
            worker.M_min_local      = worker.M_bins.minCoeff();
            worker.M_max_local      = worker.M_bins.maxCoeff();
        }
        worker.set_rate_increment();
        worker.saturation.clear();
        worker.random_walk.clear();
        worker.state_is_valid = false;

        timer::increment            = 0;
        worker.t_help_setup.toc();
    }
}


//void check_help(class_worker &worker) {
//    //Check if anybody was able to do a walk
//    timer::check_help = 0;
//    if (worker.help.MPI_COMM_HELP != MPI_COMM_NULL) {
//        struct {
//            int highest_walk;
//            int highest_idx;
//        } in, out;
//        in.highest_walk = counter::walks;
//        in.highest_idx  = worker.help.help_rank;
//        MPI_Allreduce(&in, &out, 1, MPI_2INT, MPI_MAXLOC, worker.help.MPI_COMM_HELP);
//        if (out.highest_walk > worker.help.help_walks) {
//            //Somebody made progress, then get everybody up to speed
//            //Also, broadcast an updated dos
//            if (debug_take_help) {
//                if (out.highest_idx == worker.help.help_rank) {
//                    cout << "ID: " << worker.world_ID
//                         << " Broadcast from " << out.highest_idx
//                         << " (helping: " << worker.help.helping_id << ")"
//                         << " MCS: " << counter::MCS
//                         << " Slope: " << worker.slope
//                         << " Size: " << worker.help.help_size
//                         << " Walks: " << counter::walks
//                         << " Highest walks: " << out.highest_walk
//                         << " Current walks: " << worker.help.help_walks
//                         << endl;
//                }
//            }
//            MPI_Bcast(worker.dos.data(), (int) worker.dos.size(), MPI_DOUBLE, out.highest_idx, worker.help.MPI_COMM_HELP);
//            worker.help.help_walks = out.highest_walk;
//            for (int i = 0; i < out.highest_walk - counter::walks; i++) {
//                worker.next_WL_iteration();
//                if (debug_take_help) {
//                    cout << "ID: " << worker.world_ID
//                         << " moved to " << counter::walks
//                         << endl;
//                }
//            }
//            worker.histogram.setZero();
//            worker.saturation.clear();
//            worker.random_walk.clear();
//        }
//    }
//}

//    void check_help(class_worker &worker) {
//        //Check if anybody was able to do a walk
//        timer::check_help = 0;
//        if (worker.help.MPI_COMM_HELP != MPI_COMM_NULL) {
//            struct {
//                int highest_walk;
//                int highest_idx;
//            } in, out;
//            in.highest_walk = counter::walks;
//            in.highest_idx  = worker.help.help_rank;
//            MPI_Allreduce(&in, &out, 1, MPI_2INT, MPI_MAXLOC, worker.help.MPI_COMM_HELP);
//            if (out.highest_walk > worker.help.help_walks) {
//                //Somebody made progress, then get everybody up to speed
//                //Also, broadcast an updated dos
//                if (debug_take_help) {
//                    if (out.highest_idx == worker.help.help_rank) {
//                        cout << "ID: " << worker.world_ID
//                             << " Broadcast from " << out.highest_idx
//                             << " (helping: " << worker.help.helping_id << ")"
//                             << " MCS: " << counter::MCS
//                             << " Slope: " << worker.slope
//                             << " Size: " << worker.help.help_size
//                             << " Walks: " << counter::walks
//                             << " Highest walks: " << out.highest_walk
//                             << " Current walks: " << worker.help.help_walks
//                             << endl;
//                    }
//                }
//                MPI_Bcast(worker.dos.data(), (int) worker.dos.size(), MPI_DOUBLE, out.highest_idx, worker.help.MPI_COMM_HELP);
//                worker.help.help_walks = out.highest_walk;
//                for (int i = 0; i < out.highest_walk - counter::walks; i++) {
//                    worker.next_WL_iteration();
//                    if (debug_take_help) {
//                        cout << "ID: " << worker.world_ID
//                             << " moved to " << counter::walks
//                             << endl;
//                    }
//                }
//                worker.histogram.setZero();
//                worker.saturation.clear();
//                worker.random_walk.clear();
//            }
//        }
//    }
