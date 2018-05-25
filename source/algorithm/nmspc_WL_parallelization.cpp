//
// Created by david on 2016-08-11.
//

#include <assert.h>
#include "nmspc_WL_parallelization.h"
#include "class_WL_worker.h"

#define debug_swap      0
#define debug_merge     1
#define debug_bcast     1
#define debug_divide    1
#define debug_take_help 1
#define debug_setup_help 1
#define debug_resize_local_bins 1
using namespace std;

namespace parallel {
    void swap(class_worker &worker) {
        //Use MPI Tag in the 100-200 range
        timer::swap = 0;
        worker.t_swap.tic();

        int abort;
        MPI_Allreduce(&worker.need_to_resize_global, &abort, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (abort || counter::merges < 1) {
            worker.t_swap.toc();
            return;
        }

        int swap;
        double dos_X, dos_Y;
        double E_X, E_Y, M_X, M_Y;
        int E_X_idx, E_Y_idx, M_X_idx, M_Y_idx;
        double   E_min_up, E_max_up;
        double P_swap;      //Swap probability
//        bool myTurn  = math::mod(worker.world_ID, 2) == math::mod(counter::swaps, 2);
        int even_team = math::mod(worker.world_ID / constants::team_size , 2) == 0;
        bool myTurn   = math::mod(even_team + worker.world_ID, 2) == math::mod(counter::swaps, 2);
//        int distance  = constants::team_size; //Needs to be an odd number
        int distance = constants::team_size + math::mod(constants::team_size + 1, 2); //Needs to be an odd number
        int up = math::mod(worker.world_ID + distance, worker.world_size);
        int dn = math::mod(worker.world_ID - distance, worker.world_size);
        if (debug_swap) {
            for (int w = 0; w < worker.world_size; w++) {
                if (w == worker.world_ID) {
                    cout << "ID: " << w << " Starting Swap. Myturn = " << myTurn << " dn: " << dn << " up: " << up << endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
        MPI_Request reqs[4];
        MPI_Status stats[4];
        assert(up < worker.world_size and up >= 0 and "Neighbor up out of range");
        assert(dn < worker.world_size and dn >= 0 and "Neighbor down out of range");

        if (myTurn) {
            //Receive relevant info
            MPI_Sendrecv(&worker.E, 1, MPI_DOUBLE,up,100, &E_Y, 1,MPI_DOUBLE,up, 100,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.M, 1, MPI_DOUBLE,up,101, &M_Y, 1,MPI_DOUBLE,up, 101,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            E_Y_idx = math::binary_search_exact(worker.E_bins, E_Y);
            M_Y_idx = math::binary_search_exact(worker.M_bins, M_Y);
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
            swap = swap && E_Y     <= E_max_up && E_Y     >= E_min_up;
            swap = swap && E_Y_idx != -1       && M_Y_idx != -1 && dos_X != -1;

            //Swap with probability P_swap
            if (swap) {
                P_swap = exp(worker.dos(worker.E_idx, worker.M_idx)
                             - worker.dos(E_Y_idx, M_Y_idx)
                             + dos_Y
                             - dos_X);
            } else { P_swap = 0; }
            swap = swap && rn::uniform_double_1() < fmin(1, P_swap);
            //Make sure the last worker doesn't swap!
            MPI_Send(&swap, 1, MPI_INT, up, 106, MPI_COMM_WORLD);

        } else {
            MPI_Sendrecv(&worker.E,1,MPI_DOUBLE,dn,100, &E_X, 1,MPI_DOUBLE,dn, 100,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&worker.M,1,MPI_DOUBLE,dn,101, &M_X, 1,MPI_DOUBLE,dn, 101,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            E_X_idx = math::binary_search_exact(worker.E_bins, E_X);
            M_X_idx = math::binary_search_exact(worker.M_bins, M_X);
            MPI_Isend(&worker.E_min_local,                     1, MPI_DOUBLE, dn, 102, MPI_COMM_WORLD, &reqs[0]);
            MPI_Isend(&worker.E_max_local,                     1, MPI_DOUBLE, dn, 103, MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(&worker.dos(worker.E_idx, worker.M_idx), 1, MPI_DOUBLE, dn, 104, MPI_COMM_WORLD, &reqs[2]);
            if (E_X_idx == -1 || M_X_idx == -1){
                double dummy = -1;
                MPI_Isend(&dummy                                 , 1, MPI_DOUBLE, dn, 105, MPI_COMM_WORLD, &reqs[3]);
            }else{
                MPI_Isend(&worker.dos(E_X_idx     , M_X_idx     ), 1, MPI_DOUBLE, dn, 105, MPI_COMM_WORLD, &reqs[3]);
            }
            MPI_Waitall(4,reqs,stats);
            MPI_Recv(&swap, 1, MPI_INT, dn, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
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

    void swap2(class_worker &worker) {
        //Use MPI Tag in the 100-200 range
        timer::swap = 0;
        worker.t_swap.tic();

        int abort;
        MPI_Allreduce(&worker.need_to_resize_global, &abort, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (abort || counter::merges < 1) {
            worker.t_swap.toc();
            return;
        }



        int swap;
//        double dos_X, dos_Y;
//        double E_X, E_Y, M_X, M_Y;
        double dos_partner, dos_myown;
        double E_partner, M_partner;
        int E_idx_partner, M_idx_partner;
//        int E_X_idx, E_Y_idx, M_X_idx, M_Y_idx;
        double   E_min_partner, E_max_partner;
        double P_swap;      //Swap probability
//        bool myTurn   = math::mod(worker.world_ID, 2) == math::mod(counter::swaps, 2);

        // Broadcast a randomly shuffled list to all workers like {3,2,7,1,0,5...}
        // Find myself on the list.
        // If I'm at an even position, it's my turn. Otherwise not.
        // If it is my turn, look at who's on my right, that's who I'm swapping with.
        // If it is not my turn, look left, that's who I'm swapping with.

        std::vector<int> swap_list((unsigned long)worker.world_size);
        if (worker.world_ID == 0){
            swap_list = rn::shuffled_list(0,worker.world_size - 1);
//            std::cout << "Swap list: " << swap_list << std::endl;
        }
        MPI_Bcast(swap_list.data(), worker.world_size, MPI_INT, 0, MPI_COMM_WORLD);
        int swap_partner;
        int my_idx_on_swap_list = (int) (std::find(swap_list.begin(), swap_list.end(), worker.world_ID) - swap_list.begin());
        int myTurn =  math::mod(my_idx_on_swap_list, 2); // It's my turn to throw dice if I'm at an even position on the swap_list
        if (myTurn){
            swap_partner = swap_list[math::mod(my_idx_on_swap_list + 1 , worker.world_size)];
        }else {
            swap_partner = swap_list[math::mod(my_idx_on_swap_list - 1 , worker.world_size)];
        }

        assert(swap_partner < worker.world_size and swap_partner >= 0 and "Neighbor swap_partner out of range");
        if (debug_swap) {
            for (int w = 0; w < worker.world_size; w++) {
                if (w == worker.world_ID) {
                    cout << "ID: " << w << " Starting Swap. Myturn = " << myTurn << " swap_partner: " << swap_partner << endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }

        //Send each others info
        MPI_Sendrecv(&worker.E, 1, MPI_DOUBLE,swap_partner,100, &E_partner, 1,MPI_DOUBLE,swap_partner, 100,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&worker.M, 1, MPI_DOUBLE,swap_partner,101, &M_partner, 1,MPI_DOUBLE,swap_partner, 101,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&worker.E_min_local, 1, MPI_DOUBLE,swap_partner,102, &E_min_partner, 1,MPI_DOUBLE,swap_partner, 102,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&worker.E_max_local, 1, MPI_DOUBLE,swap_partner,103, &E_max_partner, 1,MPI_DOUBLE,swap_partner, 103,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&worker.dos(worker.E_idx, worker.M_idx), 1, MPI_DOUBLE,swap_partner,104, &dos_partner, 1,MPI_DOUBLE,swap_partner, 104,MPI_COMM_WORLD, MPI_STATUS_IGNORE);


        E_idx_partner = math::binary_search_exact(worker.E_bins, E_partner);
        M_idx_partner = math::binary_search_exact(worker.M_bins, M_partner);
        //If we are doing 1d random walks we don't really care about M_idx_partner and M_idx. Set these to 0.
        if(constants::rw_dims == 1){M_idx_partner = 0;}

        double dos_send_back;

        if (E_idx_partner == -1 || M_idx_partner == -1){
            dos_send_back = -1;
        }else{
            dos_send_back = worker.dos(E_idx_partner, M_idx_partner);
        }
        MPI_Sendrecv(&dos_send_back, 1, MPI_DOUBLE,swap_partner,105, &dos_myown, 1,MPI_DOUBLE,swap_partner, 105,MPI_COMM_WORLD, MPI_STATUS_IGNORE);


        //Decide if you want to swap
        //Swap if your're inside each others windows, your respective windows, and not swapping with self.
        swap = 1;
        swap = swap && worker.E >= E_min_partner && worker.E <= E_max_partner;
        swap = swap && worker.check_in_window(worker.E);
        swap = swap && E_idx_partner != -1       && M_idx_partner != -1 && dos_partner != -1 && dos_myown != -1;
        swap = swap && worker.world_ID != swap_partner;

        //If it's your turn, throw a dice!. Swap with probability P_swap
        if (myTurn) {
            if (swap) {
                P_swap = exp(worker.dos(worker.E_idx, worker.M_idx)
                             - worker.dos(E_idx_partner, M_idx_partner)
                             + dos_partner
                             - dos_myown);
            } else { P_swap = 0; }
            swap = swap && rn::uniform_double_1() < fmin(1, P_swap);
        }

        int partner_wanna_swap;
        MPI_Sendrecv(&swap, 1, MPI_INT,swap_partner,106, &partner_wanna_swap, 1,MPI_INT,swap_partner, 106,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        swap = swap && partner_wanna_swap;

        //Do swap if yo got lucky!
        if (swap == 1) {
            MPI_Sendrecv_replace(worker.model.lattice.data(), (int) worker.model.lattice.size(), MPI_INT, swap_partner,
                                 107, swap_partner, 107, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            worker.E = E_partner;
            worker.M = M_partner;
            worker.state_is_valid = false;
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
//        int num_teams = worker.world_size/constants::team_size;

        if (worker.team.team_leader){
            if (worker.team.team_commander) {
                dos_total = worker.dos;
                E_total = worker.E_bins;
                M_total = worker.M_bins;
            }
            std::cout << "ID: " << worker.world_ID << " made it here. Num tems "<<  constants::num_teams  << std::endl;
            for (int w = 1; w < constants::num_teams; w++) {
                if (worker.team.team_commander) {
                    if (debug_merge) {
                        cout << "Receiving from " << w << endl;
                        cout.flush();
                        std::this_thread::sleep_for(std::chrono::microseconds(1000));
                    }
                    ArrayXXd dos_recv;
                    ArrayXd E_recv, M_recv;
                    mpi::recv_dynamic(dos_recv, MPI_DOUBLE, w, worker.team.MPI_COMM_LEAD);
                    mpi::recv_dynamic(E_recv,   MPI_DOUBLE, w, worker.team.MPI_COMM_LEAD);
                    mpi::recv_dynamic(M_recv,   MPI_DOUBLE, w, worker.team.MPI_COMM_LEAD);
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
                        int rows_total      = E_shared_high_idx - E_shared_low_idx;
                        int rows_total_up   = E_shared_high_up_idx - E_shared_low_up_idx;
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
//                                cout << " Inserted " << "i: " << i << endl
//                                     << " E_total(i): " << E_total(i) << endl
//                                     << " E_total   : " << E_total.transpose() << endl
//                                     << " E_recv    : " << E_recv.transpose() << endl
//                                     << " E_merge   : " << E_merge << endl
//                                     << endl;
                            }
                            continue;
                        } else if (i <= E_shared_high_idx) {
                            int j = math::binary_search_exact(E_recv, E_total(i));
                            if (j == -1) {
                                j = math::binary_search_nearest(E_recv, E_total(i));
                                printf(" Could not find E_total(%d) = %f. Closest match: E_recv(%d) = %f \n", i, E_total(i), j, E_recv(j));
//                                cout << " E_total(i): " << E_total(i) << endl
//                                     << " E_total   : " << E_total.transpose() << endl
//                                     << " E_recv    : " << E_recv.transpose() << endl
//                                     << " E_merge   : " << E_merge << endl
//                                     << endl;
                                exit(1);
                            } else {
                                //Merge taking Average
                                weight = (E_total(i) - E_total(E_shared_low_idx)) / E_span;
                                dos_merge.push_back((1 - weight) * dos_total.row(i) + weight * dos_recv.row(j));
                                E_merge.push_back(E_total(i));
                                if (debug_merge) {
//                                    cout << " Merged  i : " << i << " and j: " << j << endl
//                                         << " E_total(i): " << E_total(i) << endl
//                                         << " E_total   : " << E_total.transpose() << endl
//                                         << " E_recv    : " << E_recv.transpose() << endl
//                                         << " E_merge   : " << E_merge << endl
//                                         << endl;
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
//                            cout << " Inserted " << "j: " << j << endl
//                                 << " E_recv    : " << E_recv.transpose() << endl
//                                 << " E_merge   : " << E_merge << endl
//                                 << endl;
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
                } else if (w == worker.team.team_id) {
                    mpi::send_dynamic(worker.dos, MPI_DOUBLE, 0,worker.team.MPI_COMM_LEAD);
                    mpi::send_dynamic(worker.E_bins, MPI_DOUBLE, 0, worker.team.MPI_COMM_LEAD);
                    mpi::send_dynamic(worker.M_bins, MPI_DOUBLE, 0, worker.team.MPI_COMM_LEAD);
                }
                MPI_Barrier(worker.team.MPI_COMM_LEAD);
            }

        }

        if (setNaN) {
            if (worker.team.team_commander) {
                math::subtract_min_nonzero_nan(dos_total);
            }
        }
        if (worker.team.team_commander) {
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
        if(worker.team.team_leader){
            mpi::bcast_dynamic(worker.dos_total, MPI_DOUBLE, 0, worker.team.MPI_COMM_LEAD);
            mpi::bcast_dynamic(worker.E_bins_total, MPI_DOUBLE, 0, worker.team.MPI_COMM_LEAD);
            mpi::bcast_dynamic(worker.M_bins_total, MPI_DOUBLE, 0, worker.team.MPI_COMM_LEAD);
        }
        mpi::bcast_dynamic(worker.dos_total, MPI_DOUBLE, 0, worker.team.MPI_COMM_TEAM);
        mpi::bcast_dynamic(worker.E_bins_total, MPI_DOUBLE, 0, worker.team.MPI_COMM_TEAM);
        mpi::bcast_dynamic(worker.M_bins_total, MPI_DOUBLE, 0, worker.team.MPI_COMM_TEAM);
        if (debug_bcast) { debug_print(worker, "Broadcast OK! \n "); }
    }

    void divide_global_range_uniform(class_worker &worker) {
        //Update limits
//        int num_teams = worker.world_size / constants::team_size;
        double global_range = fabs(worker.E_max_global - worker.E_min_global);
        double local_range = global_range / (double) (constants::num_teams);
        double x = constants::overlap_factor_energy;
        //Add a little bit extra if there are too many workers (Add nothing if world_size == 2, and up to local_volume if world_size == inf)
        double overlap_range = local_range * x;// * 2.0*(world_size - 2.0 + x)/world_size;
        //Make sure the overlap covers at least a hefty percent of the global range.
        //This is needed when you have too many workers on a small global energy range.
        if (overlap_range * 2 + local_range < global_range * 0.2) {
            if (worker.world_ID == 0) {
                cout << "Forced readjustment of local energy range" << endl;
            }
            overlap_range = (global_range * 0.2 - local_range) / 2;
        }

        int team_id;
        if(worker.team.active) {
            team_id = worker.team.team_id;
        }else{
            team_id = worker.world_ID /constants::team_size;

        }
        if (team_id == 0)                              { overlap_range = overlap_range / 1; }
        else if (team_id == constants::num_teams - 1)  { overlap_range = overlap_range / 1; }
        else                                           { overlap_range = overlap_range / 2; }

        worker.E_min_local = worker.E_min_global + (team_id * local_range) - overlap_range;
        worker.E_max_local = worker.E_min_global + (team_id + 1) * local_range + overlap_range;
        worker.E_max_local = fmin(worker.E_max_local, worker.E_max_global);
        worker.E_min_local = fmax(worker.E_min_local, worker.E_min_global);
        worker.M_min_local = worker.M_min_global;
        worker.M_max_local = worker.M_max_global;
        if (debug_divide) {
            if (worker.world_ID == 0) {
                cout << "Dividing according to energy range" << endl;
                cout << flush;
                std::this_thread::sleep_for(std::chrono::microseconds(10));
            }
            MPI_Barrier(MPI_COMM_WORLD);
            for (int w = 0; w < worker.world_size; w++) {
                if (w == worker.world_ID) {
                    cout << "ID: " << worker.world_ID << endl;
                    cout << "   Global Bound= " << worker.E_min_global << " " << worker.E_max_global << endl;
                    cout << "   Local Range = " << local_range << " overlap_range = " << overlap_range << endl;
                    cout << "   Bounds E: [" << worker.E_min_local << " " << worker.E_max_local << "]" << endl;
                    cout << "          M: [" << worker.M_min_local << " " << worker.M_max_local << "]" << endl;
                    std::this_thread::sleep_for(std::chrono::microseconds(10));

                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }

    }

    void divide_global_range_dos_area(class_worker &worker) {
        //Update limits
        if (debug_divide) { debug_print(worker, "\n Dividing Area. "); }
//        int num_teams           = worker.world_size / constants::team_size;
        double global_range     = math::area(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
        double local_range      = global_range / constants::num_teams;
        double x                = constants::overlap_factor_dos_area / (1 - constants::overlap_factor_dos_area / 2);
        double overlap_range = local_range * x;
        //The overlap_range is the total range in a domain that will have overlap, either up or down.
        //Find the boundaries of the DOS domain that gives every worker the  same DOS volume to work on
        int E_min_local_idx, E_max_local_idx;
        if      (worker.team.team_id == 0)                        { overlap_range = overlap_range / 2; }
        else if (worker.team.team_id == constants::num_teams - 1) { overlap_range = overlap_range / 2; }
        else                                                      { overlap_range = overlap_range / 4; }
        E_min_local_idx = math::area_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                           worker.team.team_id * local_range - overlap_range);
        E_max_local_idx = math::area_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                           (worker.team.team_id + 1) * local_range + overlap_range);

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

        //Inherit the corresponding part of the dos
        int from = E_min_local_idx;
        int rows = E_max_local_idx - E_min_local_idx + 1;

        worker.dos                  = worker.dos_total.middleRows(from, rows);
        worker.E_bins               = worker.E_bins_total.segment(from, rows);
        worker.histogram            = ArrayXXi::Zero(worker.dos.rows(), worker.dos.cols());
        worker.random_walk.clear();
        worker.saturation.clear();
        worker.state_is_valid = false;

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
    }

    void divide_global_range_dos_volume(class_worker &worker) {
        if (debug_divide) { debug_print(worker, " \nDividing dos volume "); }
//        int num_teams           = worker.world_size / constants::team_size;
        double global_range     = math::volume(worker.dos_total, worker.E_bins_total, worker.M_bins_total);
        double local_range      = global_range / constants::num_teams;
        double x                = constants::overlap_factor_dos_vol / (1 - constants::overlap_factor_dos_vol / 2);
        double overlap_range    = local_range * x;
        //The overlap_range is the total range in a domain that will have overlap, either up or down.
        //Find the boundaries of the DOS domain that gives every worker the  same DOS volume to work on
        int E_min_local_idx, E_max_local_idx;
        if      (worker.team.team_id == 0)                          { overlap_range = overlap_range / 2 ; }
        else if (worker.team.team_id == constants::num_teams - 1)   { overlap_range = overlap_range / 2; }
        else                                                        { overlap_range = overlap_range / 4; }
        E_min_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                           worker.team.team_id * local_range - overlap_range);
        E_max_local_idx = math::volume_idx(worker.dos_total, worker.E_bins_total, worker.M_bins_total,
                                           (worker.team.team_id + 1) * local_range + overlap_range);

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

        worker.dos                  = worker.dos_total.middleRows(from, rows);
        worker.E_bins               = worker.E_bins_total.segment(from, rows);
        worker.histogram            = ArrayXXi::Zero(worker.dos.rows(), worker.dos.cols());
        worker.random_walk.clear();
        worker.saturation.clear();
        worker.state_is_valid = false;

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
    }

    void resize_global_range(class_worker &worker){
        //Compare to the other workers to find global limits
        worker.E_min_global = fmin(worker.E_min_local, worker.E_min_global);
        worker.E_max_global = fmax(worker.E_max_local, worker.E_max_global);
        switch (constants::rw_dims) {
            case 1:
                MPI_Allreduce(MPI_IN_PLACE, &worker.E_min_global,  1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &worker.E_max_global,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                worker.M_min_global = 0;
                worker.M_max_global = 0;
                break;
            case 2:
                MPI_Allreduce(MPI_IN_PLACE, &worker.E_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &worker.E_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &worker.M_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &worker.M_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                break;
            default:
                cout << "Error in check_windows(). Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
                MPI_Finalize();
                exit(1);
        }
    }

    void synchronize_sets(class_worker &worker){
        //This function collects all known energies from all sets
        vector<double> E_vector(worker.E_set.begin(),worker.E_set.end());
        vector<double> M_vector(worker.M_set.begin(),worker.M_set.end());
        vector<int> E_sizes((unsigned long)worker.world_size);
        vector<int> M_sizes((unsigned long)worker.world_size);
        int Esize = (int) E_vector.size();
        int Msize = (int) M_vector.size();
        MPI_Allgather(&Esize,1,MPI_INT, E_sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
        MPI_Allgather(&Msize,1,MPI_INT, M_sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
        vector<double> E_recv;
        vector<double> M_recv;
        for (int w = 0; w < worker.world_size; w++){
            if (w == worker.world_ID){
                E_recv = E_vector;
                M_recv = M_vector;
            }else{
                E_recv.resize((unsigned long) E_sizes[w]);
                M_recv.resize((unsigned long) M_sizes[w]);
            }
            MPI_Bcast(E_recv.data(),(int)E_sizes[w], MPI_DOUBLE,w,MPI_COMM_WORLD);
            MPI_Bcast(M_recv.data(),(int)M_sizes[w], MPI_DOUBLE,w,MPI_COMM_WORLD);
            worker.E_set.insert(E_recv.begin(),E_recv.end());
            worker.M_set.insert(M_recv.begin(),M_recv.end());
        }

    }

    void adjust_local_bins(class_worker &worker){
        // This function does rebinning of dos and histograms.
        //Snap the local boundaries to existing energies
        ArrayXd E_set_to_array;
        ArrayXd M_set_to_array;
        E_set_to_array << worker.E_set;
        M_set_to_array << worker.M_set;

        int E_idx_min = math::binary_search_nearest(E_set_to_array, worker.E_min_local);
        int E_idx_max = math::binary_search_nearest(E_set_to_array, worker.E_max_local);
        worker.E_min_local   = E_set_to_array(E_idx_min);
        worker.E_max_local   = E_set_to_array(E_idx_max);
        int M_idx_min = math::binary_search_nearest(M_set_to_array, worker.M_min_local);
        int M_idx_max = math::binary_search_nearest(M_set_to_array, worker.M_max_local);
        worker.M_min_local   = M_set_to_array(M_idx_min);
        worker.M_max_local   = M_set_to_array(M_idx_max);

        worker.E_bins                  = E_set_to_array.segment(E_idx_min, E_idx_max-E_idx_min+1);
        worker.M_bins                  = M_set_to_array.segment(M_idx_min, M_idx_max-M_idx_min+1);
        worker.histogram               = ArrayXXi::Zero(worker.E_bins.size(), worker.M_bins.size());
        worker.dos                     = ArrayXXd::Zero(worker.E_bins.size(), worker.M_bins.size());
    }

    void sync_team (class_worker &worker){
        //Check if anybody was able to do a walk
        if (worker.team.active) {
            timer::sync_team = 0;
            worker.t_sync_team.tic();
            int size = (int) worker.random_walk.size();
            ArrayXi sizes(worker.team.team_size);
            ArrayXi displ(worker.team.team_size);
            MPI_Allgather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, worker.team.MPI_COMM_TEAM);
            displ(0) = 0;
            for (int i = 1; i < worker.team.team_size; i++) {
                displ(i) = displ(i - 1) + sizes(i - 1);
            }

            vector<state> random_walk_recv((unsigned long) (sizes.sum()));
            MPI_Allgatherv(worker.random_walk.data(), size, MPI_2INT,  random_walk_recv.data(), sizes.data(),displ.data(), MPI_2INT, worker.team.MPI_COMM_TEAM);

            for (int i = 0; i < random_walk_recv.size(); i++) {
                worker.histogram(random_walk_recv[i].E_idx, random_walk_recv[i].M_idx) += 1;
                worker.dos(random_walk_recv[i].E_idx, random_walk_recv[i].M_idx) += worker.lnf;
            }
            counter::MCS            += constants::rate_sync_team * (worker.team.team_size - 1);
            timer::add_hist_volume  += constants::rate_sync_team * (worker.team.team_size - 1);
            timer::check_saturation += constants::rate_sync_team * (worker.team.team_size - 1);
            worker.random_walk.clear();
            worker.t_sync_team.toc();
        }
    }

    void setup_team(class_worker &worker){
        timer::setup_team = 0;
        worker.t_setup_team.tic();
        if (debug_setup_help) { debug_print(worker, "Reorganizing teams. "); }
        //Get a list of workers that have finished
        ArrayXi whos_finished(worker.world_size);
        MPI_Allgather(&worker.finish_line, 1, MPI_INT, whos_finished.data(), 1, MPI_INT, MPI_COMM_WORLD);


        int backfill_pos = 0;
        int team_id      = worker.world_ID / constants::team_size;
        //Note that team_id might be > constants::num_teams. In that case, simply fill them in where they're useful.
        ArrayXi team_filling = ArrayXi::Zero(constants::num_teams);
        for (int w = 0; w < worker.world_size; w++) {
            if (w == worker.world_ID){
                if (!worker.finish_line and team_id < constants::num_teams) {
                    team_filling(team_id)++;
                }
                else{
                    int iter = 0;
                    while(iter < worker.world_size*worker.world_size){
                        //Advance backfill_pos until you find a slot. Wrap around if you finish the list.
                        if(whos_finished(backfill_pos) == 1){
                            iter++;
                            backfill_pos = math::mod(backfill_pos+1, constants::num_teams);
                        }
                        else{break;}
                    }
                    //Place this guy where he is needed
                    assert(backfill_pos < constants::num_teams and backfill_pos >= 0);
                    team_id = backfill_pos;
                    team_filling(backfill_pos)++;
                    backfill_pos++;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(team_filling.data(),constants::num_teams, MPI_INT, w, MPI_COMM_WORLD);
            MPI_Bcast(&backfill_pos,1, MPI_INT, w, MPI_COMM_WORLD);
        }
        if (debug_setup_help) { debug_print(worker, " Backfilling successful. "); }

        assert(team_filling.sum() == worker.world_size);
        ArrayXi old_team_config(worker.world_size), new_team_config(worker.world_size);
        MPI_Allgather(&worker.team.team_id, 1, MPI_INT, old_team_config.data(), 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&team_id            , 1, MPI_INT, new_team_config.data(), 1, MPI_INT, MPI_COMM_WORLD);
        //If there has been no change since last time, simply return;
        if (new_team_config.cwiseEqual(old_team_config).all()) {
            if (debug_setup_help) { debug_print(worker, "No change in teams.\n"); }

            worker.t_setup_team.toc();
            return;
        }

        //Create Team communicators
        worker.team.team_id = team_id;
        setup_comm(worker);
//        if(worker.world_ID == 0) {cout << "Finished workers: " << whos_finished.transpose() << endl;}
//        if(worker.world_ID == 0) {cout << "Teams (" << new_team_size << ") : " << team_filling.transpose() << endl<<endl;}
        //Now every available guy knows who to help (rank 0 in MPI_COMM_TEAM) and the helpees know who to send their info to (Bcast).
        //Share details with the new team to begin.
        mpi::bcast_dynamic(worker.dos,       MPI_DOUBLE, 0, worker.team.MPI_COMM_TEAM);
        mpi::bcast_dynamic(worker.histogram, MPI_INT   , 0, worker.team.MPI_COMM_TEAM);
        mpi::bcast_dynamic(worker.E_bins,    MPI_DOUBLE, 0, worker.team.MPI_COMM_TEAM);
        mpi::bcast_dynamic(worker.M_bins,    MPI_DOUBLE, 0, worker.team.MPI_COMM_TEAM);
        MPI_Bcast(&worker.lnf,     1,        MPI_DOUBLE, 0, worker.team.MPI_COMM_TEAM);
        MPI_Bcast(&counter::MCS,   1,        MPI_INT   , 0, worker.team.MPI_COMM_TEAM);
        MPI_Bcast(&counter::walks, 1,        MPI_INT   , 0, worker.team.MPI_COMM_TEAM);
        ArrayXi saturation_map = Map<ArrayXi>(worker.saturation.data(), (int)worker.saturation.size());
        mpi::bcast_dynamic(saturation_map,   MPI_INT   , 0, worker.team.MPI_COMM_TEAM);
        if (!worker.team.team_leader) {
            worker.E_min_local      = worker.E_bins.minCoeff();
            worker.E_max_local      = worker.E_bins.maxCoeff();
            worker.M_min_local      = worker.M_bins.minCoeff();
            worker.M_max_local      = worker.M_bins.maxCoeff();
        }
        worker.set_rate_increment();
        worker.saturation.clear();
        worker.random_walk.clear();
        worker.state_is_valid       = false;
        worker.flag_one_over_t      = worker.lnf < 1.0 / max(1, counter::MCS) ? 1 : 0;
        worker.t_setup_team.toc();
        if (debug_setup_help) { debug_print(worker, "Success\n"); }

    }

    void setup_comm(class_worker &worker){
        //Create Team communicators
        int key, color_team, color_lead;
        worker.team.active  = true;
        if (!worker.finish_line && math::mod(worker.world_ID, constants::team_size) == 0) {
            color_team = worker.team.team_id;
            key = 0;
        } else {
            color_team = worker.team.team_id;
            key        = worker.team.team_id + 1;
        }
        MPI_Comm_split(MPI_COMM_WORLD, color_team, key, &worker.team.MPI_COMM_TEAM);
        MPI_Comm_rank(worker.team.MPI_COMM_TEAM, &worker.team.team_rank);
        MPI_Comm_size(worker.team.MPI_COMM_TEAM, &worker.team.team_size);

        worker.team.team_leader = worker.team.team_rank == 0;
        color_lead              = worker.team.team_leader ? 1 : MPI_UNDEFINED;
        MPI_Comm_split(MPI_COMM_WORLD, color_lead, key, &worker.team.MPI_COMM_LEAD);
        int team_commander = -1;
        if(worker.team.team_leader){
            MPI_Comm_rank(worker.team.MPI_COMM_LEAD, &team_commander);
        }
        worker.team.team_commander = team_commander == 0;
    }

}
