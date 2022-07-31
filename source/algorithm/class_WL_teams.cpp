//
// Created by david on 2018-05-26.
//

#include "class_WL_teams.h"
#include <algorithm/class_WL_worker.h>
#include <general/nmspc_mpi_extensions.h>

#define debug_sync_teams 0
#define debug_setup_teams 0
#define debug_setup_comms 0

class_WL_teams::class_WL_teams(class_worker &worker_ref, int num_teams_) : worker(worker_ref) {
    num_teams = num_teams_;
    set_defaults();
    setup_comms();
}

void class_WL_teams::set_defaults() {
    team_id = static_cast<int>((static_cast<double>(num_teams) / worker.world_size * worker.world_ID));
    setup_comms();
}

void class_WL_teams::sync_teams() {
    //    worker.debug_print_all<1>("size: " + to_string(worker.random_walk.size())+ "\n");
    auto size = static_cast<int>(worker.random_walk.size());
    //    worker.debug_print_all<1>("size: " + to_string(size)+ "\n");
    //    debug_print_team_commander("\n \n");
    Eigen::ArrayXi sizes(team_size);
    Eigen::ArrayXi displ(team_size);
    MPI_Allgather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, MPI_COMM_TEAM);
    displ(0) = 0;
    for(int i = 1; i < team_size; i++) { displ(i) = displ(i - 1) + sizes(i - 1); }

    std::vector<state> random_walk_recv((unsigned long) (sizes.sum()));
    MPI_Allgatherv(worker.random_walk.data(), size, MPI_2INT, random_walk_recv.data(), sizes.data(), displ.data(), MPI_2INT, MPI_COMM_TEAM);

    for(auto &rwri : random_walk_recv) {
        worker.histogram(rwri.E_idx, rwri.M_idx) += 1;
        worker.dos(rwri.E_idx, rwri.M_idx) += worker.lnf;
    }
    worker.random_walk.clear();
}

void class_WL_teams::setup_teams() {
    // Get a list of workers that have finished
    Eigen::ArrayXi workers_finished(worker.world_size);
    MPI_Allgather(&worker.finish_line, 1, MPI_INT, workers_finished.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int team_id_old = team_id;
    int team_id_new = static_cast<int>((static_cast<double>(num_teams) / worker.world_size * worker.world_ID));

    // First get a list of forbidden teams
    Eigen::ArrayXi team_forbidden = Eigen::ArrayXi::Zero(num_teams);
    for(int w = 0; w < worker.world_size; w++) {
        if(w == worker.world_ID) {
            // If you are finished, mark your team as a forbidden team.
            if(worker.finish_line == 1) { team_forbidden(team_id_new) = 1; }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(team_forbidden.data(), num_teams, MPI_INT, w, MPI_COMM_WORLD);
    }

    auto           max_filling_per_slot = std::ceil(static_cast<double>(worker.world_size) / (1 - team_forbidden).sum());
    Eigen::ArrayXi team_filling         = Eigen::ArrayXi::Zero(num_teams);
    for(int w = 0; w < worker.world_size; w++) {
        if(w == worker.world_ID) {
            // If team_id is assigned to a team that's already finished, find another team that
            //  has slots left. Unless you are a teamleader, in which case you just stay put.
            if(not team_leader) {
                while(team_forbidden(team_id_new) == 1 or team_filling(team_id_new) >= max_filling_per_slot) {
                    team_id_new = math::mod(team_id_new + 1, num_teams);
                }
            }
            assert(team_id_new < num_teams);
            team_filling(team_id_new)++;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(team_filling.data(), num_teams, MPI_INT, w, MPI_COMM_WORLD);
    }

    //    if (debug_setup_teams) { debug_print(worker, " Team assignment successful. "); }
    assert(team_filling.sum() == worker.world_size);
    Eigen::ArrayXi old_team_config(worker.world_size), new_team_config(worker.world_size);
    MPI_Allgather(&team_id_old, 1, MPI_INT, old_team_config.data(), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&team_id_new, 1, MPI_INT, new_team_config.data(), 1, MPI_INT, MPI_COMM_WORLD);
    // If there has been no change since last time, simply return;
    if(new_team_config.cwiseEqual(old_team_config).all()) {
        //        if (debug_setup_teams) { debug_print(worker, "No change in teams.\n"); }

        return;
    }

    // Create Team communicators
    team_id = team_id_new;
    setup_comms();
    //        if(worker.world_ID == 0) {std::cout << "Finished workers: " << workers_finished.transpose() << std::endl;}
    //        if(worker.world_ID == 0) {std::cout << "Teams (" << new_team_size << ") : " << team_filling.transpose() << std::endl<<std::endl;}
    // Now every available guy knows who to help (rank 0 in MPI_COMM_TEAM) and the helpees know who to send their info to (Bcast).
    // Share details with the new team to begin.
    mpi::bcast_dynamic(worker.dos, MPI_DOUBLE, 0, MPI_COMM_TEAM);
    mpi::bcast_dynamic(worker.histogram, MPI_INT, 0, MPI_COMM_TEAM);
    mpi::bcast_dynamic(worker.E_bins, MPI_DOUBLE, 0, MPI_COMM_TEAM);
    mpi::bcast_dynamic(worker.M_bins, MPI_DOUBLE, 0, MPI_COMM_TEAM);
    MPI_Bcast(&worker.lnf, 1, MPI_DOUBLE, 0, MPI_COMM_TEAM);
    MPI_Bcast(&counter::MCS, 1, MPI_INT, 0, MPI_COMM_TEAM);
    MPI_Bcast(&counter::walks, 1, MPI_INT, 0, MPI_COMM_TEAM);

    Eigen::ArrayXi saturation_map = Eigen::Map<Eigen::ArrayXi>(worker.saturation.data(), static_cast<long>(worker.saturation.size()));
    mpi::bcast_dynamic(saturation_map, MPI_INT, 0, MPI_COMM_TEAM);
    if(!team_leader) {
        worker.E_min_local = worker.E_bins.minCoeff();
        worker.E_max_local = worker.E_bins.maxCoeff();
        worker.M_min_local = worker.M_bins.minCoeff();
        worker.M_max_local = worker.M_bins.maxCoeff();
    }
    worker.set_rate_increment();
    worker.saturation.clear();
    worker.random_walk.clear();
    worker.state_is_valid  = false;
    worker.flag_one_over_t = worker.lnf < 1.0 / std::max(1, counter::MCS) ? 1 : 0;
    worker.t_setup_team.toc();
    //    if (debug_setup_teams) { debug_print(worker, "Success\n"); }
}

void class_WL_teams::setup_comms() {
    // Create Team communicators
    int key, color_team, color_lead;
    active     = true;
    color_team = team_id;
    if(team_leader) {
        key = 0;
    } else {
        key = 1;
    }

    MPI_Comm_split(MPI_COMM_WORLD, color_team, key, &MPI_COMM_TEAM);
    MPI_Comm_rank(MPI_COMM_TEAM, &team_rank);
    MPI_Comm_size(MPI_COMM_TEAM, &team_size);
    team_leader = team_rank == 0;

    // Create communicator among leaders
    color_lead = team_leader ? 1 : MPI_UNDEFINED;
    key        = 0;
    MPI_Comm_split(MPI_COMM_WORLD, color_lead, key, &MPI_COMM_LEAD);

    if(team_leader) {
        MPI_Comm_rank(MPI_COMM_LEAD, &team_leader_rank);
        MPI_Comm_size(MPI_COMM_LEAD, &team_leader_size);
    }
    team_commander = team_leader_rank == 0;

    if(team_leader) {
        team_commander_world_id = worker.world_ID;
        MPI_Bcast(&team_commander_world_id, 1, MPI_INT, 0, MPI_COMM_LEAD);
    }
    MPI_Bcast(&team_commander_world_id, 1, MPI_INT, 0, MPI_COMM_TEAM);

    MPI_Barrier(MPI_COMM_WORLD);
    if(debug_setup_comms) {
        for(int w = 0; w < worker.world_size; w++) {
            if(w == worker.world_ID) {
                std::cout << "Team properties for ID: " << worker.world_ID << "\n"
                          << "Team id               : " << team_id << "\n"
                          << "Team rank             : " << team_rank << "\n"
                          << "Team size             : " << team_size << "\n"
                          << "Team finished         : " << team_finished << "\n"
                          << "Team commander        : " << team_commander << "\n"
                          << "Team leader           : " << team_leader << "\n"
                          << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}