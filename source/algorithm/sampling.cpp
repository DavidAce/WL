//
// Created by david on 2018-05-20.
//

#include "sampling.h"
//
//void do_bootstrap(class_worker &worker){
//    //Each worker loads its segment from random iteration
//    //Then merge, then write out new bootstrapped DOS.
//    outdata out;
//    indata  in;
//    if(debug_boot){parallel::debug_print(worker,"Bootstrap: Starting\n");}
//    worker.team.team_id = worker.world_ID / constants::team_size;
//    parallel::setup_comm(worker);
//    if (worker.team.team_leader) {
//        for (int i = 0; i < constants::bootstrap_reps; i++) {
//            worker.iteration = i + constants::simulation_reps;
//            in.load_random_section(worker);
//            parallel::merge(worker, false, true);
//            out.create_iteration_folder_master(worker.iteration, worker.world_ID);
//            out.write_data_master(worker);
//        }
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    //Now compute thermodynamic quantities
//    do_thermodynamics(worker);
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    //And get error estimation
//    do_statistics(worker);
//    MPI_Barrier(MPI_COMM_WORLD);
//}