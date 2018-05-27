//
// Created by david on 8/31/16.
//
#include <algorithm/class_WL_worker.h>
#include <algorithm/class_WL_teams.h>
#include <algorithm/nmspc_WL_parallelization.h>
#include <algorithm/class_WL_statistics.h>
#include <algorithm/class_WL_thermo.h>
#include <IO/class_WL_print_data.h>
#include <IO/class_WL_read_data.h>
#include <general/nmspc_math_algorithms.h>
#include "bootstrap.h"
#define debug_boot      0
#define debug_thermo    0
#define debug_stats     0


void do_bootstrap(class_worker &worker){
    //Each worker loads its segment from random iteration
    //Then merge, then write out new bootstrapped DOS.
    outdata out;
    indata  in;
    worker.debug_print<debug_boot>("Bootstrap: Starting\n");
    worker.team->set_defaults();
//    if (worker.team->is_leader()) {
    for (int i = 0; i < constants::bootstrap_reps; i++) {
        worker.iteration = i + constants::simulation_reps;
        in.load_random_section(worker);
        parallel::merge(worker, false, true);
        out.create_iteration_folder_commander(worker, worker.iteration);
        out.write_data_commander(worker);
    }
//    }
    MPI_Barrier(MPI_COMM_WORLD);
    //Now compute thermodynamic quantities
    do_thermodynamics(worker);
    MPI_Barrier(MPI_COMM_WORLD);

    //And get error estimation
    do_statistics(worker);
    MPI_Barrier(MPI_COMM_WORLD);
}

void do_thermodynamics(class_worker &worker){
    outdata out;
    indata  in;
    class_thermodynamics thermo;
    worker.iteration = worker.world_ID;
    worker.team->debug_print_team_commander<debug_thermo>(" Computing thermodynamic quantities...\n");
    while (worker.iteration < (constants::bootstrap_reps + constants::simulation_reps)){
        if(debug_thermo){cout << "ID " << worker.world_ID << ": Thermo: Loading data..." << endl;}
        in.load_full(worker);
        if(debug_thermo){cout << "ID " << worker.world_ID << ": Thermo: Computing avgs..." << endl;}
        thermo.compute(worker);
        if(debug_thermo){cout << "ID " << worker.world_ID << ": Thermo: Writing Results..." << endl;}
        out.write_data_thermo(thermo, worker.iteration);
        worker.iteration += worker.world_size ;
    }
}

void do_statistics(class_worker &worker){
    //This is a single threaded operation
    if(worker.team->is_commander()) {
        class_stats stats;
        if(debug_thermo){cout << "ID " << worker.world_ID << ": Stats: Loading thermodynamic files..." << endl;}
        stats.load_thermo_data(worker);
        if(debug_thermo){cout << "ID " << worker.world_ID << ": Stats: Computing statistics..." << endl;}
        stats.compute(worker);
        if(debug_thermo){cout << "ID " << worker.world_ID << ": Stats: Writing final results..." << endl;}
        outdata out;
        out.write_final_data(worker, stats);
    }
}
