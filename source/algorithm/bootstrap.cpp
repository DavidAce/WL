//
// Created by david on 8/31/16.
//

#include "bootstrap.h"
#define debug_boot      1
#define debug_thermo    0
#define debug_stats     0
void do_bootstrap(class_worker &worker){
    //Each worker loads its segment from random iteration
    //Then merge, then write out new bootstrapped DOS.
    outdata out;
    indata  in;
    if(debug_boot){parallel::debug_print(worker,"Bootstrap: Starting\n");}
    worker.team.team_id = worker.world_ID /constants::team_size;
    parallel::setup_comm(worker);
    if (worker.team.team_leader) {
        for (int i = 0; i < constants::bootstrap_reps; i++) {
            worker.iteration = i + constants::simulation_reps;
            in.load_random_section(worker);
            parallel::merge(worker, false, true);
            out.create_iteration_folder_master(worker.iteration, worker.world_ID);
            out.write_data_master(worker);
        }
    }
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
    if(debug_thermo){parallel::debug_print_team_commander(worker," Computing thermodynamic quantities...\n");}
    while (worker.iteration < (constants::bootstrap_reps + constants::simulation_reps)){
        if(debug_thermo){cout << "ID: " << worker.world_ID << "  Thermo: Loading data..." << endl;}
        in.load_full(worker);
        if(debug_thermo){cout << "ID: " << worker.world_ID << "  Thermo: Computing avgs..." << endl;}
        thermo.compute(worker);
        if(debug_thermo){cout << "ID: " << worker.world_ID << "  Thermo: Computing Peaks..." << endl;}
        thermo.get_c_peak(worker);
        thermo.get_x_peak(worker);
        thermo.get_Tc_free_energy(worker);
        thermo.get_Tc_canonical_distribution(worker);
        if(debug_thermo){cout << "ID: " << worker.world_ID << "  Thermo: Writing Results..." << endl;}
        out.write_data_thermo(thermo, worker.iteration);
        worker.iteration += worker.world_size ;
    }
}

void do_statistics(class_worker &worker){
    //This is a single threaded operation
    if(worker.world_ID == 0) {
        class_stats stats(worker.world_ID, worker.world_size);
        if(debug_stats){cout << "ID: " << worker.world_ID << " Stats: Loading thermodynamic files" << endl;}
        stats.load_thermo_data(worker);
        if(debug_stats){cout << "ID: " << worker.world_ID << " Stats: Computing statistics" << endl;}
        stats.compute(worker);
        if(debug_stats){cout << "ID: " << worker.world_ID << " Stats: Writing final results" << endl;}
        outdata out;
        out.write_final_data(stats, worker.world_ID);
    }
}
