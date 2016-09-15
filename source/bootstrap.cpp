//
// Created by david on 8/31/16.
//

#include "bootstrap.h"
#define debug_boot      1
#define debug_thermo    1
#define debug_stats     1
void do_bootstrap(class_worker &worker){
    //Each worker loads its segment from random iteration
    //Then merge, then write out new bootstrapped DOS.
    outdata out;
    indata  in;
    if (worker.world_ID == 0 && debug_boot) {
        cout << "Bootstrap: Starting" << endl;
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(1000));
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < constants::bootstrap_reps; i++){
        worker.iteration = i + constants::simulation_reps;
        in.load_random_section(worker);
        mpi::merge    (worker) ;
        out.create_iteration_folder_master(worker.iteration, worker.world_ID);
        out.write_data_master(worker) ;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Now we want to compute thermodynamic quantities
    do_thermodynamics(worker);
    //And get a good error estimation
    do_statistics(worker);
    MPI_Barrier(MPI_COMM_WORLD);

}

void do_thermodynamics(class_worker &worker){
    outdata out;
    indata  in;
    class_thermodynamics thermo;
    worker.iteration = worker.world_ID;
    while (worker.iteration < (constants::bootstrap_reps + constants::simulation_reps)){
        if (debug_thermo) {
            cout << "ID: "<< worker.world_ID << " Thermo: Loading..." << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        in.load_full(worker);
        if (debug_thermo) {
            cout << "ID: "<< worker.world_ID << " Thermo: Computing avgs..." << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        thermo.compute(worker);
        if (debug_thermo) {
            cout << "ID: "<< worker.world_ID << " Thermo: Computing Peaks..." << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        thermo.get_c_peak(worker);
        thermo.get_x_peak(worker);
        thermo.get_Tc_free_energy(worker);

        if (debug_thermo) {
            cout << "ID: "<< worker.world_ID << " Thermo: Writing Data..." << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        out.write_data_thermo(thermo, worker.iteration);
        worker.iteration += worker.world_size ;
    }
    MPI_Barrier(MPI_COMM_WORLD);

}

void do_statistics(class_worker &worker){
    //This is a single threaded operation
    if(worker.world_ID == 0) {
        class_stats stats(worker.world_ID, worker.world_size);
        if (debug_stats) {
            cout << "Stats: Loading thermodynamic files" << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        stats.load_thermo_data(worker);
        if (debug_stats) {
            cout << "Stats: Computing statistics" << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        stats.compute(worker);
        if (debug_stats) {
            cout << "Writing final" << endl;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }
        outdata out;
        out.write_final_data(stats, worker.world_ID);
    }
}
