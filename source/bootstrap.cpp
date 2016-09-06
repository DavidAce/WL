//
// Created by david on 8/31/16.
//

#include "bootstrap.h"

void do_bootstrap(class_worker &worker){
    //Each worker loads its segment from random iteration
    //Then merge, then write out new bootstrapped DOS.
    outdata out(worker.world_ID);
    indata  in (worker.world_ID, worker.world_size);
    cout << "Bootstrap... " << endl;
    for (int i = 0; i < constants::bootstrap_reps; i++){
        worker.iteration = i + constants::simulation_reps;
        in.load_random_section(worker);
        mpi::merge    (worker) ;
        mpi::broadcast(worker) ;
        out.create_and_set_folder(worker.iteration);
        out.write_data_master(worker) ;
    }

    //Now we want to compute thermodynamic quantities
    do_thermodynamics(worker);
    //And get a good error estimation
    do_statistics(worker);

}

void do_thermodynamics(class_worker &worker){
    outdata out(worker.world_ID);
    indata  in (worker.world_ID, worker.world_size);
    class_thermodynamics thermo;
    worker.iteration = worker.world_ID;
    while (worker.iteration < (constants::bootstrap_reps + constants::simulation_reps)){
        in.load_full(worker);
        thermo.compute(worker);
        thermo.get_peak(worker);
        out.write_data_thermo(thermo, worker.iteration);
        worker.iteration += worker.world_size ;
    }

}

void do_statistics(class_worker &worker){
    if (worker.world_ID == 0 && debug_stats) {
        cout << "Initializing Stats" << endl;
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(1000));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    class_stats stats(worker.world_ID, worker.world_size);
    if (worker.world_ID == 0 && debug_stats) {
        cout << "Loading thermodynamic files" << endl;
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(10000));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    stats.load_thermo_data(worker);
    MPI_Barrier(MPI_COMM_WORLD);

    if (worker.world_ID == 0 && debug_stats) {
        cout << "Computing statistics" << endl;
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(10000));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    stats.compute(worker);
    MPI_Barrier(MPI_COMM_WORLD);

    if (worker.world_ID == 0 && debug_stats) {
        cout << "loading output files" << endl;
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(10000));
    }
    MPI_Barrier(MPI_COMM_WORLD);

    outdata out(worker.world_ID);
    if (worker.world_ID == 0 && debug_stats) {
        cout << "Writing final" << endl;
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(10000));
    }
    MPI_Barrier(MPI_COMM_WORLD);

    out.write_final_data(stats);
}
