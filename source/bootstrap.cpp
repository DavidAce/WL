//
// Created by david on 8/31/16.
//

#include "bootstrap.h"

void do_bootstrap(class_worker &worker){
    //Each worker loads its segment from random iteration
    //Then merge, then write out new bootstrapped DOS.
    outdata out(worker.world_ID);
    indata  in (worker.world_ID, worker.world_size);

    for (int i = 0; i < constants::bootstrap_reps; i++){
        worker.iteration = i + constants::simulation_reps;
        in.load_random_section(worker);
        mpi::merge    (worker) ;
        out.create_and_set_folder(worker.iteration);
        out.write_data_master(worker) ;
    }
    //Now we want to compute thermodynamic quantities
    do_thermodynamics(worker);
    do_statistics(worker);

}

void do_thermodynamics(class_worker &worker){
    outdata out(worker.world_ID);
    indata  in (worker.world_ID, worker.world_size);
    class_thermodynamics thermo;
    worker.iteration = worker.world_ID;
    while (worker.iteration < (constants::bootstrap_reps + constants::simulation_reps)){
        cout << "Loading iteration " << worker.iteration << endl;
        in.load_full(worker);
        thermo.compute(worker);
        out.write_data_thermo(thermo, worker.iteration);
        cout << "Wrote iteration " << worker.iteration << endl;

        worker.iteration += worker.world_size ;
    }
    MPI_Barrier(MPI_COMM_WORLD);

}

void do_statistics(class_worker &worker){
    class_stats stats(worker.world_ID, worker.world_size);
    stats.load_thermo_data();
    stats.compute();
    outdata out(worker.world_ID);
    out.write_final_data(stats);
}
