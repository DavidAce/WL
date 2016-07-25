//
// Created by david on 2016-07-24.
//
#include "constants.h"
#include "counters.h"
#include "randomNumbers.h"
#include "class_worker.h"
#include "wanglandau.h"
using namespace std;

void WangLandau(class_worker &worker){

    while(worker.finish_line == 0){
        sweep(worker);
        //update_limits();
        //check_convergence();
        //swap(worker);
        //print_status();
//        std::cout << "ID " << worker.world_ID << " MCS: " << worker.MCS << std::endl;
        if (worker.MCS > 30000){
            worker.finish_line = 1;
        }
    }
    //write_to_file();
    //
}

void sweep(class_worker &worker){

    for (int i = 0; i < constants::N ; i++){
        worker.make_MC_trial();
//        if ( worker.Eidx < 0 || worker.Eidx >= constants::bins){
//            cout << "Eidx out of bounds: " << worker.Eidx <<  endl;
//        }
//        if ( worker.Eidx_trial < 0 || worker.Eidx_trial >= constants::bins){
//            cout << "Eidx_trial out of bounds: " << worker.Eidx_trial <<  endl;
//        }
//        if ( worker.Midx < 0 || worker.Midx >= constants::bins){
//            cout << "Midx out of bounds: " << worker.Midx <<  endl;
//        }
//        if ( worker.Midx_trial < 0 || worker.Midx_trial >= constants::bins){
//            cout << "Midx_trial out of bounds: " << worker.Midx_trial <<  endl;
//        }

//        generate_new_energy();
//        check_limits();
////
        if(rn::uniform_double(0,1) < exp(worker.dos(worker.Eidx, worker.Midx) - worker.dos(worker.Eidx_trial, worker.Midx_trial))){
            worker.accept_MC_trial();

        }else{
            worker.reject_MC_trial();
        }
    }
    worker.MCS++;

}

