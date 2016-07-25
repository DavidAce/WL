//
// Created by david on 2016-07-24.
//
#include <math.h>
#include <fstream>
#include "class_worker.h"
#include "counters.h"
using namespace std;
using namespace Eigen;

//Constructor
class_worker::class_worker(): model(), finish_line(0){
    MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);           //Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);         //Get total number of threads
    rn::rng.seed((unsigned long)world_ID);
    //Initialize matrices
    switch(constants::rw_dims){
        case 1:
            histogram.conservativeResizeLike(MatrixXi::Zero(constants::bins,1));
            dos.conservativeResizeLike(MatrixXd::Zero(constants::bins,1));

            break;
        case 2:
            histogram.conservativeResizeLike(MatrixXi::Zero(constants::bins,constants::bins));
            dos.conservativeResizeLike(MatrixXd::Zero(constants::bins,constants::bins));
            break;
        default:
            std::cout << "Error in class_worker constructor. Wrong dimension for WL-random walk (rw_dims = ?) " << std::endl;
    }
    dos_temp = dos;
    histogram_temp = histogram;
    lnf = 1;
    initial_limits();
    measure_state();
    initial_spectrum();
    reset_timers();
}


void class_worker::measure_state(){
    E = model.get_E();
    //Only give value to M if we are doing 2D random walks
    switch(constants::rw_dims){
        case 1:
            M = 0;
            break;
        case 2:
            M = model.get_M();
            break;
        default:
            cout << "Error in class_worker::measure_state(). Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            exit(1);
    }
}

void class_worker::initial_limits(){
    //Measure, randomize and measure again to get 2 states
    double E1 = model.get_E();
    double M1 = model.get_M();
    model.randomize_lattice();
    double E2 = model.get_E();
    double M2 = model.get_M();
    Emin_local = fmin(E1,E2);
    Emax_local = fmax(E1,E2);
    Mmin_local = fmin(M1,M2);
    Mmax_local = fmax(M1,M2);
    //Now compare to the other workers to find global limits
    switch(constants::rw_dims){
        case 1:
            MPI_Allreduce(&Emin_local, &Emin_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
            MPI_Allreduce(&Emax_local, &Emax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
            Mmin_global = 0;
            Mmax_global = 0;
            break;
        case 2:
            MPI_Allreduce(&Emin_local, &Emin_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
            MPI_Allreduce(&Emax_local, &Emax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
            MPI_Allreduce(&Mmin_local, &Mmin_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
            MPI_Allreduce(&Mmax_local, &Mmax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
            break;
        default:
            cout << "Error in class_worker::initial_state(). Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            exit(1);

    }

}

void class_worker::initial_spectrum(){
    switch(constants::rw_dims){
        case 1:
            E_bins = VectorXd::LinSpaced(constants::bins, Emin_local, Emax_local);
            M_bins = VectorXd::Zero(constants::bins);
            Eidx = binary_search(E_bins.data(), E, E_bins.size());
            Midx = 0;
            break;
        case 2:
            E_bins = VectorXd::LinSpaced(constants::bins, Emin_local, Emax_local);
            M_bins = VectorXd::LinSpaced(constants::bins, Mmin_local, Mmax_local);
            Eidx = binary_search(E_bins.data(), E, E_bins.size());
            Midx = binary_search(M_bins.data(), M, M_bins.size());
            break;
        default:
            cout << "Error in class_worker::initial_spectrum. Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            exit(1);
    }
    X_bins = E_bins;

    cout << Midx << endl;

}

void class_worker::update_limits() {
    int i,j,k;

//    cout << "Updating Limits. E_size =  " <<  E_bins.size() << " dos_size =  " <<  dos.size() <<endl;
    if(E_trial < Emin_local || E_trial > Emax_local ) {
        Emin_local = fmin(E_trial, Emin_local);
        Emax_local = fmax(E_trial, Emax_local);
        X_bins = E_bins;
        dos_temp.fill(0);
        histogram_temp = histogram;
        E_bins = VectorXd::LinSpaced(constants::bins, Emin_local, Emax_local);
        j = 0;
        i = 0;
        k = 0;
        double dE = fabs(E_bins(1) - E_bins(0));
        for (i = 0; i < E_bins.size() ; i++){
            k = 0;
            for(j = 0; j < E_bins.size(); j++){
                if(fabs(E_bins(i) - X_bins(j)) >= dE) {
                    continue;
                }else {
                    dos_temp.row(i) += dos.row(j);
                    k++;
                }
            }
            if (k > 0){
                dos_temp.row(i) /= (double)k;
            }


        }
        dos = dos_temp;


//
//
//        while (i < E_bins.size() && j < E_bins.size()){
//            if (E_bins(i) < X_bins(j)){
//                //Take care of border
//                if(j == 0){
//                    dos.row(i).fill(0);
//                    histogram.row(i).fill(0);
//                    i++;
//                    continue;
//                }else{
//                    i++;
//                }
//            }else{
//                if (j == E_bins.size()-1){
//                    dos.row(i).fill(0);
//                    i++;
//                    continue;
//                }
//                k = 0;
//                while(E_bins(i) >= X_bins(j) && i < E_bins.size() && j < E_bins.size() ){
//                    j++;
//                    k++;
//                }
//                dos.row(i) = dos_temp.middleRows(j-k,k).colwise().mean();
//                i++;
//            }
//        }
//        for (int i = 0; i < E_bins.size() && j < E_bins.size(); i++) {
//            if (E_bins(i) < X_bins(j)) {
////                cout << "ID "<< world_ID << " Ebin(" << i << ") = " << E_bins(i);
////                cout << " Xbin(" << j << ") = " << X_bins(j) << endl;
//                if (world_ID == 0) {
//                    cout << "RESIZE BACK!!!" << endl << endl;
//                    cout << "DOS: " << endl;
//                    cout << dos << endl;
//                }
//                dos.row(i).fill(0);
//                histogram.row(i).fill(0);
//            } else {
//                while (E_bins(i) >= X_bins(j) && j < E_bins.size()) {
//                    j++;
//                }
//               // if (world_ID == 0) {
//                    cout << "j = " << j << endl;// BACK!!!" << endl << endl;
//
//                //}
//                dos.row(i) = dos.middleRows(i,j-i)/(double)(j-i);
//                histogram.row(i) = histogram.middleRows(i,j-i)/(j-i);
//
//            }
//
////        cout << "Finished loop" << endl;
////      cout << Eidx << endl;
//
//        }
        Eidx = binary_search(E_bins.data(), E, E_bins.size());

    }


//    cout << "ID "<< world_ID << " No resize. E: " << E_bins.transpose() << endl;


    if((M_trial < Mmin_local || M_trial > Mmax_local) && constants::rw_dims == 2){
        Mmin_local = fmin(M_trial, Mmin_local);
        Mmax_local = fmax(M_trial, Mmax_local);
        j = 0;
        X_bins = M_bins;
        M_bins = VectorXd::LinSpaced(constants::bins, Mmin_local, Mmax_local);
        for(int i = 0; i < M_bins.size() && j < M_bins.size() ; i++){
            if(M_bins(i) < X_bins(j)){
                dos.col(i).fill(0);
                histogram.col(i).fill(0);
            }else{
                while(M_bins(i) > X_bins(j) && j < M_bins.size()){
                    dos.col(i) += dos.col(j);
                    histogram.col(i) += histogram.col(j);
                    j++;
                }
            }
        }
        Midx = binary_search(M_bins.data(), M, M_bins.size());
    }
}



void class_worker::reset_timers() {
    MCS = 0;
    timer_append = 0;
    timer_check = 0;
    timer_print = 0;
    timer_swap = 0;
    timer_resize = 0;
}

void class_worker::make_MC_trial(){
    model.make_new_state(E,M, E_trial, M_trial);
    update_limits();


    switch (constants::rw_dims){
        case 1:
            Eidx = binary_search(E_bins.data(), E, E_bins.size());
            Midx = 0;
            Eidx_trial = binary_search(E_bins.data(), E_trial, E_bins.size());
            Midx_trial = 0;
            break;
        case 2:
            Eidx = binary_search(E_bins.data(), E, E_bins.size());
            Midx = binary_search(M_bins.data(), M, M_bins.size());

            Eidx_trial = binary_search(E_bins.data(), E_trial, E_bins.size());
            Midx_trial = binary_search(M_bins.data(), M_trial, M_bins.size());
            break;
        default:
            cout << "Error in class_worker::make_MC_trial(). Wrong dimension for WL-random walk (rw_dims = ?)" << endl;
            exit(1);
    }
//    if (world_ID == 0){
//        cout << "E: " << E << " Eidx: " << Eidx << " Enew: " << E_trial << " Eidx_trial: " << Eidx_trial << endl;
//        cout << E_bins.transpose() << endl;
//        cout << dos.transpose() << endl << endl;
//
//    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void class_worker::accept_MC_trial(){
    E = E_trial;
    M = M_trial;
    Eidx = Eidx_trial;
    Midx = Midx_trial;
    model.flip();
    dos(Eidx,Midx) += lnf;
    histogram(Eidx,Midx) += 1;


}

void class_worker::reject_MC_trial(){
    dos(Eidx,Midx) += lnf;
    histogram(Eidx,Midx) += 1;
}



//Function for printing the lattice. Very easy if L is an Eigen matrix.
std::ostream &operator<<(std::ostream &os, const class_worker &worker){
    os << "Worker(" << worker.world_ID << ") ["
       << worker.Emin_local << " "
       << worker.Emax_local << "] ["
       << worker.Mmin_local << " "
       << worker.Mmax_local << "] ["
       << worker.Emin_global << " "
       << worker.Emax_global << "] ["
       << worker.Mmin_global << " "
       << worker.Mmax_global << "]"
       << std::endl;
    return os;
}


