//
// Created by david on 2016-07-24.
//
#include <iomanip>
#include "constants.h"
#include "counters_timers.h"
#include "class_data.h"
#include "wanglandau.h"
using namespace std;

void WangLandau(class_worker &worker){
    int finish_line = 0;
    outdata out(worker.world_ID);
    while(finish_line == 0){
        sweep(worker);
        check_convergence(worker, finish_line);
        worker.split_global_spectrum();
        //swap(worker);
        print_status(worker);
        //backup_data(worker);
    }
    out.write_data(worker);
}

void sweep(class_worker &worker){
    for (int i = 0; i < constants::N ; i++){
        worker.make_MC_trial();
        worker.acceptance_criterion();
        if(worker.accept){
            worker.accept_MC_trial();

        }else{
            worker.reject_MC_trial();
        }
    }
    counter::MCS++;
}

void check_convergence(class_worker &worker, int &finish_line){
    switch(worker.flag_one_over_t){
        case 0:
            add_hist_volume(worker);
            check_saturation(worker);
            break;
        case 1:
            check_one_over_t(worker);
            break;
        default:
            cout << "Error: check_convergence has wrong flag" << endl;
            exit(2);
    }
    if (timer::check_finish_line > constants::rate_check_finish_line) {
        timer::check_finish_line = 0;
        MPI_Allreduce(&worker.finish_line, &finish_line, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    }else{
        timer::check_finish_line++;
    }
}


void add_hist_volume(class_worker &worker) {
    if (timer::add_hist_volume >= constants::rate_add_hist_volume) {
        timer::add_hist_volume = 0;
        //Find minimum histogram entry larger than 2.
        int min = math::find_min_positive(worker.histogram);
        for (int j = 0; j < worker.histogram.cols(); j++) {
            for (int i = 0; i < worker.histogram.rows(); i++) {
                if (worker.histogram(i, j) == 0) { continue; }
                worker.histogram(i, j) -= min;
            }
        }
        //Resize saturation vector
        if (worker.saturation.size() <= counter::saturation) {
            worker.saturation.conservativeResize(2 * counter::saturation);
            worker.saturation.bottomRows(counter::saturation).fill(0);
        }
        worker.saturation(counter::saturation) = worker.histogram.sum();
        counter::saturation++;
    }else{
        timer::add_hist_volume++;
    }
}

void check_saturation(class_worker &worker) {
    if (timer::check_saturation >= constants::rate_check_saturation) {
        timer::check_saturation = 0;
        int i, j;
        int idx = (int) (constants::check_saturation_from * counter::saturation);
        double Sx = 0, Sxy = 0, mX = 0, mY = 0;
        j = 0;
        //Compute means of the last 10%:
        for (i = idx; i < counter::saturation; i++) {
            mX += i;//X(i);
            mY += worker.saturation(i);
            j++;
        }
        mX /= j;
        mY /= j;
        for (i = idx; i < counter::saturation; i++) {
            Sx += pow(i - mX, 2);
            Sxy += (worker.saturation(i) - mY) * (i - mX);
        }
        double slope = Sxy / Sx;

        if (slope < 0) {
            worker.next_WL_iteration();
            if (worker.lnf < constants::minimum_lnf){
                worker.finish_line = 1;
            }
            if(worker.lnf < 1/counter::MCS){
            //if(counter::MCS > 1e5){
                worker.lnf = 1.0/counter::MCS;
                worker.flag_one_over_t = 1;         //Change to 1/t algorithm
            }
        }
    }else{
        timer::check_saturation++;
    }
}


void merge_windows (class_worker &worker) {

    //Get id of neighbours up and down.
    int dn = worker.world_ID > 0 ? worker.world_ID - 1 : worker.world_size - 1;
    int up = worker.world_ID < worker.world_size - 1 ? worker.world_ID + 1 : 0;
    int E_size_up, E_size = (int)worker.E_bins.size();
    int M_size_up, M_size = (int)worker.M_bins.size();

    //Find out the size of the data of the neighbor above
    MPI_Sendrecv(&E_size, 1, MPI_INT, dn, 0, &E_size_up, 1, MPI_INT, up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&M_size, 1, MPI_INT, dn, 0, &M_size_up, 1, MPI_INT, up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MatrixXd dos_up(E_size_up, M_size_up);
    VectorXd E_bins_up(E_size_up);
    VectorXd M_bins_up(M_size_up);
    //Receive it as well
    MPI_Sendrecv(worker.dos.data(), (int) worker.dos.size(), MPI_DOUBLE, dn, 0, dos_up.data(), E_size_up*M_size_up, MPI_DOUBLE, up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(worker.E_bins.data(), (int) worker.E_bins.size(), MPI_DOUBLE, dn, 0, E_bins_up.data(), E_size_up, MPI_DOUBLE, up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(worker.M_bins.data(), (int) worker.M_bins.size(), MPI_DOUBLE, dn, 0, M_bins_up.data(), M_size_up, MPI_DOUBLE, up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (worker.E_bins.maxCoeff() < E_bins_up.minCoeff() && up != 0 ){
        //cout << "Error in merge_windows(). Windows do not overlap." << endl;
        //exit(3);
    }
    if (worker.E_bins.maxCoeff() > E_bins_up.maxCoeff() && up != 0 ){
        //cout << "Error in merge_windows(). Overlap too large." << endl;
        //exit(4);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    Vector3d v1A, v2A, v1B, v2B, uA, uB; //Vectors connecting adjacent 3 orthogonal points on DOS
    VectorXd sum(max(E_size,E_size_up));
    sum.fill(0);
    int col = 0, num;
    int x, y; //Coordinates closest to i, j on dos of neighbor above.
    for (int i = 0; i < E_size - 1; i++) {
        if (worker.E_bins(i) < E_bins_up.minCoeff()) { continue;}
        num = 0;
        for (int j = 0; j < M_size - 1; j++) {
            if (worker.M_bins(j) > M_bins_up.maxCoeff()) { continue;}
            if (worker.M_bins(j) < M_bins_up.minCoeff()) { continue;}
            if (worker.dos(i,j) == 0) { continue;}

            x = math::upper_bound(E_bins_up.data(), worker.E_bins(i), E_size_up);
            y = math::upper_bound(M_bins_up.data(), worker.M_bins(i), M_size_up);
            v1A << worker.E_bins(i + 1) - worker.E_bins(i)  , 0, worker.dos(i + 1, j) - worker.dos(i, j);
            v1B << E_bins_up(x + 1) - E_bins_up(x)          , 0, dos_up(x + 1, y) - dos_up(x,y);
            v2A << 0, worker.M_bins(j + 1) - worker.M_bins(j)   , worker.dos(i, j+1) - worker.dos(i, j);
            v2B << 0, M_bins_up(j + 1) - M_bins_up(j)           , dos_up(i, j+1) - dos_up(i, j);
            uA = v1A.cross(v2A).normalized();
            uB = v1B.cross(v2B).normalized();
            sum(col) += uA.dot(uB);
            num++;
        }

        sum(col) /= max(1,num);
        col++;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int E_merge_idx = math::find_max_idx(sum);
    int M_merge_idx = math::find_max_idx(worker.dos.row(E_merge_idx));

    //Find the merge coordinates on the upper dos
    int E_merge_idx_up = math::upper_bound(E_bins_up.data(), worker.E_bins(E_merge_idx), E_size_up);
    int M_merge_idx_up = math::upper_bound(M_bins_up.data(), worker.M_bins(M_merge_idx), M_size_up);

    double diff;


    for (int w = 0; w < worker.world_size - 1; w++) {
        if (worker.world_ID == w) {
            //Adjust the dos above and send up
            diff = worker.dos(E_merge_idx, M_merge_idx) - dos_up(E_merge_idx_up, M_merge_idx_up);
            for (int j = 0; j < M_size_up; j++) {
                for (int i = 0; i < E_size_up; i++) {
                    if (dos_up(i, j) == 0) { continue; }
                    else {
                        dos_up(i, j) += diff;
                        dos_up(i, j) = dos_up(i, j) <= 0 ? 0 : dos_up(i, j);
                    }
                }

            }
            MPI_Send(dos_up.data(),E_size_up*M_size_up, MPI_DOUBLE, up,0,MPI_COMM_WORLD);
        }
        else if ( worker.world_ID == w+1){
            MPI_Recv(worker.dos.data(),E_size*M_size, MPI_DOUBLE, dn,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //cout << worker.dos.transpose() << endl << endl;

        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}


void check_one_over_t (class_worker &worker){
    worker.lnf = 1.0/counter::MCS;
    if (worker.lnf < constants::minimum_lnf){
        worker.finish_line = 1;
    }

}

void print_status(class_worker &worker){
    if (timer::print >= constants::rate_print_status){
        timer::print = 0;

        for (int i = 0; i < worker.world_size; i++){
            if(worker.world_ID == i){
                if (worker.world_ID == 0){cout << endl;}

                cout    << "ID: "        << i
                        << "    Walk: "  << counter::walks
                        << "    f: "     << fixed << setprecision(12) << exp(worker.lnf)
                        << "    Bins: ["  << worker.dos.rows() << "  " << worker.dos.cols() << "]"
                        << "    1/t: "   << worker.flag_one_over_t
                        << "    Fin: "   << worker.finish_line
                        << endl << flush;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0){
            cout    << "-----"
                    << "MCS: "          << counter::MCS
                    <<  "    MaxWalks: "<< fixed << setprecision(0) << ceil(log(constants::minimum_lnf)/log(constants::reduce_factor_lnf))
                    << "  -----"
                    <<endl << flush;
        }
    }else{
        timer::print++;
    }

}
