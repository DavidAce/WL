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
        split_spectrum(worker);
        //worker.resize_global_range();
        //swap(worker);
        print_status(worker);
        //backup_data(worker);
    }
    out.write_data_worker(worker);
    merge_windows(worker);
    out.write_data_master(worker);


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

void split_spectrum(class_worker &worker){
    if (timer::split_windows >= constants::rate_resize_global_range) {
        timer::split_windows = 0;
        MPI_Allreduce(&worker.need_to_resize, &worker.need_to_resize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (worker.need_to_resize == 1){

            worker.prev_WL_iteration();

            worker.resize_global_range();
            worker.divide_global_range();
            worker.resize_local_bins();
            worker.need_to_resize = 0;
        }
    }else{
        timer::split_windows++;
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
    MPI_Sendrecv(worker.dos.data()   , (int) worker.dos.size()   , MPI_DOUBLE, dn, 0, dos_up.data(), E_size_up*M_size_up, MPI_DOUBLE, up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(worker.E_bins.data(), (int) worker.E_bins.size(), MPI_DOUBLE, dn, 0, E_bins_up.data(), E_size_up, MPI_DOUBLE, up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(worker.M_bins.data(), (int) worker.M_bins.size(), MPI_DOUBLE, dn, 0, M_bins_up.data(), M_size_up, MPI_DOUBLE, up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (worker.E_bins.maxCoeff() < E_bins_up.minCoeff() && up != 0 ){
        cout << "Error in merge_windows(). Windows do not overlap." << endl;
        exit(3);
    }
    if (worker.E_bins.maxCoeff() > E_bins_up.maxCoeff() && up != 0 ){
        cout << "Error in merge_windows(). Overlap too large." << endl;
        exit(4);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    Vector3d v1A, v2A, v1B, v2B, uA, uB; //Vectors connecting adjacent 3 orthogonal points on DOS
    VectorXd sum(max(E_size,E_size_up));
    sum.fill(0);
    int col = 0, num;
    int x, y; //Coordinates closest to i, j on dos of neighbor above.
    for (int i = 0; i < E_size - 1; i++) {
        num = 0;
        for (int j = 0; j < M_size - 1; j++) {
            if (worker.E_bins(i) < E_bins_up.minCoeff()) { continue;}
            if (worker.E_bins(i) > E_bins_up.maxCoeff()) { continue;}
            if (worker.M_bins(j) > M_bins_up.maxCoeff()) { continue;}
            if (worker.M_bins(j) < M_bins_up.minCoeff()) { continue;}
            if (worker.dos(i,j) == 0) { continue;}

            x = math::upper_bound(E_bins_up.data(), worker.E_bins(i), E_size_up);
            y = math::upper_bound(M_bins_up.data(), worker.M_bins(i), M_size_up);
            v1A << worker.E_bins(i + 1) - worker.E_bins(i)  , 0, worker.dos(i + 1, j) - worker.dos(i, j);
            v1B << E_bins_up(x + 1) - E_bins_up(x)          , 0, dos_up(x + 1, y) - dos_up(x,y);
            v2A << 0, worker.M_bins(j + 1) - worker.M_bins(j)   , worker.dos(i, j+1) - worker.dos(i, j);
            v2B << 0, M_bins_up(y + 1) - M_bins_up(y)           , dos_up(x, y+1) - dos_up(x,y);
            uA = v1A.cross(v2A).normalized();
            uB = v1B.cross(v2B).normalized();
            sum(col) += uA.dot(uB);
            num++;
        }
        sum(col) /= max(1,num);
        col++;
    }
    MPI_Barrier(MPI_COMM_WORLD);


    int E_merge_idx; //Index of merging point
    int M_merge_idx; //Index of merging point
    sum.maxCoeff(&E_merge_idx);
    worker.dos.row(E_merge_idx).maxCoeff(&M_merge_idx);

    //Find the merge coordinates on the upper dos
    int E_merge_idx_up = math::upper_bound(E_bins_up.data(), worker.E_bins(E_merge_idx), E_size_up);
    int M_merge_idx_up = math::upper_bound(M_bins_up.data(), worker.M_bins(M_merge_idx), M_size_up);

    double diff;

    cout << setprecision(2)<<endl;
    for (int w = 0; w < worker.world_size - 1; w++) {
        if (worker.world_ID == w) {
            //Adjust the dos above and send up
            diff = worker.dos(E_merge_idx, M_merge_idx) - dos_up(E_merge_idx_up, M_merge_idx_up);
//            cout << "ID: " << worker.world_ID
//                 << " Diff " << diff << " (" << worker.dos(E_merge_idx, M_merge_idx) << " - " << dos_up(E_merge_idx_up, M_merge_idx_up) << ")"
//                 << " Merge [" << worker.E_bins(E_merge_idx) << "(" <<E_merge_idx<< ")"  << " " << worker.M_bins(M_merge_idx) << "(" << M_merge_idx << ")" << "]"
//                 << "   Merge up [" << E_merge_idx_up << " " << M_merge_idx_up << "]"
//                 << endl;
//            cout << "Sum: " << setprecision(2) <<sum.transpose() << endl;
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
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }



    //Now all workers have their dos aligned with the others. Now just send all to master?
    //Start by trimming off excess
    //Use dos_temp for temporary storage
    worker.dos_temp = worker.dos;
    worker.X_bins   = worker.E_bins;

    int E_merge_idx_dn;
    MPI_Sendrecv(&E_merge_idx_up, 1, MPI_INT, up, 0, &E_merge_idx_dn, 1, MPI_INT, dn, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "ID: " << worker.world_ID << " E_merge_idx: " << E_merge_idx << " E_merge_idx_dn: " << E_merge_idx_dn << endl;
    if (worker.world_ID == 0){
        worker.dos      = worker.dos_temp.topRows(E_merge_idx);
        worker.E_bins   = worker.X_bins.head(E_merge_idx);
        E_size          = (int) worker.dos.rows();

    }else if (worker.world_ID == worker.world_size - 1){
       worker.dos       = worker.dos_temp.bottomRows(E_size - E_merge_idx_dn);
       worker.E_bins    = worker.X_bins.tail(E_size - E_merge_idx_dn);
       E_size           = (int) worker.dos.rows();
    }else{
        worker.dos      = worker.dos_temp.middleRows(E_merge_idx_dn, E_merge_idx - E_merge_idx_dn );
        worker.E_bins   = worker.X_bins.segment(E_merge_idx_dn, E_merge_idx - E_merge_idx_dn );
        E_size          = (int) worker.dos.rows();
    }


    //Now send all to world_ID == 0, to do the merging
    MatrixXd dos_recv;
    VectorXd E_bins_recv;
    int E_new_size;
    for (int w = 1; w < worker.world_size; w++) {
        if (worker.world_ID == 0){
            MPI_Recv(&E_size_up, 1, MPI_INT, w, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&M_size_up, 1, MPI_INT, w, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            dos_recv.resize(E_size_up, M_size_up);
            E_bins_recv.resize(E_size_up);

            MPI_Recv(dos_recv.data()   , (int)dos_recv.size()    , MPI_DOUBLE, w, 3, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(E_bins_recv.data(), (int)E_bins_recv.size() , MPI_DOUBLE, w, 4, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            E_new_size = (int)worker.dos.rows() + (int)dos_recv.rows();

            worker.dos.conservativeResize(E_new_size, M_size_up);
            worker.E_bins.conservativeResize(E_new_size , 1);
            worker.dos.bottomRows(E_size_up)        = dos_recv;
            worker.E_bins.tail(E_size_up)           = E_bins_recv;

        }else if(worker.world_ID == w){
            MPI_Send(&E_size     , 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&M_size     , 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
            MPI_Send(worker.dos.data(), (int)worker.dos.size(), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
            MPI_Send(worker.E_bins.data()  , (int)worker.E_bins.size()  , MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if(worker.world_ID == 0){
        cout << "Total DOS "<< endl;
        cout << worker.dos << endl;
        E_size = worker.dos.rows();
        M_size = worker.dos.cols();

    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&E_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    worker.dos.conservativeResize(E_size,M_size);
    worker.E_bins.conservativeResize(E_size);
    MPI_Bcast(worker.dos.data(), E_size*M_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(worker.E_bins.data(), E_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    cout << "Merge Successful" << endl;

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
