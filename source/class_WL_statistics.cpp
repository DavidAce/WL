//
// Created by david on 9/1/16.
//

#include "class_WL_statistics.h"
#include "class_WL_read_data.h"
class_stats::class_stats(const int &id, const int &size): world_ID(id), world_size(size){
}

void class_stats::load_thermo_data(class_worker &worker){
    indata in(world_ID, world_size);
    in.load_full(*this, worker);
}

void class_stats::compute(class_worker &worker){
    double B = constants::simulation_reps + constants::bootstrap_reps;
    if (worker.world_ID == 0 && debug_compute) {
        cout << "Resizing " << endl;
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(100000));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    s_avg.resize(constants::T_num);
    c_avg.resize(constants::T_num);
    u_avg.resize(constants::T_num);
    f_avg.resize(constants::T_num);
    x_avg.resize(constants::T_num);

    s_err.resize(constants::T_num);
    c_err.resize(constants::T_num);
    u_err.resize(constants::T_num);
    f_err.resize(constants::T_num);
    x_err.resize(constants::T_num);

    dos1D_avg.resize(worker.E_bins_total.size());
    dos1D_err.resize(worker.E_bins_total.size());
    if (worker.world_ID == 0 && debug_compute) {
        cout << "Computing 1D dos " << endl;
        cout.flush();
        std::this_thread::sleep_for(std::chrono::microseconds(10000));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_ID == 0){
        int e = 0;
        while( e < worker.E_bins_total.size()){
            dos1D_avg(e) = dos1D.row(e).mean();
            dos1D_err(e) = sqrt(1/(B-1) * ((dos1D.row(e).array() - dos1D_avg(e)) * (dos1D.row(e).array() - dos1D_avg(e) )).sum()) ;
//            MPI_Gather(dos1D_avg.data() + e, 1, MPI_DOUBLE, dos1D_avg.data()+e, 1, MPI_DOUBLE,0, MPI_COMM_WORLD);
//            MPI_Gather(dos1D_err.data() + e, 1, MPI_DOUBLE, dos1D_err.data()+e, 1, MPI_DOUBLE,0, MPI_COMM_WORLD);
            e++;
        }

        int t = 0;
        while (t < constants::T_num) {
            s_avg(t) = s.row(t).mean();
            c_avg(t) = c.row(t).mean();
            u_avg(t) = u.row(t).mean();
            f_avg(t) = f.row(t).mean();
            x_avg(t) = x.row(t).mean();

            s_err(t) = sqrt(1 / (B - 1) * ((s.row(t).array() - s_avg(t)) * (s.row(t).array() - s_avg(t))).sum());
            c_err(t) = sqrt(1 / (B - 1) * ((c.row(t).array() - c_avg(t)) * (c.row(t).array() - c_avg(t))).sum());
            u_err(t) = sqrt(1 / (B - 1) * ((u.row(t).array() - u_avg(t)) * (u.row(t).array() - u_avg(t))).sum());
            f_err(t) = sqrt(1 / (B - 1) * ((f.row(t).array() - f_avg(t)) * (f.row(t).array() - f_avg(t))).sum());
            x_err(t) = sqrt(1 / (B - 1) * ((x.row(t).array() - x_avg(t)) * (x.row(t).array() - x_avg(t))).sum());

//            MPI_Gather(s_avg.data() + t, 1, MPI_DOUBLE, s_avg.data() + t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            MPI_Gather(c_avg.data() + t, 1, MPI_DOUBLE, c_avg.data() + t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            MPI_Gather(u_avg.data() + t, 1, MPI_DOUBLE, u_avg.data() + t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            MPI_Gather(f_avg.data() + t, 1, MPI_DOUBLE, f_avg.data() + t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            MPI_Gather(x_avg.data() + t, 1, MPI_DOUBLE, x_avg.data() + t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//            MPI_Gather(s_err.data() + t, 1, MPI_DOUBLE, s_err.data() + t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            MPI_Gather(c_err.data() + t, 1, MPI_DOUBLE, c_err.data() + t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            MPI_Gather(u_err.data() + t, 1, MPI_DOUBLE, u_err.data() + t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            MPI_Gather(f_err.data() + t, 1, MPI_DOUBLE, f_err.data() + t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            MPI_Gather(x_err.data() + t, 1, MPI_DOUBLE, x_err.data() + t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            t++;
        }
    }
}


