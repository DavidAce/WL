//
// Created by david on 9/1/16.
//

#include "class_WL_statistics.h"
#include "class_WL_read_data.h"
//#include <eigen3/unsupported/Eigen/CXX11/Tensor>
class_stats::class_stats(const int &id, const int &size): world_ID(id), world_size(size){
}

void class_stats::load_thermo_data(class_worker &worker){
    indata in;
    in.load_full(*this, worker);
}

void class_stats::compute(class_worker &worker){
    double B = constants::simulation_reps + constants::bootstrap_reps;
    s_avg = s.rowwise().mean();
    c_avg = c.rowwise().mean();
    u_avg = u.rowwise().mean();
    f_avg = f.rowwise().mean();
    x_avg = x.rowwise().mean();
    dos1D_avg = dos1D.rowwise().mean();
    c_peak_avg= c_peak.rowwise().mean();



    s_err = ((s.colwise() - s_avg).cwiseAbs2().rowwise().sum()/(B-1)).cwiseSqrt();
    c_err = ((c.colwise() - c_avg).cwiseAbs2().rowwise().sum()/(B-1)).cwiseSqrt();
    u_err = ((u.colwise() - u_avg).cwiseAbs2().rowwise().sum()/(B-1)).cwiseSqrt();
    f_err = ((f.colwise() - f_avg).cwiseAbs2().rowwise().sum()/(B-1)).cwiseSqrt();
    x_err = ((x.colwise() - x_avg).cwiseAbs2().rowwise().sum()/(B-1)).cwiseSqrt();
    dos1D_err  = ((dos1D.colwise() -       dos1D_avg).cwiseAbs2().rowwise().sum()/(B-1)).cwiseSqrt();
    c_peak_err = ((c_peak.colwise() - c_peak_avg).cwiseAbs2().rowwise().sum()/(B-1)).cwiseSqrt();

}


