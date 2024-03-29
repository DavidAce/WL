//
// Created by david on 9/1/16.
//

#include "class_WL_statistics.h"
#include "general/nmspc_logger.h"
#include "IO/class_WL_read_data.h"
// #include <eigen3/unsupported/Eigen/CXX11/Tensor>

void class_stats::load_thermo_data(class_worker &worker) {
    indata in;
    in.load_full(*this, worker);
}

void class_stats::compute(class_worker &worker) {
    double R  = constants::simulation_reps + constants::bootstrap_reps;
    E_avg     = E.rowwise().mean();
    M_avg     = M.rowwise().mean();
    s_avg     = s.rowwise().mean();
    c_avg     = c.rowwise().mean();
    m_avg     = m.rowwise().mean();
    u_avg     = u.rowwise().mean();
    f_avg     = f.rowwise().mean();
    x_avg     = x.rowwise().mean();
    dos1D_avg = dos1D.rowwise().mean();
    dos_avg   = math::mean_depthwise(dos);
    D_avg     = math::mean_depthwise(D);
    F_avg     = math::mean_depthwise(F);
    P_avg     = math::mean_depthwise(P);

    s_err     = ((s.colwise() - s_avg).cwiseAbs2().rowwise().sum() / R).cwiseSqrt();
    c_err     = ((c.colwise() - c_avg).cwiseAbs2().rowwise().sum() / R).cwiseSqrt();
    m_err     = ((m.colwise() - m_avg).cwiseAbs2().rowwise().sum() / R).cwiseSqrt();
    u_err     = ((u.colwise() - u_avg).cwiseAbs2().rowwise().sum() / R).cwiseSqrt();
    f_err     = ((f.colwise() - f_avg).cwiseAbs2().rowwise().sum() / R).cwiseSqrt();
    x_err     = ((x.colwise() - x_avg).cwiseAbs2().rowwise().sum() / R).cwiseSqrt();
    dos1D_err = ((dos1D.colwise() - dos1D_avg).cwiseAbs2().rowwise().sum() / R).cwiseSqrt();
    dos_err   = math::err_depthwise(dos, dos_avg);
    D_err     = math::err_depthwise(D, D_avg);
    F_err     = math::err_depthwise(F, F_avg);
    P_err     = math::err_depthwise(P, P_avg);
}
