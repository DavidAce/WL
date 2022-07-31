//
// Created by david on 9/1/16.
//

#ifndef WL_CLASS_WL_STATISTICS_H
#define WL_CLASS_WL_STATISTICS_H
#include "class_WL_worker.h"
#include <Eigen/Core>
#include <thread>
static const int debug_compute = 1;

class class_stats {
    public:
    class_stats() = default;
    void           load_thermo_data(class_worker &worker);
    void           compute(class_worker &worker);
    Eigen::ArrayXd T;

    Eigen::ArrayXXd E;
    Eigen::ArrayXXd M;

    Eigen::ArrayXXd              s;
    Eigen::ArrayXXd              c;
    Eigen::ArrayXXd              m;
    Eigen::ArrayXXd              u;
    Eigen::ArrayXXd              f;
    Eigen::ArrayXXd              x;
    Eigen::ArrayXXd              dos1D;
    std::vector<Eigen::ArrayXXd> dos;
    std::vector<Eigen::ArrayXXd> D;
    std::vector<Eigen::ArrayXXd> F;
    std::vector<Eigen::ArrayXXd> P;

    Eigen::ArrayXd E_avg;
    Eigen::ArrayXd M_avg;

    Eigen::ArrayXd  s_avg;
    Eigen::ArrayXd  c_avg;
    Eigen::ArrayXd  m_avg;
    Eigen::ArrayXd  u_avg;
    Eigen::ArrayXd  f_avg;
    Eigen::ArrayXd  x_avg;
    Eigen::ArrayXd  dos1D_avg;
    Eigen::ArrayXXd dos_avg;
    Eigen::ArrayXXd D_avg;
    Eigen::ArrayXXd F_avg;
    Eigen::ArrayXXd P_avg;

    Eigen::ArrayXd  s_err;
    Eigen::ArrayXd  c_err;
    Eigen::ArrayXd  m_err;
    Eigen::ArrayXd  u_err;
    Eigen::ArrayXd  f_err;
    Eigen::ArrayXd  x_err;
    Eigen::ArrayXd  dos1D_err;
    Eigen::ArrayXXd dos_err;
    Eigen::ArrayXXd D_err;
    Eigen::ArrayXXd F_err;
    Eigen::ArrayXXd P_err;
};

#endif // WL_CLASS_WL_STATISTICS_H
