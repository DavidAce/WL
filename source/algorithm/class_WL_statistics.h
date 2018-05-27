//
// Created by david on 9/1/16.
//

#ifndef WL_CLASS_WL_STATISTICS_H
#define WL_CLASS_WL_STATISTICS_H
#include <Eigen/Dense>
#include <Eigen/Core>
#include <thread>
#include "class_WL_worker.h"
static const int debug_compute           =	1;

using namespace Eigen;
class class_stats {
private:
    int world_ID;
    int world_size;
public:
    class_stats(const int &id, const int &size);
    void load_thermo_data(class_worker &worker);
    void compute(class_worker &worker);
    ArrayXd T;

    ArrayXXd E;
    ArrayXXd M;

    ArrayXXd s;
    ArrayXXd c;
    ArrayXXd m;
    ArrayXXd u;
    ArrayXXd f;
    ArrayXXd x;
    ArrayXXd dos1D;
    vector<ArrayXXd> dos;
    vector<ArrayXXd> D;
    vector<ArrayXXd> F;
    vector<ArrayXXd> P;

    ArrayXd E_avg;
    ArrayXd M_avg;

    ArrayXd s_avg;
    ArrayXd c_avg;
    ArrayXd m_avg;
    ArrayXd u_avg;
    ArrayXd f_avg;
    ArrayXd x_avg;
    ArrayXd dos1D_avg;
    ArrayXXd dos_avg;
    ArrayXXd D_avg;
    ArrayXXd F_avg;
    ArrayXXd P_avg;

    ArrayXd s_err;
    ArrayXd c_err;
    ArrayXd m_err;
    ArrayXd u_err;
    ArrayXd f_err;
    ArrayXd x_err;
    ArrayXd dos1D_err;
    ArrayXXd dos_err;
    ArrayXXd D_err;
    ArrayXXd F_err;
    ArrayXXd P_err;

};


#endif //WL_CLASS_WL_STATISTICS_H
