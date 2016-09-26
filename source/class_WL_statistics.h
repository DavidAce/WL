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
    ArrayXXd c_peak;
    ArrayXXd x_peak;
    ArrayXXd Tc_F;
    ArrayXXd Tc_D;
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
    ArrayXd c_peak_avg; //peak in c(T) vs T
    ArrayXd x_peak_avg; //peak in x(T) vs T
    ArrayXd Tc_F_avg; //Critical temperature average of free energy vs M
    ArrayXd Tc_D_avg; //Critical temperature average of free energy vs M
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
    ArrayXd c_peak_err;
    ArrayXd x_peak_err;
    ArrayXd Tc_F_err;
    ArrayXd Tc_D_err;
    ArrayXXd dos_err;
    ArrayXXd D_err;
    ArrayXXd F_err;
    ArrayXXd P_err;

};


#endif //WL_CLASS_WL_STATISTICS_H
