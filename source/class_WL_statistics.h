//
// Created by david on 9/1/16.
//

#ifndef WL_CLASS_WL_STATISTICS_H
#define WL_CLASS_WL_STATISTICS_H
#include <Eigen/Dense>
#include <Eigen/Core>
#include <thread>
#include "class_WL_worker.h"
static const int debug_compute           =	0;

using namespace Eigen;
class class_stats {
private:
    int world_ID;
    int world_size;
public:
    class_stats(const int &id, const int &size);
    void load_thermo_data(class_worker &worker);
    void compute(class_worker &worker);
    MatrixXd T;
    MatrixXd s;
    MatrixXd c;
    MatrixXd u;
    MatrixXd f;
    MatrixXd x;
    MatrixXd dos1D;

    MatrixXd E;
    MatrixXd M;

    ArrayXd s_avg;
    ArrayXd c_avg;
    ArrayXd u_avg;
    ArrayXd f_avg;
    ArrayXd x_avg;
    ArrayXd dos1D_avg;

    ArrayXd s_err;
    ArrayXd c_err;
    ArrayXd u_err;
    ArrayXd f_err;
    ArrayXd x_err;
    ArrayXd dos1D_err;

};


#endif //WL_CLASS_WL_STATISTICS_H
