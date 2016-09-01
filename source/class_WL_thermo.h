//
// Created by david on 9/1/16.
//

#ifndef WL_CLASS_WL_THERMO_H
#define WL_CLASS_WL_THERMO_H
//#include "matrix_extensions.h"
//#define EIGEN_MATRIXBASE_PLUGIN "matrix_extensions.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include "class_WL_worker.h"
#include "nmspc_math_algorithms.h"
using namespace Eigen;
using namespace std;
using namespace constants;
class class_thermodynamics {
private:

public:
    class_thermodynamics();
    ArrayXd T; //Temperature
    ArrayXd u; //Internal energy
    ArrayXd m; //Order paramater
    ArrayXd s; //Entropy
    ArrayXd f; //Free energy
    ArrayXd c; //Specific heat
    ArrayXd x; //Susceptibility
    ArrayXd dos_total1D; //Density of states in Energy space only
    void compute(class_worker &worker);
};


#endif //WL_CLASS_WL_THERMO_H
