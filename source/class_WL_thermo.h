//
// Created by david on 9/1/16.
//

#ifndef WL_CLASS_WL_THERMO_H
#define WL_CLASS_WL_THERMO_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include "class_WL_worker.h"
#include "nmspc_math_algorithms.h"
#include "../EMC-Lib/EMC.h"

using namespace Eigen;
using namespace std;
using namespace constants;
class class_thermodynamics {
private:

public:
    class_thermodynamics(){};
    ArrayXd T; //Temperature
    ArrayXd u; //Internal energy
    ArrayXd m; //Order paramater
    ArrayXd s; //Entropy
    ArrayXd f; //Free energy
    ArrayXd c; //Specific heat
    ArrayXd x; //Susceptibility
    ArrayXd dos_total1D; //Density of states in Energy space only
    ArrayXd c_peak;
    ArrayXd x_peak;
    ArrayXd Tc_F; //Free energy critical temperature (Tc)
    ArrayXXd D; //Canonical Distribution
    ArrayXXd P; //Probability distribution
    ArrayXXd F; //Free energy as a function of order parameter
    void compute(class_worker &worker);
    void get_c_peak(class_worker &worker);
    void get_x_peak(class_worker &worker);
    void get_Tc_free_energy(class_worker &worker);
};


#endif //WL_CLASS_WL_THERMO_H
