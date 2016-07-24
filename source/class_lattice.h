//
// Created by david on 2016-07-24.
//

#ifndef WL_CLASS_LATTICE_H
#define WL_CLASS_LATTICE_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include "constants.h"
using namespace Eigen;

class class_lattice {
private:
public:
    class_lattice() {
        randomize_lattice();
    };
    MatrixXi lattice(constants::L,constants::L);
    void randomize_lattice();

};


#endif //WL_CLASS_LATTICE_H
