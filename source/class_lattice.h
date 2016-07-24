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
        L.resize(constants::L, constants::L);
        randomize_lattice();
    };
    MatrixXi L;                                                                 //The Lattice Data structure
    void randomize_lattice();
    friend std::ostream &operator<<(std::ostream &, const class_lattice &);

};


#endif //WL_CLASS_LATTICE_H
