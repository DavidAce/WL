//
// Created by david on 2016-07-24.
//

#ifndef WL_CLASS_LATTICE_H
#define WL_CLASS_LATTICE_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include "constants.h"
using namespace Eigen;

class class_model {
private:
    int randI,randJ; //Coordinates of random position;
public:
    class_model() {
        lattice.resize(constants::L, constants::L);
        randomize_lattice();
    };
    MatrixXi lattice;                           //The Lattice Data structure
    //double E_new, M_new;                      //Store values from MC-trial
    const static bool discrete_model = true;          //Toggle model type
    void randomize_lattice();
    void flip();
    int mod(const int &,const int &);
    int sum_neighbours(const int &, const int &);

    double get_E();
    double get_M();

    void make_new_state(const double &, const double &, double &, double &);

    friend std::ostream &operator<<(std::ostream &, const class_model &);

};


#endif //WL_CLASS_LATTICE_H
