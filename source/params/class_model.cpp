//
// Created by david on 2016-07-24.
//
#include <fstream>
#include "class_model.h"

using namespace std;

void class_model::randomize_lattice() {
    for(int i = 0; i < L ; i++){
        for(int j= 0; j < L ; j++){
            lattice(i,j) = 2 * rn::uniform_integer_1() - 1;
        }
    }
}


double class_model::get_E() {
    double E = 0;
    for(int i = 0; i < L; i ++){
        for(int j = 0; j < L; j ++) {
            E += - J * lattice(i,j) * sum_neighbours(i,j);
        }
    }

    return E/2.0;
}

double class_model::get_M() {
    return (double) lattice.sum();

}


//Function for printing the lattice. Very easy if L is an Eigen matrix.
std::ostream &operator<<(std::ostream &os, const class_model &model){
    os << model.lattice << std::endl;
    return os;
}
