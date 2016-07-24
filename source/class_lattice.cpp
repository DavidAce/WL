//
// Created by david on 2016-07-24.
//
#include "class_lattice.h"
#include "class_worker.h"
#include "randomNumbers.h"

void class_lattice::randomize_lattice() {
    for(int i = 0; i < constants::L ; i++){
        for(int j= 0; j < constants::L ; j++){
            L(i,j) = 2 * rn::uniform_integer(0,1) - 1;
        }
    }
}

//Function for printing the lattice. Very easy if L is an Eigen matrix.
std::ostream &operator<<(std::ostream &os, const class_lattice &lattice){
    os << lattice.L << std::endl;
    return os;
}