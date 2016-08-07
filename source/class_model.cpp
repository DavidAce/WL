//
// Created by david on 2016-07-24.
//
#include <fstream>
#include "class_model.h"
#include "class_worker.h"
#include "randomNumbers.h"
using namespace std;

void class_model::randomize_lattice() {
    for(int i = 0; i < constants::L ; i++){
        for(int j= 0; j < constants::L ; j++){
            lattice(i,j) = 2 * rn::uniform_integer(0,1) - 1;
        }
    }
}

void class_model::flip(){
    lattice(randI,randJ) *= -1;
}

int class_model::mod(const int &a, const int &b){
    if (b < 0) {
        return class_model::mod(-a, -b);
    }
    int ret = a % b;
    if (ret < 0) {
        ret += b;
    }
    return ret;
};

int class_model::sum_neighbours(const int &i, const int &j) {
    return       lattice(mod(i+1,constants::L), mod(j  ,constants::L))
            +    lattice(mod(i  ,constants::L), mod(j+1,constants::L))
            +    lattice(mod(i-1,constants::L), mod(j  ,constants::L))
            +    lattice(mod(i  ,constants::L), mod(j-1,constants::L));
}

double class_model::get_E() {
    double E = 0;
    for(int i = 0; i < constants::L; i ++){
        for(int j = 0; j < constants::L; j ++) {
            E += - J * lattice(i,j) * sum_neighbours(i,j);
        }
    }

    return E/2.0;
}

double class_model::get_M() {
    return (double) lattice.sum();
}

void class_model::make_new_state(const double &E, const double &M, double &E_trial, double &M_trial){

    randI = rn::uniform_integer(0, constants::L - 1);
    randJ = rn::uniform_integer(0, constants::L - 1);
    //Is this correct?
    E_trial = E + J * 2*lattice(randI,randJ)*sum_neighbours(randI,randJ);
    M_trial = M - 2*lattice(randI,randJ);

}


//Function for printing the lattice. Very easy if L is an Eigen matrix.
std::ostream &operator<<(std::ostream &os, const class_model &model){
    os << model.lattice << std::endl;
    return os;
}