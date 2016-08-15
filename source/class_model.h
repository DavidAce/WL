//
// Created by david on 2016-07-24.
//

#ifndef WL_CLASS_MODEL_H
#define WL_CLASS_MODEL_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include "randomNumbers.h"
#include "math_algorithms.h"
#include "constants.h"
using namespace Eigen;
using namespace std;

class class_model {
private:
    int randI,randJ; //Coordinates of random position;
public:
    class_model() {
        //lattice.resize(constants::L, constants::L);
        randomize_lattice();
    };
    Matrix<int, constants::L, constants::L> lattice;                           //The Lattice Data structure
    //double E_new, M_new;                      //Store values from MC-trial
    const static bool discrete_model = true;          //Toggle model type
    const static int J = 1;

    void randomize_lattice();
    inline void flip(){
        lattice(randI,randJ) *= -1;
    }

    double get_E();
    double get_M();

    inline int sum_neighbours(const int &i, const int &j){
        return            lattice(math::mod(i+1,constants::L), math::mod(j  ,constants::L))
                     +    lattice(math::mod(i  ,constants::L), math::mod(j+1,constants::L))
                     +    lattice(math::mod(i-1,constants::L), math::mod(j  ,constants::L))
                     +    lattice(math::mod(i  ,constants::L), math::mod(j-1,constants::L));
    }

    inline void make_new_state(const double &E, const double &M, double &E_trial, double &M_trial){
        randI = rn::uniform_integer_L();
        randJ = rn::uniform_integer_L();
        //Is this correct?
        E_trial = E + J * 2*lattice(randI,randJ)*sum_neighbours(randI,randJ);
        M_trial = M - 2*lattice(randI,randJ);
    }

    friend std::ostream &operator<<(std::ostream &, const class_model &);

};


#endif //WL_CLASS_MODEL_H
