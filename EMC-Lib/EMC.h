//
// Created by david on 9/2/16.
//

#ifndef EMC_CLION_EMC_H
#define EMC_CLION_EMC_H
#include <Eigen/Dense>
#include <Eigen/Core>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include "constants.hpp"
#include "randomFunctions.hpp"
#include "species.hpp"
#include "evolution.hpp"
#include "population.hpp"
#include "objective_function.hpp"
using namespace Eigen;
using namespace std;
using namespace std::chrono;

//Minimization should be a function that takes a pointer to a function,
//An int with the number of parameters to fit, and two arrays with boundaries min and max for each parameter.
void minimize(objective_function &obj_fun);

#endif //EMC_CLION_EMC_H
