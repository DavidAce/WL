//
// Created by david on 9/1/16.
//
#pragma once
#include "class_WL_worker.h"
#include "general/nmspc_math_algorithms.h"
#include <Eigen/Core>

using namespace constants;
class class_thermodynamics {
    private:
    public:
    class_thermodynamics() = default;

    Eigen::ArrayXd  T;           // Temperature
    Eigen::ArrayXd  u;           // Internal energy
    Eigen::ArrayXd  m;           // Order paramater
    Eigen::ArrayXd  s;           // Entropy
    Eigen::ArrayXd  f;           // Free energy
    Eigen::ArrayXd  c;           // Specific heat
    Eigen::ArrayXd  x;           // Susceptibility
    Eigen::ArrayXd  dos_total1D; // Density of states in Energy space only
    Eigen::ArrayXXd D;           // Canonical Distribution
    Eigen::ArrayXXd P;           // Probability distribution
    Eigen::ArrayXXd F;           // Free energy as a function of order parameter
    void            compute(class_worker &worker);
};
