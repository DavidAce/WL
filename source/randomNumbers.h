//
// Created by david on 2016-07-24.
//

#ifndef WL_RANDOMNUMBERS_H
#define WL_RANDOMNUMBERS_H
#include <random>

namespace rn{
    //typedef std::mt19937 RNGType;
    //RNGType rng;
    //Random functions
    extern std::mt19937 rng;

    extern double uniform_double(const double &lowerLimit, const double &upperLimit);

    extern int uniform_integer(const int &lowerLimit, const int &upperLimit) ;

    extern double gaussian_truncated(const double &lowerLimit, const double &upperLimit, const double &mean, const double &std) ;

}

#endif //WL_RANDOMNUMBERS_H
