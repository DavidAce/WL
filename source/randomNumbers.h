//
// Created by david on 2016-07-24.
//

#ifndef WL_RANDOMNUMBERS_H
#define WL_RANDOMNUMBERS_H
#include <random>
#include <iostream>
#include "constants.h"
using namespace std;
namespace rn{
    //typedef std::mt19937 RNGType;
    //RNGType rng;
    //Random functions
    extern std::mt19937 rng;
    extern std::uniform_int_distribution<>  rand_int_L;
    extern std::uniform_int_distribution<>  rand_int_1;
    extern std::uniform_real_distribution<> rand_real_1;

    inline double uniform_double_1(){
        return rand_real_1(rng);
    }

    inline int uniform_integer_L(){
        return rand_int_L(rng);
    }

    inline int uniform_integer_1(){
        return rand_int_1(rng);
    }

    extern double gaussian_truncated(const double &lowerLimit, const double &upperLimit, const double &mean, const double &std) ;

}

#endif //WL_RANDOMNUMBERS_H
