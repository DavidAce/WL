//
// Created by david on 2016-07-24.
//

#ifndef WL_NMSPC_RANDOM_NUMBERS_H
#define WL_NMSPC_RANDOM_NUMBERS_H
#include <random>
#include <iostream>
using namespace std;
namespace rn{
    //typedef std::mt19937 RNGType;
    //RNGType rng;
    //Random functions
    extern std::mt19937 rng;

    inline int __attribute__((hot)) uniform_integer_1(){
        std::uniform_int_distribution<>  rand_int(0, 1);
        return rand_int(rng);
    }

    inline int __attribute__((hot)) uniform_integer(const int min, const int max){
        std::uniform_int_distribution<>  rand_int(min,max);
        return rand_int(rng);
    }

    inline double __attribute__((hot)) uniform_double_1(){
        std::uniform_real_distribution<>  rand_real(0,1);
        return rand_real(rng);
    }

    extern double gaussian_truncated(const double lowerLimit, const double upperLimit, const double mean, const double std) ;

}

#endif //WL_NMSPC_RANDOM_NUMBERS_H
