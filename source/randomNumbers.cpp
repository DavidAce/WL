//
// Created by david on 2016-07-24.
//

#include "randomNumbers.h"
namespace rn{
    std::mt19937 rng;

    std::uniform_int_distribution<>  rand_int_L(0, constants::L-1);
    std::uniform_int_distribution<>  rand_int_1(0, 1);
    std::uniform_real_distribution<> rand_real_1(0,1);

    double gaussian_truncated(const double &lowerLimit, const double &upperLimit, const double &mean, const double &std) {
        std::normal_distribution<double> distribution(mean,std);
        double ul = fmax(lowerLimit, upperLimit);
        double ll = fmin(lowerLimit, upperLimit);
        double number;
        while (true) {
            number = distribution(rng);
            if (number >= ll && number <= ul) {
                return number;
            }
        }
    }
}