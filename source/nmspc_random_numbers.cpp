//
// Created by david on 2016-07-24.
//

#include "nmspc_random_numbers.h"
namespace rn{
    std::mt19937 rng;

    double gaussian_truncated(const double lowerLimit, const double upperLimit, const double mean, const double std) {
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