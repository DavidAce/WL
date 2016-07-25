//
// Created by david on 2016-07-24.
//

#include "randomNumbers.h"
namespace rn{
    std::mt19937 rng;
    double uniform_double(const double &lowerLimit, const double &upperLimit) {
        std::uniform_real_distribution<> dice(lowerLimit, upperLimit);
        return dice(rng);
    }

    int uniform_integer(const int &lowerLimit, const int &upperLimit) {
        std::uniform_int_distribution<> dice(lowerLimit, upperLimit);
        return dice(rng);
    }
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