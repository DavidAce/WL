//
// Created by david on 2016-07-24.
//

#include "class_worker.h"
using namespace std;

double class_worker::uniform_double(RNGType *rn, const double &lowerLimit, const double &upperLimit) {
    uniform_real_distribution<> dice(lowerLimit, upperLimit);
    return dice(*rn);
}
int class_worker::uniform_integer(RNGType *rn, const int &lowerLimit, const int &upperLimit) {
    uniform_int_distribution<> dice(lowerLimit, upperLimit);
    return dice(*rn);
}

double class_worker::gaussian_truncated(RNGType *rn, const double &lowerLimit, const double &upperLimit, const double &mean, const double &std) {
    normal_distribution<double> distribution(mean,std);
    double ul = fmax(lowerLimit, upperLimit);
    double ll = fmin(lowerLimit, upperLimit);
    double number;
    while (true) {
        number = distribution(*rn);
        if (number >= ll && number <= ul) {
            return number;
        }
    }
}