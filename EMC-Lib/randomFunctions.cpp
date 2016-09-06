#include "randomFunctions.hpp"


namespace EMC_rnd{
	std::mt19937 rng;
	std::uniform_int_distribution<>  rand_int_1(0, 1);
	std::uniform_real_distribution<> rand_real_1(0,1);
};


//double uniform_double(RNGType *rn, const double lowerLimit, const double upperLimit) {
//	uniform_real_distribution<> dice(lowerLimit, upperLimit);
//	return dice(*rn);
//}
//int uniform_integer(RNGType *rn, const int lowerLimit, const int upperLimit) {
//	uniform_int_distribution<> dice(lowerLimit, upperLimit);
//	return dice(*rn);
//}

//long double gaussian_truncated(RNGType *rn, const long double lowerLimit, const long double upperLimit, const double mean, const double std) {
//	normal_distribution<long double> distribution(mean,std);
//	long double ul = fmax(lowerLimit, upperLimit);
//	long double ll = fmin(lowerLimit, upperLimit);
//	long double number;
//	while (true) {
//		number = distribution(*rn);
//		if (number >= ll && number <= ul) {
//			return number;
//		}
//	}
//}
//class gaussian_truncated {
//	std::default_random_engine generator;
//	std::normal_distribution<double> distribution;
//	double min;
//	double max;
//public:
//	gaussian_truncated(double mean, double stddev, double min, double max) :
//		distribution(mean, stddev), min(min), max(max)
//	{}
//
//	double operator ()(double mean, double stddev, double min, double max) {
//		while (true) {
//			double number = this->distribution(generator);
//			if (number >= this->min && number <= this->max)
//				return number;
//		}
//	}
//};