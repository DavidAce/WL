
#ifndef PERSONALITY_H   // if x.h hasn't been included yet...
#define PERSONALITY_H   //  #define this so the compiler knows it has been included
#include <Eigen/Dense>
#include <Eigen/Core>
#include <bitset>
#include <iostream>
#include <random>
#include <cfloat>
#include "constants.hpp"
#include "DNA.hpp"
#include "objective_function.hpp"

using namespace std;
using namespace EMC_constants;
using namespace Eigen;

class personality {
private:
    objective_function &obj_fun;
public:

	personality(objective_function &ref)
            :obj_fun(ref),
             born(0),
             genome(ref){
	};
	personality(objective_function &ref, bool toggle) :
            obj_fun(ref),
            born(0),
            genome(ref,toggle) {
	};
	long double H;					//Fitness, or energy
	double t; 							//temperature
	long double value;				//Actual fitting-value
	int born;							//Generation when DNA first emerged
	DNA genome;							//Contains binary and real representation of parameters
	friend ostream &operator<<(std::ostream &os, personality const &);

};
#endif
