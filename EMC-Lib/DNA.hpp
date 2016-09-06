
#ifndef DNA_H   // if x.h hasn't been included yet...
#define DNA_H   //  #define this so the compiler knows it has been included
#include <Eigen/Dense>
#include <Eigen/Core>
#include <bitset>
#include <vector>
#include <iostream>
#include <random>
#include <memory>
#include <climits>
#include <algorithm>
#include "constants.hpp"
#include "randomFunctions.hpp"
#include "objective_function.hpp"

using namespace std;
using namespace EMC_constants;
using namespace EMC_rnd;
using namespace Eigen;
class DNA {
private:
	double bin2dec(const int);
	bitset<maxbits> dec2bin(const int);

//	vector<bool> dec2bin(const int);
    objective_function &obj_fun;
public:
	DNA (objective_function &ref);
	DNA (objective_function &ref, bool toggle);
	vector< bitset<maxbits> > chromosomes; //Binary representation
	Array<long double, Dynamic,1> parameters;						  //Decimal representation

	bool operator== (const DNA& target);
	int operator()(int); 
	friend ostream &operator<<(std::ostream &os, DNA const &);
	void flip_loci(const int);
	void flip_loci(ArrayXi &);
	void copy_loci(const int, const int);

	void set_parameter(const int,const long double);
	void set_parameters(const Array<long double, Dynamic, 1>  &p);
	void update_parameters();
    void randomize_dna();
}; 

#endif
