#ifndef CONSTANTS_H   // if x.h hasn't been included yet
#define CONSTANTS_H   //  #define this so the compiler knows it has been included
#include <Eigen/Dense>
#include <Eigen/Core>
using namespace std;
using namespace Eigen;
namespace EMC_constants {
	//Evolutionary parameters
	const int M 				= 4;  				//Number of populations in a species (= threads in OpenMP)
	const int N 				= 20;				//Number of individuals per population
	const int N_best			= 4;				//Number of individuals in "hall of fame". Best individuals of all time (per population)
	extern int geneLength;          				//Number of bits per gene (The number of possible values for a parameter is 2^geneLength-1)
	extern int nGenes;								//Number of parameters in your model. This is set in datafiles.cpp, inData::inData
	extern int genomeLength;						//Number of bits for all genes, This is set in datafiles.cpp, inData::inData
    const int maxbits           = 60;

    const  int max_generations 	= (int)1e5;			//Number of generations to run the simulation

    const int print_freq        = 500;
	const int store_freq        = 10;
    const int check_freq        = 100;
	const int num_check_history = 20;
	const double Tmin 			= 0.0001;			//Minimum temperature of the ladder. Preferrably between close to 1
	const double Tmax 			= 50;				//Maximum temperature of the ladder. Preferrably around H_max
	const int seed 				= 8;				//Seed for the random number generator
	const long double log_param 		= 1.3;
	const long double log_const 		= (long double) (1/log(log_param));
	
	//Probabilities genetic operators
	const double qmig = 0.1;						//The migration probability vs evolution probability
	const double qm = 0.1;							//The mutation probability vs crossover (1-qm for crossover)
	const double qma = 0.1;//0.005;					//elite mutation rate (qm*(1-qma) for regular mutation)
	const int r_num = 10;							//Number of points to sample on snooker crossover (10-100 is ok, higher is slower but more thorough)
	
	//Probabilities for smart copy crossover
	const double P0 = 0.1;							//If parents have the  SAME bit, reverse with probability p0
	const double P1 = 0.4;							//If parents have DIFFERENT bit, reverse good bit with probability p1
	const double P2 = 0.7;							//If parents have DIFFERENT bit, reverse bad  bit with probability p2
	
	//Probability matrix for smart copy crossover (Liang and Wong, 2000)
	const double p_matrix[] = { P0*P0 + (1 - P0)*(1 - P0), 2 * P0*(1 - P0), P1*(1 - P2) + P2*(1 - P1), P1*P2 + (1 - P1)*(1 - P2) };
//	extern ArrayXd	upper_bound;							//upper boundary of parameters
//	extern ArrayXd	lower_bound;
};




//
//class boundaries {
//public:
//	boundaries(){
//		all_bounds.resize(nGenes, 2);
//		upper_bound.resize(nGenes, 1);
//		lower_bound.resize(nGenes, 1);
//	}
//	MatrixXd all_bounds;							//Lower and upper boundary of parameters
//	ArrayXd	upper_bound;							//upper boundary of parameters
//	ArrayXd	lower_bound;							//lower boundary of parameters
//};
//extern boundaries bounds;

#endif
