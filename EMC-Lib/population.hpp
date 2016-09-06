
#ifndef POPULATION_H   // if x.h hasn't been included yet...
#define POPULATION_H   //  #define this so the compiler knows it has been included
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "personality.hpp"
#include "objective_function.hpp"
#include "constants.hpp"
using namespace std;
using namespace EMC_constants;
using namespace Eigen;
//class inData;		//Forward declaration

class paramLine{
	private:
        objective_function &obj_fun;
        ArrayXd o;
        ArrayXd d;
	public:
		paramLine(objective_function &ref)
                :obj_fun(ref),
                 o(nGenes),
                 d(nGenes){};
		//Always use "Through" first, to set "o" and "d".
		void Through(ArrayXd &p1, ArrayXd &p2){ //Points in parameter space p1 and p2
			for(int i = 0;i < nGenes; i++){
				o(i) = p1(i);
				d(i) = p2(i) - p1(i);
			}
		}
		ArrayXd pointAt(double r){
			ArrayXd v(nGenes);
			for(int i = 0;i < nGenes; i++){
				v(i) = o(i)+d(i)*r;
			}
			return v;
		}
		double distance(ArrayXd &p1, ArrayXd &p2){
			return sqrt((p2-p1).square().sum());
		}
		double line_max() { //Get the largest r within upper boundary Bu and lower boundary Bl
			return (obj_fun.upper_bound - o).cwiseQuotient(d).minCoeff();
		}
	double line_min() { //Get the largest r within upper boundary Bu and lower boundary Bl
		return -(obj_fun.lower_bound - o).cwiseQuotient(-d).minCoeff();
	}
};

class population{
private:
//	void getFitness(personality& guy);
	void wakeUpPop ();
	void wakeUpGuys();
	void wakeUpBest();
    void getFitness4All();
	objective_function &obj_fun;
public:
	population(objective_function &ref)
			:obj_fun        (ref),
			 guys 	        (N,ref),
			 newguys        (N,ref),
			 bestguys       (N_best,ref),
			 snookerGuys    (r_num,ref),
             line           (ref),
			 generation(0){
        wakeUpPop();
    };
    vector<personality> guys; 			     //Make an array of N guys
	vector<personality> newguys; 			 //Make a temporary array of N guinneapigs
	vector<personality> bestguys;            //Make an array of N/10 good performers
	vector<personality> snookerGuys;         // (r_num, personality(true));//Make an array of r_num snooker guys

    paramLine line;	//for doing the snooker crossover
    int generation;  //Number of generations for this population

    void getFitness(personality &guy);
    void getFitness(const ArrayXd &point, personality &guy);

    void copy(personality &destination, const personality & source);
    void copy(DNA& destination, const DNA& source);


	int operator()() { //Return the bit at a.
		return 0;
	}
	friend ostream &operator<<(std::ostream &os, population const &);


};

#endif
