
#include "population.hpp"


using namespace std;
using namespace EMC_constants;
using namespace Eigen;


std::ostream &operator<<(std::ostream &os, population const &pop) {
	for (int i = 0; i < N; i++) {
		os << pop.guys[i].H << endl;
	}
	return os;
}

void population::getFitness(personality& guy){
	guy.genome.update_parameters();
	guy.H = obj_fun.fitness(guy.genome.parameters);
}

void population::getFitness(const Array<long double, Dynamic, 1> &p, personality &guy){
    //Set all parameters at once with an ArrayXd
    guy.genome.set_parameters(p);
    guy.H = obj_fun.fitness(guy.genome.parameters);
}

//void population::getFitness4All(){
//        for (int i = 0; i < N; i++) {
//            cout << "Getting fitness for guy " << i << " from parameters: " << guys[i].genome.parameters.transpose() << endl;
//            guys[i].H = obj_fun.fitness(guys[i].genome.parameters);
//		}
//    }

//void population::wakeUpPop(){
////	wakeUpGuys();
////	getFitness4All();
//	//Wake up newguys by copying the old guys
//	for (int i = 0; i < N; i++) {
//		copy(newguys[i], guys[i]);
//	}
//	wakeUpBest();
//}

void population::wakeUpNewGuys(){
    for (int i = 0; i < N; i++) {
        copy(newguys[i], guys[i]);
    }
}

void population::wakeUpGuys() {
	ArrayXd T(N);
	//Initialize some temperature ladder, here logarithmic.
	T = LogSpaced(N,Tmax, Tmin);
    for (int i = 0; i < N; i++) {
		guys[i].t = T(i); 		//Assign temperatures produced above
    }
}


void population::wakeUpBest() {
	int i;
	int j = 0;
//	ArrayXi copied_guys(N_best);
    std::set<int> copied_guys;
	long double lowest_H;
	int lowest_i = 0; //The position i on the temperature ladder of the guy with lowest H
    while (copied_guys.size() < N_best) {
        lowest_H = 1e10;
        //Find the best guy yet among guys
        for (i = 0; i < N; i++) {
            //Check if i is in skip-list
            if (copied_guys.find(i) != copied_guys.end()){continue;}
            if (guys[i].H < lowest_H) {
                lowest_H = guys[i].H;
                lowest_i = i;
            }
        }
        //By now we should have a winner, copy him to bestguys and add him to skiplist
        copy(bestguys[N_best - j - 1], guys[lowest_i]);
        copied_guys.insert(lowest_i);
        j++;

    }
}


void population::copy(personality &destination, const personality &source) {
	//destination.born				= source.born;
	destination.H					= source.H;
	destination.t					= source.t;
	destination.genome.parameters	= source.genome.parameters;
	destination.value				= source.value;
	destination.genome.chromosomes  = source.genome.chromosomes;

	// std::copy(destination.genome.chromosomes, destination.genome.chromosomes + nGenes, source.genome.chromosomes);
}

void population::copy(DNA &destination, const DNA &source) {
	destination.chromosomes  = source.chromosomes;

	// std::copy(destination.chromosomes, destination.chromosomes + nGenes, source.chromosomes);
}


