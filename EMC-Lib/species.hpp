#ifndef SPECIES_H   // if x.h hasn't been included yet...
#define SPECIES_H   //  #define this so the compiler knows it has been included
#include <iostream>
#include <vector>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iomanip>
#include <string>
#include <fstream>
#include <time.h>

#include "population.hpp"
#include "constants.hpp"
#include "mymath.hpp"
#include "objective_function.hpp"


typedef std::chrono::high_resolution_clock Clock;
class counters {
public:
	int store_counter;
    int store_last_since;
	int generation;
	double simulation_time;
	Clock::time_point simulation_tic;
	Clock::time_point simulation_toc;
	double evolution_time;
    Clock::time_point evolution_tic;
    Clock::time_point evolution_toc;
};



class species{
private:
    objective_function &obj_fun;
	void zero_counters();
public:
    species(objective_function &ref): obj_fun(ref), pop(M,ref){
        zero_counters();
    }
    vector<population> pop;		//Array of separate populations to evolve independently
	counters count;

    long double champion_fitness();
	Array<long double, Dynamic,1> champion_value();
	int champion_number();
	long double latest_history_diff();
    Array<long double,Dynamic,1> fitness_history;
    void store_best_fitness();
    bool below_tolerance();
    void print_progress();
	void print_progress(bool);
	void copy(personality &, personality &);
};


#endif
