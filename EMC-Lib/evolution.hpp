
#ifndef EVOLUTION_H
#define EVOLUTION_H  

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "constants.hpp"
#include "randomFunctions.hpp"
#include "mymath.hpp"
#include "population.hpp"
#include "species.hpp"
#include "objective_function.hpp"
using namespace std;
using namespace Eigen;
using namespace EMC_constants;
using namespace EMC_rnd;
class population;	//Forward declaration
class personality;	//Forward declaration

extern void evolve (population &);
extern void exchange(population &);
extern void migration(species &);
extern void insertguy(population &, int , int );
extern void find_elite(population &);

extern void roulette_select(personality [], int [], double *, double );
extern void bitselector_smartCopy(population &, int [], int []);

extern void mutation            (population &);
extern void mutation_elite      (population &);
extern void crossover           (population &);
extern void crossover_elite     (population &);
extern void crossover_smartCopy (population &);
extern void crossover_snooker   (population &);

#endif