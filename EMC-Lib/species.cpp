
#include "species.hpp"


using namespace std;
using namespace EMC_constants;
using namespace Eigen;

ArrayXd species::champion_value() {
	int best_index = 0;
	double best_H = pop[0].bestguys[N_best - 1].H;
	for (int i = 1; i < M; i++) {
		if (pop[i].bestguys[N_best - 1].H < best_H) {
			best_index = i;
			best_H = pop[i].bestguys[N_best - 1].H;
		}
	}
	return pop[best_index].bestguys[N_best - 1].genome.parameters;
}

double species::champion_fitness() {
	double best_H = pop[0].bestguys[N_best - 1].H;
	for (int i = 1; i < M; i++) {
		if (pop[i].bestguys[N_best - 1].H < best_H) {
			best_H = pop[i].bestguys[N_best - 1].H;
		}
	}
	return best_H;
}

int species::champion_number() {
	int best_index = 0;
	double best_H = pop[0].bestguys[N_best - 1].H;
	for (int i = 1; i < M; i++) {
		if (pop[i].bestguys[N_best - 1].H < best_H) {
			best_index = i;
			best_H = pop[i].bestguys[N_best - 1].H;
		}
	}
	return best_index;
}

void species::store_best_fitness() {
    if (mod(count.generation, store_freq) == 0) {
        if(fitness_history.size() <= count.store_counter) {
            fitness_history.conservativeResize(2 * count.store_counter);
        }
        fitness_history(count.store_counter) = champion_fitness();
        count.store_counter++;


    }
}

double species::latest_history_diff(const int & num) {
    int n = min(num, count.store_counter);
    double diff = 1;
    int from = count.store_counter - n;

    if (count.store_counter > 2) {
//        cout << n << " " << from << endl;
        diff = fitness_history.segment(from, n).mean() - fitness_history.segment(from, n).minCoeff();
    }
    return fabs(diff);
}

bool species::below_tolerance() {
    if (mod(count.generation, check_freq) == 0) {
        if (latest_history_diff(10) < obj_fun.tolerance) {
            return true;
        }
    }
    return false;
}

void species::print_progress(){
	count.generation++;
	if (mod(count.generation, print_freq) == 0) {
		cout << fixed << setprecision(12);
		cout << "Generation... " << setw(7) << count.generation << " | Current Best: ";
        cout << setw(12) << champion_value() << "   (H = " << champion_fitness() << ") ";
		cout << setw(12) << "diff: " << latest_history_diff(10);
        cout << endl;
		cout << flush;

	}
}



void species::zero_counters() {
	count.evolution_time = 0;
	count.simulation_time = 0;
	count.store_counter = 0;
	count.generation = 0;
}

void species::copy(personality &destination, personality &source) {
	destination.born = source.born;
	destination.H = source.H;
	destination.t = source.t;
	destination.genome.parameters = source.genome.parameters;
	destination.value = source.value;
	destination.genome.chromosomes = source.genome.chromosomes;
}


