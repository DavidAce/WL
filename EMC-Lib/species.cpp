
#include "species.hpp"


using namespace std;
using namespace EMC_constants;
using namespace Eigen;

Array<long double, Dynamic,1> species::champion_value() {
	int best_index = 0;
	long double best_H = pop[0].bestguys[N_best - 1].H;
	for (int i = 1; i < M; i++) {
		if (pop[i].bestguys[N_best - 1].H < best_H) {
			best_index = i;
			best_H = pop[i].bestguys[N_best - 1].H;
		}
	}
	return pop[best_index].bestguys[N_best - 1].genome.parameters;
}

long double species::champion_fitness() {
	long double best_H = pop[0].bestguys[N_best - 1].H;
	for (int i = 1; i < M; i++) {
		if (pop[i].bestguys[N_best - 1].H < best_H) {
			best_H = pop[i].bestguys[N_best - 1].H;
		}
	}
	return best_H;
}

int species::champion_number() {
	int best_index = 0;
	long double best_H = pop[0].bestguys[N_best - 1].H;
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
            fitness_history.conservativeResize(2 * max(1,count.store_counter));
        }
        long double champ = champion_fitness();
        if (count.store_counter > 0){
            if(champ == fitness_history(count.store_counter-1)){
                count.store_last_since++;
            }else{
                fitness_history(count.store_counter) = champ;
                count.store_counter++;
                count.store_last_since = 0;
            }
            if(count.store_last_since == 5){ //Override
                fitness_history(count.store_counter) = champ;
                count.store_counter++;
                count.store_last_since = 0;
            }
        }else{
            fitness_history(count.store_counter) = champ;
            count.store_counter++;
            count.store_last_since = 0;
        }
    }
}

long double species::latest_history_diff() {

    int n = min(EMC_constants::num_check_history, count.store_counter);
    long double diff = 1;
    int from = count.store_counter - n;
    if (from > 0){
        diff = (fitness_history.segment(from-1, n) - fitness_history.segment(from,n)).mean();
    }
    return fabs(diff);
}

bool species::below_tolerance() {
    if (mod(count.generation, check_freq) == 0) {
        if (latest_history_diff() < obj_fun.tolerance) {
            return true;
        }
    }
    return false;
}

void species::print_progress(){
	count.generation++;
	if (mod(count.generation, print_freq) == 0) {
        IOFormat fmt(19, 0, "  ", "  ", "", "", "", "");
        cout << fixed << setprecision(19);
		cout << "Generation... " << setw(7) << count.generation << " | Current Best: ";
        cout << champion_value().transpose().eval().format(fmt) << "   | H = ";
        cout << setw(22)  << champion_fitness() << " | " << " diff: ";
		cout << setw(22)  << latest_history_diff();
        cout << endl;
		cout << flush;

	}
}

void species::print_progress(bool){
    IOFormat fmt(19, 0, "  ", "  ", "", "", "", "");
    cout << fixed << setprecision(19);
    cout << "Generation... " << setw(7) << count.generation << " | Current Best: ";
    cout << champion_value().transpose().eval().format(fmt) << "   | H = ";
    cout << setw(22)  << champion_fitness() << " | " << " diff: ";
    cout << setw(22)  << latest_history_diff();
    cout << endl;
    cout << flush;
}


void species::zero_counters() {
	count.evolution_time = 0;
	count.simulation_time = 0;
	count.store_counter = 0;
    count.store_last_since = 0;
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


