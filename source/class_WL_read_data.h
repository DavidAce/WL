//
// Created by david on 9/1/16.
//

#ifndef WL_CLASS_WL_READ_DATA_H
#define WL_CLASS_WL_READ_DATA_H
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <string>
#include "nmspc_random_numbers.h"
#include "nmspc_WL_constants.h"
#include "nmspc_WL_parallel_algorithms.h"
#include "class_WL_worker.h"
#include "class_WL_statistics.h"

const int debug_load_full_thermo = 1;
const int debug_load_full_statis = 1;
using namespace Eigen;
using namespace std;
class indata {
private:
    string folder;
public:
    indata();
    void load_random_section(class_worker &worker);
    void load_full(class_worker &worker);
    void load_full(class_stats &stats, class_worker &worker);
    ArrayXXd read_file(string file);

};


#endif //WL_CLASS_WL_READ_DATA_H
