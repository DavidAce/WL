//
// Created by david on 9/1/16.
//

#ifndef WL_CLASS_WL_READ_DATA_H
#define WL_CLASS_WL_READ_DATA_H
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include "nmspc_random_numbers.h"
#include "nmspc_WL_constants.h"
#include "class_WL_worker.h"
#include "class_WL_statistics.h"
#define MAXBUFSIZE  ((int) 1e6)

using namespace Eigen;
using namespace std;
class indata {
private:
    string folder;
    int   world_ID;
    int   world_size;
//    int   iteration;

public:
    indata(int &id, int &size);
    void load_random_section(class_worker &worker);
    void load_full(class_worker &worker);
    void load_full(class_stats &stats);

    MatrixXd read_file(string file);

};


#endif //WL_CLASS_WL_READ_DATA_H
