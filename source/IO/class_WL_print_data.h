//
// Created by david on 2016-07-26.
//

#ifndef WL_CLASS_PRINT_DATA_H
#define WL_CLASS_PRINT_DATA_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cstdlib>
#include <stdexcept>
#include <sys/stat.h>
#include <limits.h>
#include <errno.h>
#include <iostream>
#include "algorithm/class_WL_worker.h"
#include "algorithm/class_WL_thermo.h"
#include "algorithm/class_WL_statistics.h"



using namespace std;
using namespace Eigen;

class outdata {
private:
    string      folder;
    int         iteration = 0;
    int         precision = 12;
public:
    outdata() = default;  //Default constructor (does not set folder! make sure to set it yourself!

    void create_folder(string folder_name);
    void create_set_folder(string folder_name);
    void create_iteration_folder_master(const int &iter, const int &id);
    void create_iteration_folder_worker(const int &iter);
    void set_foldername_to_iteration(const int &iter);
    int mkdir_p(const char *path);
    //File streams

    void write_data_worker(class_worker &);
    void write_data_master(class_worker &);
    void write_data_thermo(class_thermodynamics &, const int &iter);
    void write_final_data(class_stats &stats, const int &id);
    void write_sample(class_worker &);
//    template<typename Derived>
//    void write_to_file(const MatrixBase<Derived> &data, string filename){
//        ofstream file(filename,ios::out | ios::trunc);
//        file << fixed << showpoint << setprecision(precision);
//        string      _coeffSeparator = "	";
//        IOFormat fmt(StreamPrecision, 0, _coeffSeparator);
//        file << data.format(fmt) << endl;
//        file.close();
//    }
    template<typename Derived>
    void write_to_file(const EigenBase<Derived> &data, string filename){
        ofstream file(filename,ios::out | ios::trunc);
        file << fixed << showpoint << setprecision(precision);
        string      _coeffSeparator = "	";
        IOFormat fmt(StreamPrecision, 0, _coeffSeparator);
        file << data.derived().format(fmt) << endl;
        file.close();
    }

    void write_to_file(const double data, string filename){
        ofstream file(filename,ios::out | ios::trunc);
        file << fixed << showpoint << setprecision(precision) << data;
        file.close();
    }

};

#endif //WL_CLASS_PRINT_DATA_H
