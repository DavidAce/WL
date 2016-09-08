//
// Created by david on 2016-07-26.
//

#ifndef WL_CLASS_PRINT_DATA_H
#define WL_CLASS_PRINT_DATA_H
#include <Eigen/Core>
#include <Eigen/Dense>

#include <experimental/filesystem>
#include <cstdlib>

#include <iostream>
#include "class_WL_worker.h"
#include "class_WL_thermo.h"
#include "class_WL_statistics.h"
using namespace std;
using namespace Eigen;
namespace fs = std::experimental::filesystem;

class outdata {
private:
    string      folder;
    fs::path path;
    int   world_ID;
    int   iteration;
    int   precision = 10;
public:
    outdata(const int &id, const int &iter);
    outdata(const int &id);
    outdata();
    void create_folder();
    void set_folder(const int &iter);
    void create_and_set_folder(const int &iter);

    //File streams

    void write_data_worker(class_worker &);
    void write_data_master(class_worker &);
    void write_data_thermo(class_thermodynamics &, const int &iter);
    void write_final_data(class_stats &stats);
//    void write_data_worker_binary(class_worker &);
//    void write_data_master_binary(class_worker &);
    template<typename Derived>
    void write_to_file(const MatrixBase<Derived> &data, string &filename){
        ofstream file(filename,ios::out | ios::trunc);
        file << fixed << showpoint << setprecision(precision);
        string      _coeffSeparator = "	";
        IOFormat fmt(StreamPrecision, DontAlignCols, _coeffSeparator);
        file << data.format(fmt) << endl;
        file.close();
    }
    template<typename Derived>
    void write_to_file(const ArrayBase<Derived> &data, string &filename){
        ofstream file(filename,ios::out | ios::trunc);
        file << fixed << showpoint << setprecision(precision);
        string      _coeffSeparator = "	";
        IOFormat fmt(StreamPrecision, DontAlignCols, _coeffSeparator);;
        file << data.format(fmt) << endl;
        file.close();
    }

};



//void gather_data(class_worker &worker);



#endif //WL_CLASS_PRINT_DATA_H
