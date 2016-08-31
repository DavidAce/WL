//
// Created by david on 2016-07-26.
//

#ifndef WL_CLASS_PRINT_DATA_H
#define WL_CLASS_PRINT_DATA_H
#include <Eigen/Core>
#include <Eigen/Dense>


#include <iostream>
#include "class_WL_worker.h"
using namespace std;
using namespace Eigen;

namespace patch{
    template < typename T > std::string to_string( const T& n ){
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

class outdata {
private:
    string      folder;
    string      name_dos;
    string      name_E_bins;
    string      name_M_bins;

    ofstream    file_dos;
    ofstream    file_E_bins;
    ofstream    file_M_bins;
    int   world_ID;
    int   iteration;

public:
    outdata(int &, int &);
    void create_directories();
    //File streams

    void write_data_worker(class_worker &);
    void write_data_master(class_worker &);
    template<typename Derived>
    void print_to_file(const Eigen::MatrixBase<Derived>& data, ofstream &file){
        string      _coeffSeparator = "	";
        IOFormat fmt(StreamPrecision, DontAlignCols, _coeffSeparator);;
        file << data.format(fmt) << endl;
        file.close();
    }

};



//void gather_data(class_worker &worker);



#endif //WL_CLASS_PRINT_DATA_H
