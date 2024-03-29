//
// Created by david on 9/1/16.
//

#ifndef WL_CLASS_WL_READ_DATA_H
#define WL_CLASS_WL_READ_DATA_H
#include <Eigen/Core>
#include <iostream>

const int debug_load_full_thermo = 1;
const int debug_load_full_statis = 1;

class indata {
    private:
    std::string folder;

    public:
    indata();
    void            load_random_section(class_worker &worker);
    void            load_your_section(class_worker &worker);
    void            load_full(class_worker &worker);
    void            load_full(class_stats &stats, class_worker &worker);
    Eigen::ArrayXXd read_file(std::string file);
};

#endif // WL_CLASS_WL_READ_DATA_H
