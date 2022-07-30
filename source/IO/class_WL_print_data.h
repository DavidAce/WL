//
// Created by david on 2016-07-26.
//

#ifndef WL_CLASS_PRINT_DATA_H
#define WL_CLASS_PRINT_DATA_H
#include "IO/std_filesystem.h"
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>

class class_thermodynamics;
class class_stats;
class class_worker;

class outdata {
    private:
    std::string folder;
    int         iteration = 0;
    int         precision = 12;

    public:
    outdata() = default; // Default constructor (does not set folder! make sure to
                         // set it yourself!

    void create_folder(std::string folder_name);
    void create_set_folder(std::string folder_name);
    void create_iteration_folder_commander(class_worker &, int iter);
    void create_iteration_folder_worker(const int iter);
    void set_foldername_to_iteration(const int iter);
    int  mkdir_p(const char *path);
    // File streams

    void write_data_team_leader(class_worker &);
    void write_data_commander(class_worker &);
    void write_data_thermo(class_thermodynamics &, const int iter);
    void write_final_data(class_worker &worker, class_stats &stats);
    void write_sample(class_worker &);

    template<typename Derived>
    void write_to_file(const Eigen::EigenBase<Derived> &data, const std::string &filename) {
        auto filepath = fs::absolute(fs::path(filename));

        if(not fs::is_directory(filepath.parent_path())) {
            std::cout << "Creating directory: " << filepath.parent_path() << std::endl;
            fs::create_directories(filepath.parent_path());
        }
        std::ofstream file(filepath, std::ios::out | std::ios::trunc);
        if(not file.is_open()) throw std::runtime_error("Failed to open file: " + filepath.string());
        file << std::fixed << std::showpoint << std::setprecision(precision);
        std::string     _coeffSeparator = "	";
        Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, _coeffSeparator);
        file << data.derived().format(fmt) << std::endl;
        file.close();
    }

    void write_to_file(const double data, const std::string &filename) {
        auto filepath = fs::absolute(fs::path(filename));
        if(not fs::is_directory(filepath.parent_path())) {
            std::cout << "Creating directory: " << filepath.parent_path() << std::endl;
            fs::create_directories(filepath.parent_path());
        }
        std::ofstream file(filepath, std::ios::out | std::ios::trunc);
        file << std::fixed << std::showpoint << std::setprecision(precision) << data;
        file.close();
    }
};

#endif // WL_CLASS_PRINT_DATA_H
