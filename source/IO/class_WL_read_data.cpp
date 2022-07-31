//
// Created by david on 9/1/16.
//

#include <algorithm/class_WL_statistics.h>
#include <algorithm/class_WL_teams.h>
#include <algorithm/class_WL_worker.h>
#include <algorithm/nmspc_WL_parallelization.h>
#include <fstream>
#include <general/nmspc_random_numbers.h>
#include <params/nmspc_WL_constants.h>
#include <string>

#include "class_WL_read_data.h"
#ifdef __linux__
    #define os 0
#elif _WIN32
    #define os 1
#else
    #define os 2
#endif

indata::indata() {
    // Get names for indata storage
    switch(os) {
        case 0: folder = "outdata/"; break;
        case 1: folder = "..\\outdata\\"; break;
        default: folder = "outdata/"; break;
    }
}

void indata::load_random_section(class_worker &worker) {
    if(worker.team->is_leader()) {
        int         iteration   = rn::uniform_integer(0, constants::simulation_reps - 1);
        std::string name_dos    = folder + std::to_string(iteration) + std::string("/dos") + std::to_string(worker.team->get_team_id()) + std::string(".dat");
        std::string name_E_bins = folder + std::to_string(iteration) + std::string("/E") + std::to_string(worker.team->get_team_id()) + std::string(".dat");
        std::string name_M_bins = folder + std::to_string(iteration) + std::string("/M") + std::to_string(worker.team->get_team_id()) + std::string(".dat");
        worker.dos              = read_file(name_dos);
        worker.E_bins           = read_file(name_E_bins);
        worker.M_bins           = read_file(name_M_bins);
    }
}
void indata::load_your_section(class_worker &worker) {
    if(worker.team->is_leader()) {
        int         iteration   = worker.iteration;
        std::string name_dos    = folder + std::to_string(iteration) + std::string("/dos") + std::to_string(worker.team->get_team_id()) + std::string(".dat");
        std::string name_E_bins = folder + std::to_string(iteration) + std::string("/E") + std::to_string(worker.team->get_team_id()) + std::string(".dat");
        std::string name_M_bins = folder + std::to_string(iteration) + std::string("/M") + std::to_string(worker.team->get_team_id()) + std::string(".dat");
        worker.dos              = read_file(name_dos);
        worker.E_bins           = read_file(name_E_bins);
        worker.M_bins           = read_file(name_M_bins);
    }
}

void indata::load_full(class_worker &worker) {
    std::string name_dos    = folder + std::to_string(worker.iteration) + std::string("/dos.dat");
    std::string name_E_bins = folder + std::to_string(worker.iteration) + std::string("/E.dat");
    std::string name_M_bins = folder + std::to_string(worker.iteration) + std::string("/M.dat");
    std::cout << "ID " << worker.world_ID << ": reading file: " << name_dos << std::endl;
    worker.dos_total    = read_file(name_dos);
    worker.E_bins_total = read_file(name_E_bins);
    worker.M_bins_total = read_file(name_M_bins);
}

void indata::load_full(class_stats &stats, class_worker &worker) {
    int         reps   = constants::bootstrap_reps + constants::simulation_reps;
    std::string name_T = folder + std::to_string(0) + std::string("/T.dat");
    stats.T            = read_file(name_T);
    stats.E.resize(worker.E_bins_total.size(), reps);
    stats.M.resize(worker.M_bins_total.size(), reps);
    stats.s.resize(constants::T_num, reps);
    stats.f.resize(constants::T_num, reps);
    stats.c.resize(constants::T_num, reps);
    stats.m.resize(constants::T_num, reps);
    stats.u.resize(constants::T_num, reps);
    stats.x.resize(constants::T_num, reps);
    stats.dos1D.resize(worker.E_bins_total.size(), reps);
    if(debug_load_full_thermo) { std::cout << "Finished Resizing" << std::endl; }

    for(int i = 0; i < reps; i++) {
        if(debug_load_full_thermo) { std::cout << "Loading rep " << std::endl; }
        std::string name_E     = folder + std::to_string(i) + std::string("/E.dat");
        std::string name_M     = folder + std::to_string(i) + std::string("/M.dat");
        std::string name_s     = folder + std::to_string(i) + std::string("/s.dat");
        std::string name_f     = folder + std::to_string(i) + std::string("/f.dat");
        std::string name_c     = folder + std::to_string(i) + std::string("/c.dat");
        std::string name_m     = folder + std::to_string(i) + std::string("/m.dat");
        std::string name_u     = folder + std::to_string(i) + std::string("/u.dat");
        std::string name_x     = folder + std::to_string(i) + std::string("/x.dat");
        std::string name_dos1D = folder + std::to_string(i) + std::string("/dos1D.dat");
        std::string name_dos   = folder + std::to_string(i) + std::string("/dos.dat");
        std::string name_D     = folder + std::to_string(i) + std::string("/D.dat");
        std::string name_F     = folder + std::to_string(i) + std::string("/F.dat");
        std::string name_P     = folder + std::to_string(i) + std::string("/P.dat");
        stats.E.col(i)         = read_file(name_E);
        stats.M.col(i)         = read_file(name_M);
        stats.s.col(i)         = read_file(name_s);
        stats.f.col(i)         = read_file(name_f);
        stats.c.col(i)         = read_file(name_c);
        stats.m.col(i)         = read_file(name_m);
        stats.u.col(i)         = read_file(name_u);
        stats.x.col(i)         = read_file(name_x);
        stats.dos1D.col(i)     = read_file(name_dos1D);
        stats.dos.push_back(read_file(name_dos));
        stats.D.push_back(read_file(name_D));
        stats.F.push_back(read_file(name_F));
        stats.P.push_back(read_file(name_P));
    }
}

Eigen::ArrayXXd indata::read_file(std::string filename) {
    std::ifstream            infile;
    std::vector<std::string> lines;
    std::string              line;
    std::string              number;
    infile.open(filename);

    if(!infile.is_open()) {
        std::cout << "Could not open file with name: " << filename << std::endl;
        MPI_Finalize();
        exit(5);
    }
    unsigned long int rows = 0, cols = 0;
    // Load all the lines first and count rows.
    while(getline(infile, line)) {
        lines.push_back(line);
        rows++;
    }
    // Now count the number of elements on the  first line
    std::stringstream stream(lines[0]);
    while(!stream.eof()) {
        stream >> number;
        cols++;
    }
    // Now make your matrix and fill it with the contents of line[].
    Eigen::ArrayXXd result(rows, cols);
    int             j;
    for(unsigned long int i = 0; i < rows; i++) {
        std::stringstream new_stream;
        new_stream << lines[i];
        j = 0;
        while(!new_stream.eof()) {
            new_stream >> number;
            result(static_cast<long>(i), j++) = std::stod(number);
        }
    }
    return result;
}