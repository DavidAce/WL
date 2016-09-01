//
// Created by david on 2016-07-26.
//
#include <iomanip>
#include "class_WL_print_data.h"
#ifdef __linux__
#define os 0
#elif _WIN32
#define os 1
#else
#define os 2
#endif



void outdata::write_data_worker(class_worker &worker) {
    set_folder(worker.iteration);
    string name_dos = folder + string("dos") + to_string(world_ID) + string(".dat");
    string name_E_bins = folder + string("E") + to_string(world_ID) + string(".dat");
    string name_M_bins = folder + string("M") + to_string(world_ID) + string(".dat");
    write_to_file(worker.dos, name_dos);
    write_to_file(worker.E_bins, name_E_bins);
    write_to_file(worker.M_bins, name_M_bins);
}

void outdata::write_data_master(class_worker &worker){
    if (worker.world_ID == 0) {
        set_folder(worker.iteration);
        string name_dos    = folder + string("dos.dat");
        string name_E_bins = folder + string("E.dat");
        string name_M_bins = folder + string("M.dat");
        write_to_file(worker.dos_total, name_dos);
        write_to_file(worker.E_bins_total, name_E_bins);
        write_to_file(worker.M_bins_total, name_M_bins);
    }
}

void outdata::write_data_thermo(class_thermodynamics &thermo, const int &iter){
    iteration = iter;
    set_folder(iteration);
    string name_T     = folder    + string("T.dat");
    string name_s     = folder + string("s.dat");
    string name_c     = folder + string("c.dat");
    string name_u     = folder + string("u.dat");
    string name_f     = folder + string("f.dat");
    string name_x     = folder + string("x.dat");
    string name_dos1D = folder + string("dos1D.dat");
    write_to_file(thermo.T, name_T);
    write_to_file(thermo.s, name_s);
    write_to_file(thermo.c, name_c);
    write_to_file(thermo.u, name_u);
    write_to_file(thermo.f, name_f);
    write_to_file(thermo.x, name_x);
}

void outdata::write_final_data(class_stats &stats){
    folder = "outdata/final/";
    string name_T     = folder    + string("T.dat");
    string name_s     = folder + string("s.dat");
    string name_c     = folder + string("c.dat");
    string name_u     = folder + string("u.dat");
    string name_f     = folder + string("f.dat");
    string name_x     = folder + string("x.dat");
    string name_dos1D = folder + string("dos1D.dat");
    write_to_file(stats.T, name_T);
    write_to_file(stats.s, name_s);
    write_to_file(stats.c, name_c);
    write_to_file(stats.u, name_u);
    write_to_file(stats.f, name_f);
    write_to_file(stats.x, name_x);

    name_s     = folder + string("s_err.dat");
    name_c     = folder + string("c_err.dat");
    name_u     = folder + string("u_err.dat");
    name_f     = folder + string("f_err.dat");
    name_x     = folder + string("x_err.dat");
    name_dos1D = folder + string("dos1D_err.dat");

    write_to_file(stats.s_err, name_s);
    write_to_file(stats.c_err, name_c);
    write_to_file(stats.u_err, name_u);
    write_to_file(stats.f_err, name_f);
    write_to_file(stats.x_err, name_x);
}



void outdata::set_folder(const int &iter){
    iteration = iter;
    //Set folder for out data storage
    switch (os) {
        case 0:
            folder = "outdata/" + to_string(iteration) + string("/");
            break;
        case 1:
            folder = "..\\outdata\\" + to_string(iteration) + "\\";
            break;
        default:
            folder = "outdata/" + to_string(iteration) + string("/");
            break;
    }
}


void outdata::create_folder(const int &iter){
    iteration = iter;
    //Create folder for out data storage
    ostringstream mdir;
    switch (os) {
        case 0:
            mdir << "mkdir -p " << folder;
            break;
        case 1:
            mdir << "mkdir " << folder ;
            break;
        default:
            mdir << "mkdir -p " << folder ;
            break;
    }
    cout << "Creating folder " << folder << endl;
    string mdir_str = mdir.str();
    if (std::system(mdir_str.c_str())) {}
}

void outdata::create_and_set_folder(const int &iter){
    iteration = iter;
    set_folder(iter);
    create_folder(iter);
}

//Constructor
outdata::outdata(const int &id, const int &iter): world_ID(id), iteration(iter) {
    set_folder(iter);
    create_folder(iter);
}

//Constructor (does not set folder! make sure to set it yourself!
outdata::outdata(const int &id): world_ID(id) {

}