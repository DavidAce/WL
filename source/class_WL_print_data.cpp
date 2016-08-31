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
    name_dos = folder + string("/dos") + to_string(world_ID) + string(".dat");
    file_dos.open(name_dos, ofstream::out | ofstream::trunc);
    file_dos << fixed << showpoint << setprecision(10);

    name_E_bins = folder + string("/E") + to_string(world_ID) + string(".dat");
    file_E_bins.open(name_E_bins, ofstream::out | ofstream::trunc);
    file_E_bins << fixed << showpoint << setprecision(10);

    name_M_bins = folder + string("/M") + to_string(world_ID) + string(".dat");
    file_M_bins.open(name_M_bins, ofstream::out | ofstream::trunc);
    file_M_bins << fixed << showpoint << setprecision(10);
    print_to_file(worker.dos, file_dos);
    print_to_file(worker.E_bins, file_E_bins);
    print_to_file(worker.M_bins, file_M_bins);
}


void outdata::write_data_master(class_worker &worker){
    if (worker.world_ID == 0) {
        name_dos = folder + string("/dos.dat");
        file_dos.open(name_dos, ofstream::out | ofstream::trunc);
        file_dos << fixed << showpoint << setprecision(10);

        name_E_bins = folder + string("/E.dat");
        file_E_bins.open(name_E_bins, ofstream::out | ofstream::trunc);
        file_E_bins << fixed << showpoint << setprecision(10);

        name_M_bins = folder + string("/M.dat");
        file_M_bins.open(name_M_bins, ofstream::out | ofstream::trunc);
        file_M_bins << fixed << showpoint << setprecision(10);


        print_to_file(worker.dos_total, file_dos);
        print_to_file(worker.E_bins_total, file_E_bins);
        print_to_file(worker.M_bins_total, file_M_bins);
    }
}

void outdata::create_directories(){
    //Create folder for out data storage
    ostringstream mdir;
    switch (os) {
        case 0:
            mdir << "mkdir -p outdata/" << to_string(iteration) << "/" ;
            folder = "outdata/" + to_string(iteration);
            break;
        case 1:
            mdir << "mkdir ..\\outdata\\" << to_string(iteration) << "\\" ;
            folder = "..\\outdata\\" + to_string(iteration);
            break;
        default:
            mdir << "mkdir -p outdata/" << to_string(iteration) << "/" ;
            folder = "outdata/" + to_string(iteration);
            break;
    }
    string mdir_str = mdir.str();
    if (std::system(mdir_str.c_str())) {}
}

//Constructor
outdata::outdata(int &id, int &iter): world_ID(id), iteration(iter) {
    create_directories();
}
