//
// Created by david on 2016-07-26.
//
#include <iomanip>
#include "class_data.h"
#ifdef __linux__
#define os 0
#elif _WIN32
#define os 1
#else
#define os 2
#endif


void outdata::write_data(class_worker &worker){
    print_to_file(worker.dos, file_dos);
    print_to_file(worker.E_bins, file_E_bins);
    print_to_file(worker.M_bins, file_M_bins);

}

void outdata::create_directories(){
    //Create folder for out data storage
    string mdir;
    switch (os) {
        case 0:
            mdir = "mkdir -p outdata";
            folder = "outdata";
            break;
        case 1:
            mdir = "mkdir ..\\outdata";
            folder = "..\\outdata";
            break;
        default:
            mdir = "mkdir -p outdata";
            folder = "outdata";
            break;
    }
    if (std::system(mdir.c_str())) {}
}

outdata::outdata() {
    create_directories();

    name_dos    = folder + "/dos.dat";
    file_dos.open(name_dos.c_str() , ofstream::out | ofstream::trunc);
    file_dos << fixed << showpoint << setprecision(10);

    name_E_bins = folder + "/E.dat";
    file_E_bins.open(name_E_bins.c_str() , ofstream::out | ofstream::trunc);
    file_E_bins << fixed << showpoint << setprecision(10);

    name_M_bins = folder + "/M.dat";
    file_M_bins.open(name_M_bins.c_str() , ofstream::out | ofstream::trunc);
    file_M_bins << fixed << showpoint << setprecision(10);
}

outdata::outdata(int &i): world_ID(i) {
    create_directories();
    name_dos = folder + "/dos" + patch::to_string(world_ID) + ".dat";
    file_dos.open(name_dos.c_str() , ofstream::out | ofstream::trunc);
    file_dos << fixed << showpoint << setprecision(10);

    name_E_bins = folder + "/E" + patch::to_string(world_ID) + ".dat";
    file_E_bins.open(name_E_bins.c_str() , ofstream::out | ofstream::trunc);
    file_E_bins << fixed << showpoint << setprecision(10);

    name_M_bins = folder + "/M" + patch::to_string(world_ID) + ".dat";
    file_M_bins.open(name_M_bins.c_str() , ofstream::out | ofstream::trunc);
    file_M_bins << fixed << showpoint << setprecision(10);
}
