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


//void gather_data(class_worker &worker){
//    //Gather data into worker 0
//    MatrixXd dos_receive;
//
//    //Resize dos_receive to be able to receive the matrix
//    if (worker.world_ID == 0){
//        dos_receive = worker.dos;
//    }
//
//    for (int i = 1; i < worker.world_size; i++){
//        if(worker.world_ID == 0){
//            MPI_Recv(dos_receive.data(), (int)dos_receive.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
//
//        }else{
//            MPI_Send(worker.dos.data(), (int) worker.dos.size(), MPI_INT, 0, 0, MPI_COMM_WORLD );
//        }
//    }
//}

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
