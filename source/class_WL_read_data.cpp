//
// Created by david on 9/1/16.
//

#include "class_WL_read_data.h"
#include "class_WL_worker.h"

#ifdef __linux__
#define os 0
#elif _WIN32
#define os 1
#else
#define os 2
#endif

indata::indata(int &id, int &size): world_ID(id), world_size(size) {
    //Create folder for out data storage
    switch (os) {
        case 0:
            folder = "outdata/";
            break;
        case 1:
            folder = "..\\outdata\\";
            break;
        default:
            folder = "outdata/";
            break;
    }
}


void indata::load_random_section(class_worker &worker) {
    int iteration = rn::uniform_integer(0,constants::simulation_reps-1);
    string name_dos    = folder + to_string(iteration) + string("/dos")   + to_string(world_ID) + string(".dat");
    string name_E_bins = folder + to_string(iteration) + string("/E")     + to_string(world_ID) + string(".dat");
    string name_M_bins = folder + to_string(iteration) + string("/M")     + to_string(world_ID) + string(".dat");
    worker.dos  = read_file(name_dos);
    worker.E_bins = read_file(name_E_bins);
    worker.M_bins = read_file(name_M_bins);
}


void indata::load_full(class_worker &worker) {
    string name_dos    = folder + to_string(worker.iteration) + string("/dos.dat");
    string name_E_bins = folder + to_string(worker.iteration) + string("/E.dat");
    string name_M_bins = folder + to_string(worker.iteration) + string("/M.dat");
    worker.dos_total   = read_file(name_dos);
    worker.E_bins_total = read_file(name_E_bins);
    worker.M_bins_total = read_file(name_M_bins);
    cout <<" Size " << worker.E_bins_total.size() << endl;
}

void indata::load_full(class_stats &stats, class_worker &worker) {
    int reps = constants::bootstrap_reps + constants::simulation_reps;
    stats.T.resize(constants::T_num, 1);
    stats.s.resize(constants::T_num, reps);
    stats.f.resize(constants::T_num, reps);
    stats.c.resize(constants::T_num, reps);
    stats.u.resize(constants::T_num, reps);
    stats.x.resize(constants::T_num, reps);
    stats.dos1D.resize(worker.E_bins_total.size(), reps);

    cout << worker.E_bins_total.size() << endl;
    if (world_ID == 0) {
        string name_T, name_s, name_f, name_c, name_u, name_x, name_dos1D;
        name_T = folder + to_string(0) + string("/T.dat");
        stats.T = read_file(name_T);

        for (int i = 0; i < reps; i++) {
            name_s = folder + to_string(i) + string("/s.dat");
            name_f = folder + to_string(i) + string("/f.dat");
            name_c = folder + to_string(i) + string("/c.dat");
            name_u = folder + to_string(i) + string("/u.dat");
            name_x = folder + to_string(i) + string("/x.dat");
            name_dos1D = folder + to_string(i) + string("/dos1D.dat");
            stats.s.col(i) = read_file(name_s);
            stats.f.col(i) = read_file(name_f);
            stats.c.col(i) = read_file(name_c);
            stats.u.col(i) = read_file(name_u);
            stats.x.col(i) = read_file(name_x);
            stats.dos1D.col(i) = read_file(name_dos1D);
        }
    }
    stats.E = worker.E_bins_total;
    stats.M = worker.M_bins_total;

    MPI_Bcast(stats.T.data(),(int)stats.T.size(),MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(stats.s.data(),(int)stats.s.size(),MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(stats.f.data(),(int)stats.f.size(),MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(stats.c.data(),(int)stats.c.size(),MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(stats.u.data(),(int)stats.u.size(),MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(stats.x.data(),(int)stats.x.size(),MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(stats.dos1D.data(),(int)stats.dos1D.size(),MPI_DOUBLE, 0, MPI_COMM_WORLD);

}




MatrixXd indata::read_file(string filename) {
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];
    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    if (!infile.is_open()){cout << "Could not open file with name: " << filename << endl;}
    string line;
    while (! infile.eof()){
        getline(infile, line);
        int temp_cols = 0;
        double num;
        stringstream stream(line);
        string number;
        char* end;

        while(! stream.eof()) {
            stream >> number;
            buff[cols * rows + temp_cols++] = strtod(number.c_str(),&end);
        }
        if (temp_cols == 0) {
            continue;
        }
        if (cols == 0) {
            cols = temp_cols;
        }
        rows++;
    }
    infile.close();

    rows--;
    // Populate matrix with numbers.
    MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result(i, j) = buff[cols * i + j];
        }
    }

    return result;
};

