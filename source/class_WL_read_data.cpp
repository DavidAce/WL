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

indata::indata(){
    //Get names for indata storage
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
    string name_dos    = folder + to_string(iteration) + string("/dos")   + to_string(worker.world_ID) + string(".dat");
    string name_E_bins = folder + to_string(iteration) + string("/E")     + to_string(worker.world_ID) + string(".dat");
    string name_M_bins = folder + to_string(iteration) + string("/M")     + to_string(worker.world_ID) + string(".dat");
    worker.dos  = read_file(name_dos);
    worker.E_bins = read_file(name_E_bins);
    worker.M_bins = read_file(name_M_bins);
}


void indata::load_full(class_worker &worker) {
    string name_dos     = folder + to_string(worker.iteration) + string("/dos.dat");
    string name_E_bins  = folder + to_string(worker.iteration) + string("/E.dat");
    string name_M_bins  = folder + to_string(worker.iteration) + string("/M.dat");
    worker.dos_total    = read_file(name_dos);
    worker.E_bins_total = read_file(name_E_bins);
    worker.M_bins_total = read_file(name_M_bins);
}

void indata::load_full(class_stats &stats, class_worker &worker) {
    int reps = constants::bootstrap_reps + constants::simulation_reps;
    string name_T = folder + to_string(0) + string("/T.dat");
    stats.T = read_file(name_T);
    stats.E.resize(worker.E_bins_total.size(), reps);
    stats.M.resize(worker.M_bins_total.size(), reps);
    stats.s.resize(constants::T_num, reps);
    stats.f.resize(constants::T_num, reps);
    stats.c.resize(constants::T_num, reps);
    stats.u.resize(constants::T_num, reps);
    stats.x.resize(constants::T_num, reps);
    stats.dos1D.resize(worker.E_bins_total.size(), reps);
    stats.c_peak.resize(2, reps);
    cout << "Finished Resizing" << endl;
    for (int i = 0; i < reps; i++) {
        cout << "Loading rep " << i << endl;
        string name_E       = folder + to_string(i) + string("/E.dat");
        string name_M       = folder + to_string(i) + string("/M.dat");
        string name_s       = folder + to_string(i) + string("/s.dat");
        string name_f       = folder + to_string(i) + string("/f.dat");
        string name_c       = folder + to_string(i) + string("/c.dat");
        string name_u       = folder + to_string(i) + string("/u.dat");
        string name_x       = folder + to_string(i) + string("/x.dat");
        string name_dos1D   = folder + to_string(i) + string("/dos1D.dat");
        string name_c_peak  = folder + to_string(i) + string("/c_peak.dat");
        string name_dos     = folder + to_string(i) + string("/dos.dat");
        string name_D       = folder + to_string(i) + string("/D.dat");
        string name_F       = folder + to_string(i) + string("/F.dat");
        stats.E.col(i)      = read_file(name_E);
        stats.M.col(i)      = read_file(name_M);
        stats.s.col(i)      = read_file(name_s);
        stats.f.col(i)      = read_file(name_f);
        stats.c.col(i)      = read_file(name_c);
        stats.u.col(i)      = read_file(name_u);
        stats.x.col(i)      = read_file(name_x);
        stats.dos1D.col(i)  = read_file(name_dos1D);
        stats.c_peak.col(i) = read_file(name_c_peak);
        stats.dos.push_back(read_file(name_dos));
        stats.D.push_back(read_file(name_D));
        stats.F.push_back(read_file(name_F));
    }
}

ArrayXXd indata::read_file(string filename) {
    ifstream infile;
    vector<string> lines;
    string line;
    string number;
    infile.open(filename);

    if (!infile.is_open()) {
        cout << "Could not open file with name: " << filename << endl;
        exit(5);
    }
    unsigned long int rows = 0, cols = 0;
    //Load all the lines first and count rows.
    while (getline(infile,line)) {
        lines.push_back(line);
        rows++;
    }
    //Now count the number of elements on the  first line
    stringstream stream(lines[0]);
    while (!stream.eof()) {
        stream >> number;
        cols++;
    }
    //Now make your matrix and fill it with the contents of line[].
    ArrayXXd result(rows, cols);
    for (unsigned long int i = 0; i < rows; i++) {
        stream.clear();
        stream << lines[i];
        for (unsigned long int j = 0; j < cols; j++) {
            stream >> number;
            result(i, j) = std::stod(number);
        }
    }
    return result;
}

