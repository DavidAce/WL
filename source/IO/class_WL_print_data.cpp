//
// Created by david on 2016-07-26.
//
#include "class_WL_print_data.h"
#include "algorithm/class_WL_statistics.h"
#include "algorithm/class_WL_teams.h"
#include "algorithm/class_WL_thermo.h"
#include "algorithm/class_WL_worker.h"
#include <iomanip>
#ifdef __linux__
    #define os 0
#elif _WIN32
    #define os 1
#else
    #define os 2
#endif

void outdata::write_data_team_leader(class_worker &worker) {
    if(worker.team->is_leader()) {
        set_foldername_to_iteration(worker.iteration);
        string name_dos    = folder +std::string("dos") + to_string(worker.team->get_team_id()) +std::string(".dat");
        string name_E_bins = folder +std::string("E") + to_string(worker.team->get_team_id()) +std::string(".dat");
        string name_M_bins = folder +std::string("M") + to_string(worker.team->get_team_id()) +std::string(".dat");
        for(int id = 0; id < worker.team->get_team_size(); id++) {
            if(id == worker.team->get_team_id()) {
                write_to_file(worker.dos, name_dos);
                write_to_file(worker.E_bins, name_E_bins);
                write_to_file(worker.M_bins, name_M_bins);
            }
            MPI_Barrier(worker.team->get_MPI_COMM_LEAD());
        }
    }
}

void outdata::write_data_commander(class_worker &worker) {
    if(worker.team->is_commander()) {
        set_foldername_to_iteration(worker.iteration);
        string name_dos    = folder +std::string("dos.dat");
        string name_E_bins = folder +std::string("E.dat");
        string name_M_bins = folder +std::string("M.dat");
        write_to_file(worker.dos_total, name_dos);
        write_to_file(worker.E_bins_total, name_E_bins);
        write_to_file(worker.M_bins_total, name_M_bins);
    }
}

void outdata::write_data_thermo(class_thermodynamics &thermo, const int iter) {
    create_iteration_folder_worker(iter);
    string name_T     = folder +std::string("T.dat");
    string name_s     = folder +std::string("s.dat");
    string name_c     = folder +std::string("c.dat");
    string name_m     = folder +std::string("m.dat");
    string name_u     = folder +std::string("u.dat");
    string name_f     = folder +std::string("f.dat");
    string name_x     = folder +std::string("x.dat");
    string name_D     = folder +std::string("D.dat");
    string name_F     = folder +std::string("F.dat");
    string name_P     = folder +std::string("P.dat");
    string name_dos1D = folder +std::string("dos1D.dat");
    write_to_file(thermo.T, name_T);
    write_to_file(thermo.s, name_s);
    write_to_file(thermo.c, name_c);
    write_to_file(thermo.m, name_m);
    write_to_file(thermo.u, name_u);
    write_to_file(thermo.f, name_f);
    write_to_file(thermo.x, name_x);
    write_to_file(thermo.D, name_D);
    write_to_file(thermo.F, name_F);
    write_to_file(thermo.P, name_P);
    write_to_file(thermo.dos_total1D, name_dos1D);
}

void outdata::write_final_data(class_worker &worker, class_stats &stats) {
    if(worker.team->is_commander()) {
        folder = "outdata/final/";
        create_folder(folder);
        string name_E     = folder +std::string("E.dat");
        string name_M     = folder +std::string("M.dat");
        string name_T     = folder +std::string("T.dat");
        string name_s     = folder +std::string("s.dat");
        string name_c     = folder +std::string("c.dat");
        string name_m     = folder +std::string("m.dat");
        string name_u     = folder +std::string("u.dat");
        string name_f     = folder +std::string("f.dat");
        string name_x     = folder +std::string("x.dat");
        string name_dos1D = folder +std::string("dos1D.dat");
        string name_dos   = folder +std::string("dos.dat");
        string name_D     = folder +std::string("D.dat");
        string name_F     = folder +std::string("F.dat");
        string name_P     = folder +std::string("P.dat");

        write_to_file(stats.E_avg, name_E);
        write_to_file(stats.M_avg, name_M);
        write_to_file(stats.T, name_T);
        write_to_file(stats.s_avg, name_s);
        write_to_file(stats.c_avg, name_c);
        write_to_file(stats.m_avg, name_m);
        write_to_file(stats.u_avg, name_u);
        write_to_file(stats.f_avg, name_f);
        write_to_file(stats.x_avg, name_x);
        write_to_file(stats.dos1D_avg, name_dos1D);
        write_to_file(stats.dos_avg, name_dos);
        write_to_file(stats.D_avg, name_D);
        write_to_file(stats.F_avg, name_F);
        write_to_file(stats.P_avg, name_P);

        name_s     = folder +std::string("s_err.dat");
        name_c     = folder +std::string("c_err.dat");
        name_m     = folder +std::string("m_err.dat");
        name_u     = folder +std::string("u_err.dat");
        name_f     = folder +std::string("f_err.dat");
        name_x     = folder +std::string("x_err.dat");
        name_dos1D = folder +std::string("dos1D_err.dat");
        name_dos   = folder +std::string("dos_err.dat");
        name_D     = folder +std::string("D_err.dat");
        name_F     = folder +std::string("F_err.dat");
        name_P     = folder +std::string("P_err.dat");

        write_to_file(stats.s_err, name_s);
        write_to_file(stats.c_err, name_c);
        write_to_file(stats.m_err, name_m);
        write_to_file(stats.u_err, name_u);
        write_to_file(stats.f_err, name_f);
        write_to_file(stats.x_err, name_x);
        write_to_file(stats.dos1D_err, name_dos1D);
        write_to_file(stats.dos_err, name_dos);
        write_to_file(stats.D_err, name_D);
        write_to_file(stats.F_err, name_F);
        write_to_file(stats.P_err, name_P);
    }
}

void outdata::write_sample(class_worker &worker) {
    std::string name_sample_lattice = folder + "lattice_" + to_string(iteration) + ".dat";
    std::string name_sample_energy  = folder + "energy_" + to_string(iteration) + ".dat";
    std::string name_sample_magnet  = folder + "magnet_" + to_string(iteration) + ".dat";
    write_to_file(worker.model.lattice, name_sample_lattice);
    write_to_file(worker.model.get_E(), name_sample_energy);
    write_to_file(worker.model.get_M(), name_sample_magnet);
    iteration++;
}

void outdata::set_foldername_to_iteration(const int iter) {
    iteration = iter;
    // Set folder for out data storage
    switch(os) {
        case 0: folder = "outdata/" + to_string(iteration) +std::string("/"); break;
        case 1: folder = "..\\outdata\\" + to_string(iteration) + "\\"; break;
        default: folder = "outdata/" + to_string(iteration) +std::string("/"); break;
    }
}

int outdata::mkdir_p(const char *path) {
    /* Adapted from http://stackoverflow.com/a/2336245/119527 */
    const size_t len = strlen(path);
    char         _path[PATH_MAX];
    char        *p;

    errno = 0;

    /* Copy string so its mutable */
    if(len > sizeof(_path) - 1) {
        errno = ENAMETOOLONG;
        return -1;
    }
    strcpy(_path, path);

    /* Iterate the string */
    for(p = _path + 1; *p; p++) {
        if(*p == '/') {
            /* Temporarily truncate */
            *p = '\0';

            if(mkdir(_path, S_IRWXU) != 0) {
                if(errno != EEXIST) return -1;
            }

            *p = '/';
        }
    }

    if(mkdir(_path, S_IRWXU) != 0) {
        if(errno != EEXIST) return -1;
    }

    return 0;
}

void outdata::create_folder(string folder_name) {
    if(mkdir_p(folder_name.c_str()) == 0) {
        cout << "Set folder: " << folder_name << endl;

    } else {
        cout << "Failed to set folder: " << folder_name << endl;
    }
}

void outdata::create_set_folder(string folder_name) {
    folder = folder_name;
    create_folder(folder);
}

void outdata::create_iteration_folder_commander(class_worker &worker, int iter) {
    if(worker.team->is_commander()) {
        set_foldername_to_iteration(iter);
        create_folder(folder);
    }
}

void outdata::create_iteration_folder_worker(const int iter) {
    set_foldername_to_iteration(iter);
    create_folder(folder);
}
