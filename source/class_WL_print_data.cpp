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
    set_foldername_to_iteration(worker.iteration);
    string name_dos     = folder + string("dos") + to_string(worker.world_ID) + string(".dat");
    string name_E_bins  = folder + string("E") + to_string(worker.world_ID) + string(".dat");
    string name_M_bins  = folder + string("M") + to_string(worker.world_ID) + string(".dat");
    write_to_file(worker.dos, name_dos);
    write_to_file(worker.E_bins, name_E_bins);
    write_to_file(worker.M_bins, name_M_bins);
}

void outdata::write_data_master(class_worker &worker){
    if (worker.world_ID == 0) {
        set_foldername_to_iteration(worker.iteration);
        string name_dos    = folder + string("dos.dat");
        string name_E_bins = folder + string("E.dat");
        string name_M_bins = folder + string("M.dat");
        write_to_file(worker.dos_total, name_dos);
        write_to_file(worker.E_bins_total, name_E_bins);
        write_to_file(worker.M_bins_total, name_M_bins);
    }
}

void outdata::write_data_thermo(class_thermodynamics &thermo, const int &iter){
    create_iteration_folder_worker(iter);
    string name_T       = folder + string("T.dat");
    string name_s       = folder + string("s.dat");
    string name_c       = folder + string("c.dat");
    string name_m       = folder + string("m.dat");
    string name_u       = folder + string("u.dat");
    string name_f       = folder + string("f.dat");
    string name_x       = folder + string("x.dat");
    string name_dos1D   = folder + string("dos1D.dat");
    string name_c_peak  = folder + string("c_peak.dat");
    string name_x_peak  = folder + string("x_peak.dat");
    string name_Tc_F    = folder + string("Tc_F.dat");
    string name_D       = folder + string("D.dat");
    string name_F       = folder + string("F.dat");
    string name_P       = folder + string("P.dat");
    write_to_file(thermo.T, name_T);
    write_to_file(thermo.s, name_s);
    write_to_file(thermo.c, name_c);
    write_to_file(thermo.m, name_m);
    write_to_file(thermo.u, name_u);
    write_to_file(thermo.f, name_f);
    write_to_file(thermo.x, name_x);
    write_to_file(thermo.dos_total1D, name_dos1D);
    write_to_file(thermo.c_peak, name_c_peak);
    write_to_file(thermo.x_peak, name_x_peak);
    write_to_file(thermo.Tc_F, name_Tc_F);
    write_to_file(thermo.D, name_D);
    write_to_file(thermo.F, name_F);
    write_to_file(thermo.P, name_P);
}

void outdata::write_final_data(class_stats &stats, const int &id){
    if (id == 0) {
        folder = "outdata/final/";
        create_folder(folder);
        string name_E       = folder + string("E.dat");
        string name_M       = folder + string("M.dat");
        string name_T       = folder + string("T.dat");
        string name_s       = folder + string("s.dat");
        string name_c       = folder + string("c.dat");
        string name_m       = folder + string("m.dat");
        string name_u       = folder + string("u.dat");
        string name_f       = folder + string("f.dat");
        string name_x       = folder + string("x.dat");
        string name_dos1D   = folder + string("dos1D.dat");
        string name_c_peak  = folder + string("c_peak.dat");
        string name_x_peak  = folder + string("x_peak.dat");
        string name_Tc_F    = folder + string("Tc_F.dat");
        string name_dos     = folder + string("dos.dat");
        string name_D       = folder + string("D.dat");
        string name_F       = folder + string("F.dat");
        string name_P       = folder + string("P.dat");

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
        write_to_file(stats.c_peak_avg, name_c_peak);
        write_to_file(stats.x_peak_avg, name_x_peak);
        write_to_file(stats.Tc_F_avg, name_Tc_F);
        write_to_file(stats.dos_avg, name_dos);
        write_to_file(stats.D_avg, name_D);
        write_to_file(stats.F_avg, name_F);
        write_to_file(stats.P_avg, name_P);

        name_s          = folder + string("s_err.dat");
        name_c          = folder + string("c_err.dat");
        name_m          = folder + string("m_err.dat");
        name_u          = folder + string("u_err.dat");
        name_f          = folder + string("f_err.dat");
        name_x          = folder + string("x_err.dat");
        name_dos1D      = folder + string("dos1D_err.dat");
        name_c_peak     = folder + string("c_peak_err.dat");
        name_x_peak     = folder + string("x_peak_err.dat");
        name_Tc_F       = folder + string("Tc_F_err.dat");
        name_dos        = folder + string("dos_err.dat");
        name_D          = folder + string("D_err.dat");
        name_F          = folder + string("F_err.dat");
        name_P          = folder + string("P_err.dat");

        write_to_file(stats.s_err, name_s);
        write_to_file(stats.c_err, name_c);
        write_to_file(stats.m_err, name_m);
        write_to_file(stats.u_err, name_u);
        write_to_file(stats.f_err, name_f);
        write_to_file(stats.x_err, name_x);
        write_to_file(stats.dos1D_err, name_dos1D);
        write_to_file(stats.c_peak_err, name_c_peak);
        write_to_file(stats.x_peak_err, name_x_peak);
        write_to_file(stats.Tc_F_err, name_Tc_F);
        write_to_file(stats.dos_err, name_dos);
        write_to_file(stats.D_err, name_D);
        write_to_file(stats.F_err, name_F);
        write_to_file(stats.P_err, name_P);
    }
}



void outdata::set_foldername_to_iteration(const int &iter){
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

int outdata::mkdir_p(const char *path)
{
    /* Adapted from http://stackoverflow.com/a/2336245/119527 */
    const size_t len = strlen(path);
    char _path[PATH_MAX];
    char *p;

    errno = 0;

    /* Copy string so its mutable */
    if (len > sizeof(_path)-1) {
        errno = ENAMETOOLONG;
        return -1;
    }
    strcpy(_path, path);

    /* Iterate the string */
    for (p = _path + 1; *p; p++) {
        if (*p == '/') {
            /* Temporarily truncate */
            *p = '\0';

            if (mkdir(_path, S_IRWXU) != 0) {
                if (errno != EEXIST)
                    return -1;
            }

            *p = '/';
        }
    }

    if (mkdir(_path, S_IRWXU) != 0) {
        if (errno != EEXIST)
            return -1;
    }

    return 0;
}

void outdata::create_folder(string folder_name){

    if (mkdir_p(folder_name.c_str()) == 0)
    {
        cout << "Set folder: " << folder_name << endl;

    }else{
        cout << "Failed to set folder: " << folder_name << endl;
    }
}

void outdata::create_iteration_folder_master(const int &iter, const int &id){
    if(id == 0){
        set_foldername_to_iteration(iter);
        create_folder(folder);
    }
}

void outdata::create_iteration_folder_worker(const int &iter){
    set_foldername_to_iteration(iter);
    create_folder(folder);
}


//Default constructor (does not set folder! make sure to set it yourself!
outdata::outdata() {

}