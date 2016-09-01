//
// Created by david on 9/1/16.
//

#include "class_WL_statistics.h"
#include "class_WL_read_data.h"
class_stats::class_stats(const int &id, const int &size): world_ID(id), world_size(size){
}

void class_stats::load_thermo_data(){
    indata in(world_ID, world_size);
    in.load_full(*this);
}

void class_stats::compute(){

    int t = world_ID;
    double B = constants::simulation_reps + constants::bootstrap_reps;
    s_avg.resize(constants::T_num);
    c_avg.resize(constants::T_num);
    u_avg.resize(constants::T_num);
    f_avg.resize(constants::T_num);
    x_avg.resize(constants::T_num);
    dos1D_avg.resize(constants::T_num);
    s_err.resize(constants::T_num);
    c_err.resize(constants::T_num);
    u_err.resize(constants::T_num);
    f_err.resize(constants::T_num);
    x_err.resize(constants::T_num);
    dos1D_err.resize(constants::T_num);


    while (t < constants::T_num){
        s_avg(t) = s.row(t).mean();
        c_avg(t) = c.row(t).mean();
        u_avg(t) = u.row(t).mean();
        f_avg(t) = f.row(t).mean();
        x_avg(t) = x.row(t).mean();
        dos1D_avg(t) = dos1D.row(t).mean();

        s_err(t) = sqrt(1/(B-1) * math::nansquared(s.row(t).array() - s_avg(t))) ;
//        c_err(t) = sqrt(1/(B-1) * math::nansquared(c.row(t).array() - c_avg(t)) );
//        u_err(t) = sqrt(1/(B-1) * math::nansquared(u.row(t).array() - u_avg(t)) );
//        f_err(t) = sqrt(1/(B-1) * math::nansquared(f.row(t).array() - f_avg(t)) );
//        x_err(t) = sqrt(1/(B-1) * math::nansquared(x.row(t).array() - x_avg(t)) );
//        dos1D_err(t) = sqrt(1/(B-1) * math::nansquared(dos1D.row(t).array() - dos1D_avg(t)) );
        t += world_size ;
    }

}


