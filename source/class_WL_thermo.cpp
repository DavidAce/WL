//
// Created by david on 9/1/16.
//

#include "class_WL_thermo.h"

void class_thermodynamics::compute(class_worker &worker) {
    T = ArrayXd::LinSpaced(T_num, T_min, T_max);
    ArrayXd beta_vec = T.cwiseInverse();
    ArrayXXd betaE   = (beta_vec.matrix()*worker.E_bins_total.matrix().transpose()).array();
    double lambda    = math::nanmaxCoeff(worker.dos_total);

    dos_total1D     = lambda + math::nansum_rowwise((worker.dos_total-lambda).exp()).log();

    ArrayXd  lambdaT = math::nanmaxCoeff_rowwise((-betaE).rowwise() + dos_total1D.transpose() );
    D                = (((-betaE).colwise() - lambdaT).rowwise() + dos_total1D.transpose() ).exp();


    //Energy-dependent quantities:
    ArrayXXd weight1D   = D;
    ArrayXd  Z          = math::nansum_rowwise(weight1D);
    ArrayXd eAvg        = math::nansum_rowwise(weight1D.rowwise() * worker.E_bins_total.transpose())/Z;
    ArrayXd eSqAvg      = math::nansum_rowwise(weight1D.rowwise() * worker.E_bins_total.cwiseAbs2().transpose())/Z;
    u                   = eAvg / constants::N;
    c                   = (beta_vec.cwiseAbs2() * (eSqAvg - eAvg.cwiseAbs2())) /  constants::N;
    s                   = (Z.log() + lambdaT + beta_vec * eAvg) /  constants::N;
    s                   = s - s.minCoeff();
    f                   = u - T.cwiseProduct(s);


    //Order-parameter-dependent quantities:
    ArrayXd weight;
    ArrayXd  mAvg   = ArrayXd::Zero(constants::T_num, 1);
    ArrayXd  mSqAvg = ArrayXd::Zero(constants::T_num, 1);
    Z               = ArrayXd::Zero(constants::T_num, 1);

    lambdaT  = math::nanmaxCoeff_rowwise((-betaE) + math::nanmaxCoeff(worker.dos_total));
    for (int j = 0; j < worker.M_bins_total.size(); j++){
        weight   = math::nansum_rowwise((((-betaE).colwise() - lambdaT).rowwise() + worker.dos_total.col(j).transpose()  ).exp());
        mAvg    += (weight)  * fabs(worker.M_bins_total(j));
        mSqAvg  += (weight)  * worker.M_bins_total(j) * worker.M_bins_total(j) ;
        Z       += (weight);
    }
    mAvg   /= Z;
    mSqAvg /= Z;
    m       = mAvg                          /  constants::N;
    x       = (mSqAvg - mAvg.cwiseAbs2())   /  constants::N;
    double ZT;
    P.resize(T_num, worker.M_bins_total.size());
    for (int t = 0; t < T_num ; t++){
        ZT       = math::nansum(math::nansum_rowwise((worker.dos_total.colwise() - (beta_vec(t)*worker.E_bins_total + lambdaT(t)) ).exp()));
        P.row(t) = math::nansum_colwise( ( worker.dos_total.colwise() - (beta_vec(t)*worker.E_bins_total + lambdaT(t)) ).exp()) / ZT;
    }
    int middle = (int)(worker.M_bins_total.size()-1)/2;
    F = P.log().array().colwise() /(-beta_vec)/ constants::N;
    for (int t = 0; t < T_num ; t++){
        F.row(t) -= F(t,middle);
    }

}


long double temperature_to_specific_heat(objective_function &obj_fun, Array<long double, Dynamic,1> &input){
    auto &dos1D     = obj_fun.aux[0];
    auto &E_bins    = obj_fun.aux[1];
    auto &T         = input(0);
    auto beta = (double) (1 / T);
    auto lambda = math::nanmaxCoeff(dos1D - E_bins*beta);
    auto weight = (dos1D - beta*E_bins - lambda).exp();
    auto ene = math::nansum(weight * E_bins);
    auto eSq = math::nansum(weight * E_bins.cwiseAbs2());
    auto Z   = math::nansum(weight);
    auto eAvg    = ene / Z;
    auto eSqAvg  = eSq / Z;
    return -(beta * beta * (eSqAvg - eAvg*eAvg)) / constants::N;
}

void class_thermodynamics::get_c_peak(class_worker &worker){
    ArrayXd lower_bound(1);
    ArrayXd upper_bound(1);
    lower_bound << constants::T_min;
    upper_bound << constants::T_max;
    double tolerance = 1e-8;
    objective_function obj_fun(temperature_to_specific_heat,lower_bound, upper_bound, tolerance, this->dos_total1D, worker.E_bins_total, worker.M_bins_total);
    obj_fun.id      = worker.world_ID;
    obj_fun.threads = 1;
    minimize(obj_fun);

    double temp  = (double) obj_fun.fitness(obj_fun.optimum);
    c_peak.resize(2);
    c_peak << obj_fun.optimum(0), -temp;


}

long double temperature_to_susceptibility(objective_function &obj_fun, Array<long double, Dynamic,1> &input){
    auto &dos_total = obj_fun.aux[0];
    auto &E_bins    = obj_fun.aux[1];
    auto &M_bins    = obj_fun.aux[2];
    auto &T         =  input(0);
    double beta     = (double) (1 / T);
    double  mAvg   = 0;
    double  mSqAvg = 0;
    double  Z      = 0;
    double lambda  = math::nanmaxCoeff((-beta*E_bins + math::nanmaxCoeff(dos_total)));
    for (int j = 0; j < M_bins.size(); j++){
        double weight   = math::nansum((-beta*E_bins + dos_total.col(j) -   lambda   ).exp().eval());
        mAvg            += (weight)  * fabs(M_bins(j));
        mSqAvg          += (weight)  * M_bins(j) * M_bins(j) ;
        Z               += (weight);
    }
    mAvg   /= Z;
    mSqAvg /= Z;
    return      -(mSqAvg - mAvg*mAvg)   /  constants::N;
}


void class_thermodynamics::get_x_peak(class_worker &worker){
    ArrayXd lower_bound(1);
    ArrayXd upper_bound(1);
    lower_bound << constants::T_min;
    upper_bound << constants::T_max;
    double tolerance = 1e-8;
    objective_function obj_fun(temperature_to_susceptibility,lower_bound, upper_bound, tolerance, worker.dos_total, worker.E_bins_total, worker.M_bins_total);
    obj_fun.id      = worker.world_ID;
    obj_fun.threads = 1;

    minimize(obj_fun);
    double temp  = (double) obj_fun.fitness(obj_fun.optimum);
    x_peak.resize(2);
    x_peak << obj_fun.optimum(0), -temp;


}


long double temperature_to_free_energy(objective_function &obj_fun, Array<long double, Dynamic,1> &input){
    auto &dos_total =  obj_fun.aux[0];
    auto &E_bins    =  obj_fun.aux[1];
    auto &M_bins    =  obj_fun.aux[2];

    auto &t         =  input(0);
    double beta     = (double) (1 / t);
    double lambda   = math::nanmaxCoeff((-beta*E_bins + math::nanmaxCoeff(dos_total)));
    double ZT       = math::nansum(math::nansum_rowwise(  (dos_total.colwise() - (beta*E_bins.col(0) + lambda) ).exp()));
    //Important to transpose below, or else P cant be properly initialized!! (Eigen will complain)
    ArrayXd P       = math::nansum_colwise( ( dos_total.colwise() - (beta*E_bins.col(0) + lambda) ).exp()).transpose() / ZT;

    int mid      = (int)((M_bins.size()-1)/2);
    int mid_mid  = (int)(mid/2);
    ArrayXd F       = P.log().array() /(-beta);/// constants::N;
    F              -= F(mid);
    if (F.hasNaN()){
        F.fill(std::numeric_limits<double>::infinity());
    }
    return math::nansum(F.segment(mid-mid_mid, mid).abs());
}


void class_thermodynamics::get_Tc_free_energy(class_worker &worker) {
    ArrayXd lower_bound(1);
    ArrayXd upper_bound(1);
    lower_bound << constants::T_min;
    upper_bound << constants::T_max;
    double tolerance = 1e-8;
    objective_function obj_fun(temperature_to_free_energy, lower_bound, upper_bound, tolerance, worker.dos_total,
                               worker.E_bins_total, worker.M_bins_total);
    obj_fun.id      = worker.world_ID;
    obj_fun.threads = 1;
    minimize(obj_fun);

    Tc_F.resize(1);
    Tc_F(0) = (double) obj_fun.optimum(0);
}
