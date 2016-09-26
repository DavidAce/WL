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


void resize_bins(ArrayXXd &dos, ArrayXXd &E_bins, ArrayXXd & M_bins, const int &E_new_size, const int &M_new_size) {
    // This function does rebinning of dos and histograms.
    // If E_set contains more than the default number of bins, then enlarge E_bins, otherwise shrink it!
    // If M_set contains more ----" " ---
    int x, y, i, j;
    double dE, dM, dR, dx, dy, w;

    //Check if we need more bins
    int E_old_size = (int) E_bins.size();
    int M_old_size = (int) M_bins.size();
    double E_max_local = E_bins.maxCoeff();
    double E_min_local = E_bins.minCoeff();
    double M_max_local = M_bins.maxCoeff();
    double M_min_local = M_bins.minCoeff();

    ArrayXXd dos_new        = ArrayXXd::Zero(E_new_size, M_new_size);
    dE              = fabs(E_max_local - E_min_local) / E_new_size;    //New spacing in E_bins
    dM              = fabs(M_max_local - M_min_local) / M_new_size;    //New spacing in M_bins
    dR              = sqrt(dE * dE + dM * dM);                         //Diagonal distance ?
    ArrayXXd E_old = E_bins;
    ArrayXXd M_old = M_bins;
    E_bins          = ArrayXd::LinSpaced(E_new_size, E_min_local, E_max_local);
    M_bins          = ArrayXd::LinSpaced(M_new_size, M_min_local, M_max_local);
    ArrayXXd weight = ArrayXXd::Zero(E_new_size, M_new_size);
    //Coarsen the histogram and dos.
    for (y = 0; y < M_old_size; y++) {
        for (x = 0; x < E_old_size; x++) {
            for (j = 0; j < M_new_size; j++) {
                for (i = 0; i < E_new_size; i++) {
                    dx = fabs(E_bins(i) - E_old(x));   //Distance between old and new bins
                    dy = fabs(M_bins(j) - M_old(y));   //Distance between old and new bins
                    //Distance between old and new bins should not exceed dE/2
                    //This is so to avoid double counting
                    if (dx >= dE/2 || dy >= dM/2) { continue; }
//                      w          = fabs(1.0 - 4*dx*dy/dE/dM);
                    w          = fabs(1.0 - sqrt(dx * dx + dy * dy) / dR);

                    weight(i,j)              += w;
                    dos_new(i,j)             +=  dos(x,y) * w;
//                    dos(x,y)                 -= 0 ; //This entry has been used
                }
            }

        }
    }
    //We have now inserted all old entries to the new dos and histogram, and we only need to divide by the weight.
    for (j = 0; j < M_new_size; j++) {
        for (i = 0; i < E_new_size; i++) {
            dos_new(i,j)        = weight(i,j) > 0  ? dos_new(i,j)/weight(i,j) : 0;

        }
    }
    dos             = dos_new;
}


double temperature_to_specific_heat(objective_function &obj_fun, ArrayXd &input){
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
    ArrayXd initial_conditions(3);
    initial_conditions.fill(0);
    lower_bound << constants::T_min;
    upper_bound << constants::T_max;
    double tolerance = 1e-8;
    objective_function obj_fun(temperature_to_specific_heat,lower_bound, upper_bound, tolerance,initial_conditions, this->dos_total1D, worker.E_bins_total, worker.M_bins_total);
    obj_fun.id      = worker.world_ID;
    obj_fun.name    = "c";
    minimize(obj_fun);

    double temp  = (double) obj_fun.fitness(obj_fun.optimum);
    c_peak.resize(2);
    c_peak << obj_fun.optimum(0), -temp;
    cout << "ID: " << worker.world_ID << " Tc_c = " << c_peak(0) << endl;

}

double temperature_to_susceptibility(objective_function &obj_fun, ArrayXd &input){
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
    ArrayXd initial_conditions(3);
    initial_conditions.fill(0);
    objective_function obj_fun(temperature_to_susceptibility,lower_bound, upper_bound, tolerance,initial_conditions, worker.dos_total, worker.E_bins_total, worker.M_bins_total);
    obj_fun.id      = worker.world_ID;
    obj_fun.name    = "x";

    minimize(obj_fun);
    double temp  = (double) obj_fun.fitness(obj_fun.optimum);
    x_peak.resize(2);
    x_peak << obj_fun.optimum(0), -temp;
    cout << "ID: " << worker.world_ID << " Tc_x = " << x_peak(0) << endl;


}


double temperature_to_free_energy(objective_function &obj_fun, ArrayXd &input){
    auto &dos_total =  obj_fun.aux[0];
    auto &E_bins    =  obj_fun.aux[1];
    auto &M_bins    =  obj_fun.aux[2];
    auto &T         = input(0);
    double lambda   = (-E_bins/ T + dos_total.maxCoeff()).maxCoeff();
    //Important to transpose below, or else P can't be properly initialized!! (Eigen will complain)
    ArrayXd F       = (((dos_total.colwise() - (E_bins.col(0)/ T + lambda) ).exp()).colwise().sum().transpose());
//    F /= F.sum();
    int mid      = (int)((M_bins.size()-1)/2);
    int mid_mid  = (int)(mid/2);
    F            = F.log().array() * (-T);
    F              -= F(mid);
    double result = (F.segment(mid - mid_mid, mid).abs()).sum() + 1/fabs(F(0));
    if (std::isinf(result)){result = 1e6;}
    if (std::isnan(result)){result = 1e6;}
    return result;

//    return (F.segment(mid-mid_mid, mid).abs()).sum();
}



void class_thermodynamics::get_Tc_free_energy(class_worker &worker) {
    ArrayXd lower_bound(1);
    ArrayXd upper_bound(1);
    ArrayXd initial_conditions(3);
    initial_conditions.fill(0);
    lower_bound << constants::T_min;
    upper_bound << constants::T_max;
    double tolerance = 1e-6;
    ArrayXXd dos_no_nan = math::NaN_to_Zero(worker.dos_total);
    objective_function obj_fun(temperature_to_free_energy, lower_bound, upper_bound, tolerance,initial_conditions, dos_no_nan,
                               worker.E_bins_total, worker.M_bins_total);
    obj_fun.id      = worker.world_ID;
    obj_fun.name    = "F";
    minimize(obj_fun);
    Tc_F.resize(1);
    Tc_F(0) = (double) obj_fun.optimum(0);
    cout << "ID: " << worker.world_ID << " Tc_F = " << Tc_F(0) << endl;
}


double temperature_to_canonical_distribution(objective_function &obj_fun, ArrayXd &input){
    auto &dos1D     = obj_fun.aux[0];
    auto &E_bins    = obj_fun.aux[1];
    auto &T         = input(0);
    auto beta = (double) (1 / T);
    auto lambda = math::nanmaxCoeff(dos1D - E_bins*beta);
    auto weight = (dos1D - beta*E_bins - lambda).exp();
    auto Z   = math::nansum(weight);
    return (weight/Z).maxCoeff();


}


void class_thermodynamics::get_Tc_canonical_distribution(class_worker &worker) {
    ArrayXd lower_bound(1);
    ArrayXd upper_bound(1);
    ArrayXd initial_conditions(3);
    initial_conditions.fill(0);
    lower_bound << constants::T_min;
    upper_bound << constants::T_max;
    double tolerance = 1e-6;
    ArrayXXd dos_no_nan = math::NaN_to_Zero(worker.dos_total);
    objective_function obj_fun(temperature_to_canonical_distribution, lower_bound,  upper_bound, tolerance,initial_conditions, this->dos_total1D, worker.E_bins_total, worker.M_bins_total);
    obj_fun.id      = worker.world_ID;
    obj_fun.name    = "D";
    minimize(obj_fun);
    Tc_D.resize(1);
    Tc_D(0) = (double) obj_fun.optimum(0);
    cout << "ID: " << worker.world_ID << " Tc_D = " << Tc_D(0) << endl;
}
