//
// Created by david on 9/1/16.
//

#include "class_WL_thermo.h"

void class_thermodynamics::compute(class_worker &worker) {
    int i, j;// , j;
    double beta, E, M, weight, ene, eSq, eAvg, eSqAvg, mag, mSq, mAvg, mSqAvg, Z;
    T = ArrayXd::LinSpaced(T_num, T_min, T_max);
    ArrayXd beta_vec = T.cwiseInverse();
    double lambda = math::nanmaxCoeff(worker.dos_total); //For each energy
    dos_total1D = lambda + math::nansum_rowwise((worker.dos_total-lambda).exp()).log();
    ArrayXd lambdaT = math::nanmaxCoeff_rowwise((-beta_vec.matrix()*worker.E_bins_total.matrix().transpose()).array().rowwise() + dos_total1D.transpose() );
    D = (((-beta_vec.matrix()*worker.E_bins_total.matrix().transpose()).array().colwise() - lambdaT).rowwise() + dos_total1D.transpose() ).exp();



    F.resize(T_num,worker.E_bins_total.size());
    for (i = 0; i < T.size(); i++) {
        F.row(i) = (dos_total1D - worker.E_bins_total * beta_vec(i) - lambdaT).exp();
    }
//    cout << "D: " << D << endl;

    u.resize(T_num);
    m.resize(T_num);
    s.resize(T_num);
    f.resize(T_num);
    c.resize(T_num);
    x.resize(T_num);
    for (int t = 0; t < T_num; t++) {
        ene = 0;
        eSq = 0;
        mag = 0;
        mSq = 0;
        Z   = 0;
        beta = 1 / T(t);
        lambda = 0;
        //Find lambda
        for (i = 0; i < worker.E_bins_total.size(); i++) {
            for (j = 0; j < worker.M_bins_total.size(); j++) {
                if(isnan(worker.dos_total(i,j))){continue;}
                lambda = fmax(lambda ,worker.dos_total(i, j) - worker.E_bins_total(i) * beta);
            }
        }


        for (j = 0; j < worker.M_bins_total.size(); j++) {
            //w = 0;
            //w = exp(worker.dos_total(i,j) - beta*E - lambda);  //// g(E)*exp(-E/T) = the weight of this energy
            M = worker.M_bins_total(j);
            for (i = 0; i < worker.E_bins_total.size(); i++) {
                if(isnan(worker.dos_total(i,j))){continue;}
                E = worker.E_bins_total(i);//- worker.E_bins_total[0];
                //w = exp(worker.dos_total(i,j) - beta*E - lambda);  //// g(E)*exp(-E/T) = the weight of this energy
                weight  = exp(worker.dos_total(i, j) - beta*E - lambda);  //// g(E)*exp(-E/T) = the weight of this energy
                ene     += weight*E;
                mag     += weight*fabs(M);
                mSq     += weight*M*M;
                eSq     += weight*E*E;
                Z       += weight;
            }
        }

        eAvg    = ene / Z;
        eSqAvg  = eSq / Z;
        mAvg    = mag / Z;
        mSqAvg  = mSq / Z;

        u(t)    = (eAvg) / constants::N;
        m(t)    = (mAvg) /  constants::N;
        s(t)    = (log(Z) + lambda + beta * eAvg) /  constants::N;
        c(t)    = (beta * beta * (eSqAvg - eAvg*eAvg)) /  constants::N;
        x(t)    = (mSqAvg - mAvg*mAvg)/constants::N;

    }

    s = s - s.minCoeff();
    f = u - T.cwiseProduct(s);
}

long double temperature_to_specific_heat(objective_function &obj_fun, Array<long double, Dynamic,1> &input){
    int i, j;
    double beta, E, weight, ene, eSq, eAvg, eSqAvg, Z, lambda;

    ene = 0;
    eSq = 0;

    Z   = 0;
    beta = 1 / input(0);
    lambda = 0;
    //Find lambda
    for (i = 0; i < obj_fun.aux[1].size(); i++) {
        for (j = 0; j < obj_fun.aux[2].size(); j++) {
            if(isnan(obj_fun.aux[0](i,j))){continue;}
            lambda = fmax(lambda ,obj_fun.aux[0](i, j) - obj_fun.aux[1](i) * beta);
        }
    }

    for (j = 0; j < obj_fun.aux[2].size(); j++) {
        for (i = 0; i < obj_fun.aux[1].size(); i++) {
            if(isnan(obj_fun.aux[0](i,j))){continue;}
            E = obj_fun.aux[1](i);//- obj_fun.aux[1][0];
            weight  = exp(obj_fun.aux[0](i, j) - beta*E - lambda);  //// g(E)*exp(-E/T) = the weight of this energy
            ene     += weight*E;
            eSq     += weight*E*E;
            Z       += weight;
        }
    }

    eAvg    = ene / Z;
    eSqAvg  = eSq / Z;
    return -(beta * beta * (eSqAvg - eAvg*eAvg)) / constants::N;

}


void class_thermodynamics::get_peak(class_worker &worker){
    ArrayXd lower_bound(1);
    ArrayXd upper_bound(1);
    lower_bound << constants::T_min;
    upper_bound << constants::T_max;
    double tolerance = 1e-8;
    objective_function obj_fun(temperature_to_specific_heat,lower_bound, upper_bound, tolerance, worker.dos_total, worker.E_bins_total,
                               worker.M_bins_total);
    minimize(obj_fun);

    double c_peak;
    c_peak = (double) obj_fun.fitness(obj_fun.optimum);
    peak.resize(2);
    peak << obj_fun.optimum(0), -c_peak;


}