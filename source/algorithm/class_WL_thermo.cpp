//
// Created by david on 9/1/16.
//

#include "class_WL_thermo.h"

void class_thermodynamics::compute(class_worker &worker) {
    T                        = Eigen::ArrayXd::LinSpaced(T_num, T_min, T_max);
    Eigen::ArrayXd  beta_vec = T.cwiseInverse();
    Eigen::ArrayXXd betaE    = (beta_vec.matrix() * worker.E_bins_total.matrix().transpose()).array();
    double          lambda   = math::nanmaxCoeff(worker.dos_total);

    dos_total1D              = lambda + math::nansum_rowwise((worker.dos_total - lambda).exp()).log();

    Eigen::ArrayXd lambdaT   = math::nanmaxCoeff_rowwise((-betaE).rowwise() + dos_total1D.transpose());
    D                        = (((-betaE).colwise() - lambdaT).rowwise() + dos_total1D.transpose()).exp();

    // Energy-dependent quantities:
    Eigen::ArrayXXd weight1D = D;
    Eigen::ArrayXd  Z        = math::nansum_rowwise(weight1D);
    Eigen::ArrayXd  eAvg     = math::nansum_rowwise(weight1D.rowwise() * worker.E_bins_total.transpose()) / Z;
    Eigen::ArrayXd  eSqAvg   = math::nansum_rowwise(weight1D.rowwise() * worker.E_bins_total.cwiseAbs2().transpose()) / Z;
    u                        = eAvg / constants::N();
    c                        = (beta_vec.cwiseAbs2() * (eSqAvg - eAvg.cwiseAbs2())) / constants::N();
    s                        = (Z.log() + lambdaT + beta_vec * eAvg) / constants::N();
    s                        = s - s.minCoeff();
    f                        = u - T.cwiseProduct(s);

    // Order-parameter-dependent quantities:
    Eigen::ArrayXd weight;
    Eigen::ArrayXd mAvg   = Eigen::ArrayXd::Zero(constants::T_num, 1);
    Eigen::ArrayXd mSqAvg = Eigen::ArrayXd::Zero(constants::T_num, 1);
    Z                     = Eigen::ArrayXd::Zero(constants::T_num, 1);

    lambdaT               = math::nanmaxCoeff_rowwise((-betaE) + math::nanmaxCoeff(worker.dos_total));
    for(int j = 0; j < worker.M_bins_total.size(); j++) {
        weight = math::nansum_rowwise((((-betaE).colwise() - lambdaT).rowwise() + worker.dos_total.col(j).transpose()).exp());
        mAvg += (weight) *fabs(worker.M_bins_total(j));
        mSqAvg += (weight) *worker.M_bins_total(j) * worker.M_bins_total(j);
        Z += (weight);
    }
    mAvg /= Z;
    mSqAvg /= Z;
    m = mAvg / constants::N();
    x = (mSqAvg - mAvg.cwiseAbs2()) / constants::N();
    double ZT;
    P.resize(T_num, worker.M_bins_total.size());
    for(int t = 0; t < T_num; t++) {
        ZT       = math::nansum(math::nansum_rowwise((worker.dos_total.colwise() - (beta_vec(t) * worker.E_bins_total + lambdaT(t))).exp()));
        P.row(t) = math::nansum_colwise((worker.dos_total.colwise() - (beta_vec(t) * worker.E_bins_total + lambdaT(t))).exp()) / ZT;
    }
    int middle = (int) (worker.M_bins_total.size() - 1) / 2;
    F          = P.log().array().colwise() / (-beta_vec) / constants::N();
    for(int t = 0; t < T_num; t++) { F.row(t) -= F(t, middle); }
}
