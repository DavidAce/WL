//
// Created by david on 2016-08-05.
//

#include "nmspc_math_algorithms.h"
#include "nmspc_WL_constants.h"

namespace math {

    int mod2(const int &a, const int &b) {
        if (b < 0) {
            return mod(-a, -b);
        }
        int ret = a % b;
        if (ret < 0) {
            ret += b;
        }
        return ret;
    }

    double volume(const MatrixXd &dos, const MatrixXd &E, const MatrixXd &M) {
        double vol = 0;
        switch (constants::rw_dims){
            case 1:
                for (int j = 0; j < M.size(); j++) {
                    for (int i = 0; i < E.size() - 1; i++) {
                        if (dos(i, j) == 0 || dos(i,j) == std::numeric_limits<double>::quiet_NaN()) { continue; }
                        vol += (E(i + 1) - E(i)) * dos(i, j);
                    }
                }
                break;
            case 2:
                for (int j = 0; j < M.size() - 1; j++) {
                    for (int i = 0; i < E.size() - 1; i++) {
                        if (dos(i, j) == 0 || dos(i,j) == std::numeric_limits<double>::quiet_NaN()) { continue; }
                        vol += (E(i + 1) - E(i)) * (M(j + 1) - M(j)) * dos(i, j);
                    }
                }
                break;
            default:
                std::cout << "Wrong rw-dimension constants::rw_dims " << std::endl;
                exit(5);

        }
        return vol;
    }


    int volume_idx(const MatrixXd &dos, const MatrixXd &E, const MatrixXd &M, const double &vol_limit) {
        double vol = 0;
        switch(constants::rw_dims){
            case 1:
                for (int i = 0; i < E.size() - 1; i++) {
                    for (int j = 0; j < M.size(); j++) {
                        if (dos(i, j) == 0 || dos(i,j) == std::numeric_limits<double>::quiet_NaN()) { continue; }
                        vol += (E(i + 1) - E(i)) * dos(i, j);
                    }
                    if (vol >= vol_limit){
                        return i;
                    }
                }
                break;
            case 2:
                for (int i = 0; i < E.size() - 1; i++) {
                    for (int j = 0; j < M.size() - 1; j++) {
                        if (dos(i, j) == 0 || dos(i,j) == std::numeric_limits<double>::quiet_NaN()) { continue; }
                        vol += (E(i + 1) - E(i)) * (M(j + 1) - M(j)) * dos(i, j);
                    }
                    if (vol >= vol_limit){
                        return i;
                    }
                }
                break;
            default:
                std::cout << "Wrong rw-dimension constants::rw_dims " << std::endl;
                exit(5);
        }

        return (int)E.size() - 1;
    }

    Vector3d gradient_vector(const MatrixXd &dos, const VectorXd &E, const VectorXd &M, const int &i, const int &j){
        Vector3d v_up, v_rt, v_dn, v_lf; //Vectors connecting adjacent 3 orthogonal points on DOS
        switch(constants::rw_dims){
            case 1:
                v_dn << E(i + 1) - E(i), 0, dos(i + 1, j) - dos(i, j);
                return v_dn.normalized();
            case 2:
                if (j <= 0 || j >= M.size()-1){v_up << 0,0,0; return v_up;} //Detect if on edge
                v_up << E(i - 1) - E(i), 0, dos(i - 1, j) - dos(i, j);
                v_dn << E(i + 1) - E(i), 0, dos(i + 1, j) - dos(i, j);
                v_rt << 0, M(j + 1) - M(j), dos(i, j + 1) - dos(i, j);
                v_lf << 0, M(j - 1) - M(j), dos(i, j - 1) - dos(i, j);
                return (v_up.cross(v_rt) + v_rt.cross(v_dn) + v_dn.cross(v_lf) + v_lf.cross(v_dn)).normalized();

        }
    }

    int  find_matching_slope(const MatrixXd &dos1, const MatrixXd &dos2,
                             const VectorXd &E1  , const VectorXd &E2,
                             const VectorXd &M1  , const VectorXd &M2){
        Vector3d u1, u2;                //Vectors connecting adjacent 3 orthogonal points on DOS
        //Vector3d v_up, v_rt, v_dn, v_lf; //Vectors connecting adjacent 3 orthogonal points on DOS
        VectorXd sum(E1.size());
        sum.fill(0);
        int col = 1, num;
        int E_merge_idx;    //Index of merging point

        if (dos1.rows() < 3){
            std::cout << "Rows too few to merge!!" << std::endl;
            exit(35);
        }
        int x, y; //Coordinates closest to i, j on dos of neighbor above.
        for (int i = 1; i < E1.size() - 1; i++) {
            num = 0;
            for (int j = 0; j < M1.size(); j++) {
                if (E1(i) <= E2.minCoeff()) { continue; }
                if (E1(i) >= E2.maxCoeff()) { continue; }
                if (M1(j) >  M2.maxCoeff()) { continue; }
                if (M1(j) <  M2.minCoeff()) { continue; }
                x = math::binary_search(E2.data(), E1(i), E2.size());
                y = math::binary_search(M2.data(), M1(j), M2.size());
                if (x <= 0 || x >= E2.size() - 1) { continue; }

                u1 = gradient_vector(dos1, E1, M1, i, j);
                u2 = gradient_vector(dos2, E2, M2, x, y);
                sum(col) += u1.dot(u2);
                num++;
            }
            sum(col) /= std::fmax(1, num);
            col++;
        }
        sum.maxCoeff(&E_merge_idx);
        E_merge_idx = sum.sum() == 0 ? (int)E1.size() - 1 : E_merge_idx;
        return E_merge_idx;
    }


}
