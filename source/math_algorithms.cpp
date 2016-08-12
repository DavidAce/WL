//
// Created by david on 2016-08-05.
//

#include "math_algorithms.h"
namespace math {
//    int find_max_idx(const Ref<const VectorXd>  &list) {
//        //int find_max_idx(const Ref<const Array<double,Dynamic,1>> &list){
//        int i, idx = 0;
//        double peak = list(0);
//        for (i = 0; i < list.size(); i++) {
//            if (list(i) > peak) {
//                peak = list(i);
//                idx = i;
//            }
//        }
//        return idx;
//    }

//    int find_min_positive(const MatrixXi &H) {
//        int min = 1000000000;
//        for (int j = 0; j < H.cols(); j++) {
//            for (int i = 0; i < H.rows(); i++) {
//                if (H(i, j) < min) {
//                    if (H(i, j) > 0) {
//                        min = H(i, j);
//                    }
//                }
//            }
//        }
//        return min;
//    }
//
//    double find_min_positive(const MatrixXd &M) {
//        double min = 1000000000;
//        for (int j = 0; j < M.cols(); j++) {
//            for (int i = 0; i < M.rows(); i++) {
//                if (M(i, j) < min) {
//                    if (M(i, j) > 0) {
//                        min = M(i, j);
//                    }
//                }
//            }
//        }
//        return min;
//    }

    int mod(const int &a, const int &b) {
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
        for (int j = 0; j < M.size() - 1; j++) {
            for (int i = 0; i < E.size() - 1; i++) {
                if (dos(i, j) == 0) { continue; }
                vol += (E(i + 1) - E(i)) * (M(j + 1) - M(j)) * dos(i, j);
            }
        }
        return vol;
    }


    int volume_idx(const MatrixXd &dos, const MatrixXd &E, const MatrixXd &M, const double &vol_limit) {
        double vol = 0;
        for (int i = 0; i < E.size() - 1; i++) {
            for (int j = 0; j < M.size() - 1; j++) {
                if (dos(i, j) == 0) { continue; }
                vol += (E(i + 1) - E(i)) * (M(j + 1) - M(j)) * dos(i, j);
            }
            if (vol > vol_limit){
                return i;
            }
        }
        return (int)E.size() - 1;
    }

    int     find_matching_slope(const MatrixXd &dos1, const MatrixXd &dos2,
                                const VectorXd &E1  , const VectorXd &E2,
                                const VectorXd &M1  , const VectorXd &M2){
        Vector3d u1, u2;                //Vectors connecting adjacent 3 orthogonal points on DOS
        Vector3d v_up, v_rt, v_dn, v_lf; //Vectors connecting adjacent 3 orthogonal points on DOS
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
            for (int j = 1; j < M1.size() - 1; j++) {
                if (E1(i) <= E2.minCoeff()) { continue; }
                if (E1(i) >= E2.maxCoeff()) { continue; }
                if (M1(j) >= M2.maxCoeff()) { continue; }
                if (M1(j) <= M2.minCoeff()) { continue; }
                x = math::upper_bound(E2.data(), E1(i), E2.size());
                y = math::upper_bound(M2.data(), M1(i), M2.size());
                if (x <= 0 || x >= E2.size() - 1) { continue; }
                if (y <= 0 || y >= M2.size() - 1) { continue; }

                v_up << E1(i - 1) - E1(i), 0, dos1(i - 1, j) - dos1(i, j);
                v_dn << E1(i + 1) - E1(i), 0, dos1(i + 1, j) - dos1(i, j);
                v_rt << 0, M1(j + 1) - M1(j), dos1(i, j + 1) - dos1(i, j);
                v_lf << 0, M1(j - 1) - M1(j), dos1(i, j - 1) - dos1(i, j);
                u1 = (v_up.cross(v_rt) + v_rt.cross(v_dn) + v_dn.cross(v_lf) + v_lf.cross(v_dn)).normalized();

                v_up << E2(x - 1) - E2(x), 0, dos2(x - 1, y) - dos2(x, y);
                v_dn << E2(x + 1) - E2(x), 0, dos2(x + 1, y) - dos2(x, y);
                v_rt << 0, M2(y + 1) - M2(y), dos2(x, y + 1) - dos2(x, y);
                v_lf << 0, M2(y - 1) - M2(y), dos2(x, y - 1) - dos2(x, y);

                u2 = (v_up.cross(v_rt) + v_rt.cross(v_dn) + v_dn.cross(v_lf) + v_lf.cross(v_dn)).normalized();
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
