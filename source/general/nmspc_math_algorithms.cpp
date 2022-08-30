//
// Created by david on 2016-08-05.
//

#include "nmspc_math_algorithms.h"
#include "params/nmspc_WL_constants.h"
#include <Eigen/Geometry>

namespace math {

    Eigen::ArrayXXd mean_depthwise(std::vector<Eigen::ArrayXXd> &array3D) {
        auto            rows = (int) array3D[0].rows();
        auto            cols = (int) array3D[0].cols();
        auto            dept = (int) array3D.size();
        Eigen::ArrayXXd result(rows, cols);
        Eigen::ArrayXd  depth_array(dept);
        for(auto j = 0; j < cols; j++) {
            for(auto i = 0; i < rows; i++) {
                for(auto k = 0; k < dept; k++) { depth_array(k) = array3D[k](i, j); }
                if(depth_array.hasNaN()) { result(i, j) = std::numeric_limits<double>::quiet_NaN(); }
                result(i, j) = depth_array.mean();
            }
        }
        return result;
    }

    Eigen::ArrayXXd err_depthwise(std::vector<Eigen::ArrayXXd> &array3D, Eigen::ArrayXXd &avg) {
        auto            rows = array3D[0].rows();
        auto            cols = array3D[0].cols();
        auto            dept = array3D.size();
        Eigen::ArrayXXd result(rows, cols);
        Eigen::ArrayXd  depth_array(dept);
        for(auto j = 0; j < cols; j++) {
            for(auto i = 0; i < rows; i++) {
                for(unsigned long int k = 0; k < dept; k++) { depth_array(static_cast<long>(k)) = array3D[k](i, j); }
                if(depth_array.hasNaN()) {
                    result(i, j) = std::numeric_limits<double>::quiet_NaN();
                } else {
                    result(i, j) = sqrt((depth_array - avg(i, j)).cwiseAbs2().sum() / static_cast<double>(dept - 1));
                }
                //                result(i,j) = depth_array.mean();
            }
        }
        return result;
    }

    int mod2(const int &a, const int &b) {
        if(b < 0) { return mod(-a, -b); }
        int ret = a % b;
        if(ret < 0) { ret += b; }
        return ret;
    }

    double volume(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M) {
        double vol = 0;
        switch(constants::rw_dims) {
            case 1:
                for(int j = 0; j < M.size(); j++) {
                    for(int i = 0; i < E.size() - 1; i++) {
                        if(dos(i, j) == 0 || std::isnan(dos(i, j))) { continue; }
                        vol += (E(i + 1) - E(i)) * dos(i, j);
                    }
                }
                break;
            case 2:
                for(int j = 0; j < M.size() - 1; j++) {
                    for(int i = 0; i < E.size() - 1; i++) {
                        if(dos(i, j) == 0 || std::isnan(dos(i, j))) { continue; }
                        vol += (E(i + 1) - E(i)) * (M(j + 1) - M(j)) * dos(i, j);
                    }
                }
                break;
            default: std::cout << "Wrong rw-dimension constants::rw_dims " << std::endl; exit(5);
        }
        return vol;
    }

    int volume_idx(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M, const double &vol_limit) {
        double vol = 0;
        switch(constants::rw_dims) {
            case 1:
                for(int i = 0; i < E.size() - 1; i++) {
                    for(int j = 0; j < M.size(); j++) {
                        if(dos(i, j) == 0 || std::isnan(dos(i, j))) { continue; }
                        vol += (E(i + 1) - E(i)) * dos(i, j);
                    }
                    if(vol >= vol_limit) { return i; }
                }
                break;
            case 2:
                for(int i = 0; i < E.size() - 1; i++) {
                    for(int j = 0; j < M.size() - 1; j++) {
                        if(dos(i, j) == 0 || std::isnan(dos(i, j))) { continue; }
                        vol += (E(i + 1) - E(i)) * (M(j + 1) - M(j)) * dos(i, j);
                    }
                    if(vol >= vol_limit) { return i; }
                }
                break;
            default: std::cout << "Wrong rw-dimension constants::rw_dims " << std::endl; exit(5);
        }
        return (int) E.size() - 1;
    }

    double area(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M) {
        double area = 0;
        switch(constants::rw_dims) {
            case 1:
                for(int j = 0; j < M.size(); j++) {
                    for(int i = 0; i < E.size() - 1; i++) {
                        if(dos(i, j) == 0 || std::isnan(dos(i, j))) { continue; }
                        area += (E(i + 1) - E(i));
                    }
                }
                break;
            case 2:
                for(int j = 0; j < M.size() - 1; j++) {
                    for(int i = 0; i < E.size() - 1; i++) {
                        if(dos(i, j) == 0 || std::isnan(dos(i, j))) { continue; }
                        area += (E(i + 1) - E(i)) * (M(j + 1) - M(j));
                    }
                }
                break;
            default: std::cout << "Wrong rw-dimension constants::rw_dims " << std::endl; exit(5);
        }
        return area;
    }

    int area_idx(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M, const double &area_limit) {
        double area = 0;
        switch(constants::rw_dims) {
            case 1:
                for(int i = 0; i < E.size() - 1; i++) {
                    for(int j = 0; j < M.size(); j++) {
                        if(dos(i, j) == 0 || std::isnan(dos(i, j))) { continue; }
                        area += (E(i + 1) - E(i));
                    }
                    if(area >= area_limit) { return i; }
                }
                break;
            case 2:
                for(int i = 0; i < E.size() - 1; i++) {
                    for(int j = 0; j < M.size() - 1; j++) {
                        if(dos(i, j) == 0 || std::isnan(dos(i, j))) { continue; }
                        area += (E(i + 1) - E(i)) * (M(j + 1) - M(j));
                    }
                    if(area >= area_limit) { return i; }
                }
                break;
            default: std::cout << "Wrong rw-dimension constants::rw_dims " << std::endl; exit(5);
        }
        return (int) E.size() - 1;
    }

    bool on_the_edge_up(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M, const long &i, const long &j) {
        return i == 0 || dos(std::max(i - 1, 0l), j) == 0 || std::isnan(dos(std::max(i - 1, 0l), j));
    }

    bool on_the_edge_dn(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M, const long &i, const long &j) {
        return i == E.size() - 1 || dos(std::min(i + 1, E.size() - 1), j) == 0 || std::isnan(dos(std::min(i + 1, E.size() - 1), j));
    }

    bool on_the_edge_lf(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M, const long &i, const long &j) {
        return j == 0 || dos(i, std::max(j - 1, 0l)) == 0 || std::isnan(dos(i, std::max(j - 1, 0l)));
    }

    bool on_the_edge_rt(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M, const long &i, const long &j) {
        return j == M.size() - 1 || dos(i, std::min(j + 1, M.size() - 1)) == 0 || std::isnan(dos(i, std::min(j + 1, M.size() - 1)));
    }

    Eigen::Vector3d gradient_vector(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M, const long &i, const long &j) {
        Eigen::Vector3d v_up, v_rt, v_dn, v_lf; // Vectors connecting adjacent 3 orthogonal points on DOS
        switch(constants::rw_dims) {
            case 1: v_dn << E(i + 1) - E(i), 0, dos(i + 1, j) - dos(i, j); return v_dn.normalized();
            case 2:
                // Down is in the +x direction, and right in the +y direction (right hand rule)
                if(on_the_edge_up(dos, E, M, i, j)) {
                    v_up << -1, 0, 0;
                } else {
                    v_up << E(i - 1) - E(i), 0, dos(i - 1, j) - dos(i, j);
                } // Detect if on edge
                if(on_the_edge_dn(dos, E, M, i, j)) {
                    v_dn << 1, 0, 0;
                } else {
                    v_dn << E(i + 1) - E(i), 0, dos(i + 1, j) - dos(i, j);
                } // Detect if on edge
                if(on_the_edge_rt(dos, E, M, i, j)) {
                    v_rt << 0, 1, 0;
                } else {
                    v_rt << 0, M(j + 1) - M(j), dos(i, j + 1) - dos(i, j);
                } // Detect if on edge
                if(on_the_edge_lf(dos, E, M, i, j)) {
                    v_lf << 0, -1, 0;
                } else {
                    v_lf << 0, M(j - 1) - M(j), dos(i, j - 1) - dos(i, j);
                } // Detect if on edge
                return (v_up.cross(v_rt) + v_rt.cross(v_dn) + v_dn.cross(v_lf) + v_lf.cross(v_dn)).normalized();
                //                return v_tot.normalized();
            default: throw std::runtime_error("constants::rw_dims must be 1 or 2");
        }
    }

    int find_matching_slope(const Eigen::ArrayXXd &dos1, const Eigen::ArrayXXd &dos2, const Eigen::ArrayXd &E1, const Eigen::ArrayXd &E2,
                            const Eigen::ArrayXd &M1, const Eigen::ArrayXd &M2) {
        Eigen::Vector3d u1, u2; // Vectors connecting adjacent 3 orthogonal points on DOS
        //        Eigen::ArrayXd sum(E1.size());
        Eigen::ArrayXXd sum(E1.size(), M1.size());
        sum.fill(0);
        int E_merge_idx; // Index of merging point

        if(dos1.rows() < 1 || dos2.rows() < 1) {
            std::cout << "Rows too few to merge!!" << std::endl;
            exit(35);
        }

        long x, y; // Coordinates closest to i, j on dos of neighbor above.
        for(int i = 0; i < E1.size(); i++) {
            if(E1(i) < E2.minCoeff()) { continue; }
            if(E1(i) > E2.maxCoeff()) { continue; }
            x = math::binary_search_nearest(E2, E1(i));
            for(int j = 0; j < M1.size(); j++) {
                if(dos1(i, j) == 0 || std::isnan(dos1(i, j))) { continue; }
                y = math::binary_search_nearest(M2, M1(j));
                if(dos2(x, y) == 0 || std::isnan(dos2(x, y))) { continue; }

                u1 = gradient_vector(dos1, E1, M1, i, j);
                u2 = gradient_vector(dos2, E2, M2, x, y);
                //                if (u1.norm() == 0 || u2.norm() == 0){continue;}
                sum(i, j) = u1.dot(u2);
            }
        }

        sum.rowwise().mean().maxCoeff(&E_merge_idx);
        //        std::cout << "sum [" << E_merge_idx << "] = " << std::endl <<  sum.rowwise().mean() << std::endl << std::endl;

        //        E_merge_idx = sum.sum() == 0 ? (int)E1.size() - 1 : E_merge_idx;
        return E_merge_idx;
    }

}
