//
// Created by david on 2016-08-05.
//

#include "nmspc_math_algorithms.h"
#include "nmspc_WL_constants.h"


namespace math {

    ArrayXXd mean_depthwise(std::vector<ArrayXXd> &array3D) {
        auto rows = (int)array3D[0].rows();
        auto cols = (int)array3D[0].cols();
        auto dept = (int)array3D.size();
        ArrayXXd result(rows,cols);
        ArrayXd  depth_array(dept);
        for(auto j = 0 ; j < cols; j++ ){
            for(auto i = 0; i < rows; i++){
                for(auto k = 0; k < dept; k++){
                    depth_array(k) = array3D[k](i,j);
                }
                if(depth_array.hasNaN()){
                    result(i,j) = std::numeric_limits<double>::quiet_NaN();
                }
                result(i,j) = depth_array.mean();
            }
        }
        return result;
    }

    ArrayXXd err_depthwise(std::vector<ArrayXXd> &array3D, ArrayXXd &avg) {
        auto rows = array3D[0].rows();
        auto cols = array3D[0].cols();
        auto dept = array3D.size();
        ArrayXXd result(rows,cols);
        ArrayXd  depth_array(dept);
        for(auto j = 0 ; j < cols; j++ ){
            for(auto i = 0; i < rows; i++){
                for(unsigned long int k = 0; k < dept; k++){
                    depth_array(k) = array3D[k](i,j);
                }
                if (depth_array.hasNaN()){
                    result(i,j) = std::numeric_limits<double>::quiet_NaN();
                }else{
                    result(i,j) = sqrt((depth_array - avg(i,j)).cwiseAbs2().sum()/(dept-1));
                }
//                result(i,j) = depth_array.mean();
            }
        }
        return result;
    }

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

    double volume(const ArrayXXd &dos, const ArrayXd &E, const ArrayXd &M) {
        double vol = 0;
//        ArrayXd dE;
//        ArrayXd dM;
        switch (constants::rw_dims){
            case 1:
//                dE = E.head(E.size()-1) - E.tail(E.size()-1);
//                vol = nansum(dos.topLeftCorner(dE.size(),1) * dE);
                for (int j = 0; j < M.size(); j++) {
                    for (int i = 0; i < E.size() - 1; i++) {
                        if (dos(i, j) == 0 || std::isnan(dos(i,j))) { continue; }
                        vol += (E(i + 1) - E(i)) * dos(i, j);
                    }
                }
                break;
            case 2:
//                dE = E.head(E.size()-1) - E.tail(E.size()-1);
//                dM = M.head(M.size()-1) - M.tail(M.size()-1);
//                vol = nansum(dos.topLeftCorner(dE.size(), dM.size()) * (dE.matrix() * dM.matrix().transpose()).array());
                for (int j = 0; j < M.size() - 1; j++) {
                    for (int i = 0; i < E.size() - 1; i++) {
                        if (dos(i, j) == 0 || std::isnan(dos(i,j))) { continue; }
//                        vol += (E(i + 1) - E(i)) * (M(j + 1) - M(j)) * dos(i, j);
                        vol += (E(i + 1) - E(i)) * (M(j + 1) - M(j)) ;
                    }
                }
                break;
            default:
                std::cout << "Wrong rw-dimension constants::rw_dims " << std::endl;
                exit(5);

        }
        return vol;
    }


    int volume_idx(const ArrayXXd &dos, const ArrayXd &E, const ArrayXd &M, const double &vol_limit) {
        double vol = 0;
        switch(constants::rw_dims){
            case 1:
                for (int i = 0; i < E.size() - 1; i++) {
                    for (int j = 0; j < M.size(); j++) {
                        if (dos(i, j) == 0 || std::isnan(dos(i,j))) { continue; }
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
                        if (dos(i, j) == 0 || std::isnan(dos(i,j))) { continue; }
//                        vol += (E(i + 1) - E(i)) * (M(j + 1) - M(j)) * dos(i, j);
                        vol += (E(i + 1) - E(i)) * (M(j + 1) - M(j)) ;
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

    bool  on_the_edge_up(const ArrayXXd &dos, const ArrayXd &E, const ArrayXd &M, const int &i, const int &j){
        return i == 0 ||  dos(std::max(i-1,0), j) == 0||  std::isnan(dos(std::max(i-1,0), j));
    }

    bool  on_the_edge_dn(const ArrayXXd &dos, const ArrayXd &E, const ArrayXd &M, const int &i, const int &j){
        return i == E.size()-1 || dos(std::min(i+1, (int)E.size()-1), j) == 0 || std::isnan(dos(std::min(i+1, (int)E.size()-1), j));
    }

    bool  on_the_edge_lf(const ArrayXXd &dos, const ArrayXd &E, const ArrayXd &M, const int &i, const int &j){
        return  j == 0 || dos(i, std::max(j-1,0)) == 0 || std::isnan(dos(i, std::max(j-1,0)));
    }

    bool  on_the_edge_rt(const ArrayXXd &dos, const ArrayXd &E, const ArrayXd &M, const int &i, const int &j){
        return j == M.size()-1 || dos(i, std::min(j+1, (int)M.size()-1)) == 0 ||std::isnan(dos(i, std::min(j+1, (int)M.size()-1)));
    }

    Vector3d gradient_vector(const ArrayXXd &dos, const ArrayXd &E, const ArrayXd &M, const int &i, const int &j){
        Vector3d v_up, v_rt, v_dn, v_lf; //Vectors connecting adjacent 3 orthogonal points on DOS
        switch(constants::rw_dims){
            case 1:
                v_dn << E(i + 1) - E(i), 0, dos(i + 1, j) - dos(i, j);
                return v_dn.normalized();
            case 2:
                //Down is in the +x direction, and right in the +y direction (right hand rule)
                if (on_the_edge_up(dos,E,M,i,j)){v_up << -1, 0, 0;}else{v_up << E(i - 1) - E(i), 0, dos(i - 1, j) - dos(i, j);} //Detect if on edge
                if (on_the_edge_dn(dos,E,M,i,j)){v_dn <<  1, 0, 0;}else{v_dn << E(i + 1) - E(i), 0, dos(i + 1, j) - dos(i, j);} //Detect if on edge
                if (on_the_edge_rt(dos,E,M,i,j)){v_rt <<  0, 1, 0;}else{v_rt << 0, M(j + 1) - M(j), dos(i, j + 1) - dos(i, j);} //Detect if on edge
                if (on_the_edge_lf(dos,E,M,i,j)){v_lf <<  0,-1, 0;}else{v_lf << 0, M(j - 1) - M(j), dos(i, j - 1) - dos(i, j);} //Detect if on edge
                return (v_up.cross(v_rt) + v_rt.cross(v_dn) + v_dn.cross(v_lf) + v_lf.cross(v_dn)).normalized();
//                return v_tot.normalized();

        }
    }

    int  find_matching_slope(const ArrayXXd &dos1, const ArrayXXd &dos2,
                             const ArrayXd &E1  , const ArrayXd &E2,
                             const ArrayXd &M1  , const ArrayXd &M2){
        Vector3d u1, u2;                //Vectors connecting adjacent 3 orthogonal points on DOS
//        ArrayXd sum(E1.size());
        ArrayXXd sum(E1.size(), M1.size());
        sum.fill(0);
        int E_merge_idx;    //Index of merging point

        if (dos1.rows() < 3 || dos2.rows() < 3){
            std::cout << "Rows too few to merge!!" << std::endl;
            exit(35);
        }

        int x, y; //Coordinates closest to i, j on dos of neighbor above.
        for (int i = 0; i < E1.size(); i++) {
            if (E1(i) <  E2.minCoeff()) { continue; }
            if (E1(i) >  E2.maxCoeff()) { continue; }
            x = math::binary_search(E2, E1(i));
            for (int j = 0; j < M1.size(); j++) {
                if (dos1(i,j) == 0 || std::isnan(dos1(i,j)))  { continue; }
                y = math::binary_search(M2, M1(j));
                if (dos2(x,y) == 0 ||std::isnan(dos2(x,y)))  { continue; }

                u1 = gradient_vector(dos1, E1, M1, i, j);
                u2 = gradient_vector(dos2, E2, M2, x, y);
//                if (u1.norm() == 0 || u2.norm() == 0){continue;}
                sum(i,j) = u1.dot(u2);
            }
        }

        sum.rowwise().mean().maxCoeff(&E_merge_idx);
//        std::cout << "sum [" << E_merge_idx << "] = " << std::endl <<  sum.rowwise().mean() << std::endl << std::endl;

//        E_merge_idx = sum.sum() == 0 ? (int)E1.size() - 1 : E_merge_idx;
        return E_merge_idx;
    }


}
