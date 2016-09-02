//
// Created by david on 2016-08-05.
//

#ifndef WL_NMSPC_MATH_ALGORITHMS_H
#define WL_NMSPC_MATH_ALGORITHMS_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace Eigen;
namespace math{
    //Find index of maximum element in an Eigen-type array
    //extern int find_max_idx(const Ref<const ArrayXd> &list);
//    extern int find_max_idx(const Ref<const VectorXd>  &list);

    //Find the minimum element larger than zero in a matrix
//    extern int      find_min_positive(const MatrixXi &);
//    extern double   find_min_positive(const MatrixXd &);
    extern int      mod2(const int &,const int &);
    inline int      mod (const int &x, const int &y){
        return x >= 0 ? x%y : x%y + y;
    }
    extern double   volume(const MatrixXd &,const MatrixXd &,const MatrixXd &);
    int             volume_idx(const MatrixXd &,const MatrixXd &,const MatrixXd &, const double &);
    extern Vector3d gradient_vector(const MatrixXd &dos, const VectorXd &E, const VectorXd &M, const int &i,
                                    const int &j);
    extern int     find_matching_slope(const MatrixXd &dos1, const MatrixXd &dos2,
                                       const VectorXd &E1  , const VectorXd &E2,
                                       const VectorXd &M1  , const VectorXd &M2);

    template <typename Derived>
    double typical_spacing(const MatrixBase<Derived> &matrix){
        //Compute the median of the differences to get a typical spacing
        if (matrix.size() < 2){
            return 1;
        }
        ArrayXd arr  = matrix.array();
        VectorXd diff = arr.tail(arr.size() - 1) - arr.head(arr.size()-1);
        std::sort(diff.data(), diff.data()+diff.size());
        int last = (int) diff.size() -1;
        if (mod(last, 2)){
            //diff has even number of elements
            return (diff(last/2) + diff(last/2 + 1)) / 2.0;

        }else{
            //diff has odd number of elements
            return diff(last/2);

        }


    }

    template <typename Derived>
    double find_min_positive(const MatrixBase<Derived> &matrix){
        auto min = 1000000000;
        for (int j = 0; j < matrix.cols(); j++) {
            for (int i = 0; i < matrix.rows(); i++) {
                if (isnan(matrix(i,j))){continue;}
                if (matrix(i, j) < min) {
                    if (matrix(i, j) > 0) {
                        min = matrix(i, j);
                    }
                }
            }
        }
        return min;
    }

    template <typename Derived>
    void subtract_min_nonzero(MatrixBase<Derived> &matrix){
        auto min_positive = find_min_positive(matrix);
        for (int j = 0; j < matrix.cols();j++){
            for(int i = 0; i < matrix.rows(); i++){
                if (matrix(i,j) == 0){continue;}
                matrix(i,j) = fmax(0, matrix(i,j) - min_positive);
            }
        }
    }

    template <typename Derived>
    void subtract_min_nonzero_nan(MatrixBase<Derived> &matrix){
        matrix = (matrix.array() <= 0).select(std::numeric_limits<double>::quiet_NaN(), matrix);
        auto min_positive = find_min_positive(matrix);
        std::cout << std::endl << "min Positive " << min_positive << std::endl<<std::endl;
        for (int j = 0; j < matrix.cols();j++){
            for(int i = 0; i < matrix.rows(); i++){
                if (isnan(matrix(i,j))){continue;}
                matrix(i,j) = fmax(0, matrix(i,j) - min_positive);
            }
        }
    }

    template <typename Derived, typename T>
    void add_to_nonzero(MatrixBase<Derived> &matrix, const T &x){
        for (int j = 0; j < matrix.cols();j++){
            for(int i = 0; i < matrix.rows(); i++){
                if (matrix(i,j) == 0){continue;}
                matrix(i,j) = fmax(0, matrix(i,j) + x);
            }
        }
    }


    template <typename Derived1, typename Derived2>
    void remove_nan_rows(MatrixBase<Derived1> &mat, MatrixBase<Derived2> &vec){
        MatrixXd  mat_temp;
        VectorXd  vec_temp;
        mat_temp = mat;
        vec_temp = vec;

        int k = 0;
        for (int j = 0; j < mat.rows(); j++){
            if ((mat.row(j).array() == mat.row(j).array()).any()){
                mat_temp.row(k) = mat.row(j);
                vec_temp(k)     = vec(j);
//                std::cout << k << " " << mat_temp.row(k) << " " << vec_temp(k) << std::endl;
                k++;
            }
        }
        mat = mat_temp.topRows(k);
        vec = vec_temp.head(k);
    }

    template <typename Derived>
    inline double nansquared(const ArrayBase<Derived> & array)  {
        const auto& no_nan_this = (array == array).select(0, array);
        return no_nan_this.dot(no_nan_this);
    }

    template <typename Derived>
    inline double nansum(const ArrayBase<Derived> & array)  {
        const auto& no_nan_this = (array == array).select(0, array);
        return no_nan_this.sum();
    }

    template <typename Derived>
    inline double nansum(const MatrixBase<Derived> & matrix)  {
        const auto& no_nan_this = (matrix.array() == matrix.array()).select(0, matrix);
        return no_nan_this.sum();
    }

    template <typename Derived>
    inline double nanmaxCoeff(const ArrayBase<Derived> & array)  {
        const auto& no_nan_this = (array == array).select(0, array);
        return no_nan_this.maxCoeff();
    }

    template <typename Derived>
    inline double nanmaxCoeff(const MatrixBase<Derived> & matrix)  {
        const auto& no_nan_this = (matrix.array() == matrix.array()).select(0, matrix);
        return no_nan_this.maxCoeff();
    }


    //Finds the element nearest x in a C-style array
    template <typename List_type, typename T, typename size_type>
    inline int binary_search(const List_type &list , const T& x, const size_type &size){
        //Now find the point in list closest to x, from below
        auto low  = std::lower_bound(list, list + size, x);
        if (low-list >= size ){
            low--;
        }
       return  low-list;

    }

    template<typename T> struct TD;
    template <typename List_type, typename T, typename T_idx, typename size_type>
    inline int binary_search(const List_type &list , const T& x, const size_type &size, const T &y, const T_idx &y_idx){
        //Now find the point in list closest to x, from below
        if (x == y){
            return y_idx;
        }
        double *low;
        if (x > y){
            low  = std::lower_bound(list + y_idx, list + size, x);
        }else if(x < y){
            low  = std::lower_bound(list, list + y_idx, x);
        }
        while (low-list >= size ){
            low--;
        }
        return  low-list;

    }


}

#endif //WL_NMSPC_MATH_ALGORITHMS_H
