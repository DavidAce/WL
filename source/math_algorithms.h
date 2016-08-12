//
// Created by david on 2016-08-05.
//

#ifndef WL_MATH_ALGORITHMS_H
#define WL_MATH_ALGORITHMS_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

using namespace Eigen;
namespace math{
    //Find index of maximum element in an Eigen-type array
    //extern int find_max_idx(const Ref<const ArrayXd> &list);
//    extern int find_max_idx(const Ref<const VectorXd>  &list);

    //Find the minimum element larger than zero in a matrix
//    extern int      find_min_positive(const MatrixXi &);
//    extern double   find_min_positive(const MatrixXd &);
    extern int      mod(const int &,const int &);
    extern double   volume(const MatrixXd &,const MatrixXd &,const MatrixXd &);
    int             volume_idx(const MatrixXd &,const MatrixXd &,const MatrixXd &, const double &);

    extern int     find_matching_slope(const MatrixXd &dos1, const MatrixXd &dos2,
                                       const VectorXd &E1  , const VectorXd &E2,
                                       const VectorXd &M1  , const VectorXd &M2);
    template <typename Derived>
    int find_min_positive(const MatrixBase<Derived> &matrix){
        int min = 1000000000;
        for (int j = 0; j < matrix.cols(); j++) {
            for (int i = 0; i < matrix.rows(); i++) {
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
        double min_positive = find_min_positive(matrix);
        for (int j = 0; j < matrix.cols();j++){
            for(int i = 0; i < matrix.rows(); i++){
                if (matrix(i,j) == 0){continue;}
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
    //Finds the element nearest x in a C-style array
    template <typename List_type, typename T, typename size_type>
    int binary_search(const List_type &list , const T& x, const size_type &size){
        //Now find the point in list closest to x, from below
        if (size <= 1){return 0;}
        auto low  = std::upper_bound (list, list + size, x);
        int index = low-list-1;
        if (index + 1 < size && fabs(list[index]- x) > fabs(list[index+1] < x)){
            index++;
        }
        return  index;
    }

    //Finds the element nearest x FROM ABOVE in a C-style array
    template <typename List_type, typename T, typename size_type>
    int upper_bound(const List_type &list , const T& x, const size_type &size){
        //Now find the point in list closest to x, from below.
        if (size <= 1){return 0;}
        auto low  = std::upper_bound (list, list + size, x);
        return  low-list-1;
    }

    template <typename List_type, typename T, typename size_type>
    int lower_bound(const List_type &list , const T& x, const size_type &size){
        //Now find the point in list closest to x, from below.
        if (size <= 1){return 0;}
        auto low  = std::lower_bound(list, list + size, x);
        return  low-list-1;
    }



}

#endif //WL_MATH_ALGORITHMS_H
