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
//    extern int find_max_idx(const Ref<const ArrayXd>  &list);

    //Find the minimum element larger than zero in a array
//    extern int      find_min_positive(const arrayXi &);
//    extern double   find_min_positive(const ArrayXd &);
    extern int      mod2(const int &,const int &);
    inline int      mod (const int &x, const int &y){
        return x >= 0 ? x%y : x%y + y;
    }
    extern double   volume(const ArrayXXd &dos,const ArrayXd &E,const ArrayXd &M);
    int             volume_idx(const ArrayXXd &dos,const ArrayXd &E,const ArrayXd &M, const double &vol_limit);
    extern Vector3d gradient_vector(const ArrayXXd &dos, const ArrayXd &E, const ArrayXd &M, const int &i, const int &j);
    extern int     find_matching_slope(const ArrayXXd &dos1, const ArrayXXd &dos2,
                                       const ArrayXd &E1  , const ArrayXd &E2,
                                       const ArrayXd &M1  , const ArrayXd &M2);

    template <typename Derived>
    double typical_spacing(const ArrayBase<Derived> &array){
        //Compute the median of the differences to get a typical spacing
        if (array.size() < 2){
            return 1;
        }
        ArrayXd arr  = array;
        ArrayXd diff = arr.tail(arr.size() - 1) - arr.head(arr.size()-1);
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
    double find_min_positive(const ArrayBase<Derived> &array){
        auto min = 1000000000;
        for (int j = 0; j < array.cols(); j++) {
            for (int i = 0; i < array.rows(); i++) {
                if (isnan(array(i,j))){continue;}
                if (array(i, j) < min) {
                    if (array(i, j) > 0) {
                        min = array(i, j);
                    }
                }
            }
        }
        return min;
    }

    template <typename Derived>
    void subtract_min_nonzero(ArrayBase<Derived> &array){
        auto min_positive = find_min_positive(array);
        for (int j = 0; j < array.cols();j++){
            for(int i = 0; i < array.rows(); i++){
                if (array(i,j) == 0){continue;}
                array(i,j) = fmax(0, array(i,j) - min_positive);
            }
        }
    }

    template <typename Derived>
    void subtract_min_nonzero_nan(ArrayBase<Derived> &array){
        array = (array <= 0).select(std::numeric_limits<double>::quiet_NaN(), array);
        auto min_positive = find_min_positive(array);
        for (int j = 0; j < array.cols();j++){
            for(int i = 0; i < array.rows(); i++){
                if (isnan(array(i,j))){continue;}
                array(i,j) = fmax(0, array(i,j) - min_positive);
            }
        }
    }

    template <typename Derived, typename T>
    void add_to_nonzero(ArrayBase<Derived> &array, const T &x){
        for (int j = 0; j < array.cols();j++){
            for(int i = 0; i < array.rows(); i++){
                if (array(i,j) == 0){continue;}
                array(i,j) = fmax(0, array(i,j) + x);
            }
        }
    }


    template <typename Derived1, typename Derived2>
    void remove_nan_rows(ArrayBase<Derived1> &mat, ArrayBase<Derived2> &vec){
        ArrayXXd  mat_temp;
        ArrayXd   vec_temp;
        mat_temp = mat;
        vec_temp = vec;

        int k = 0;
        for (int j = 0; j < mat.rows(); j++){
            if ((mat.row(j) == mat.row(j)).any()){
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
        const auto& no_nan_this = (array == array).select(array,0);
        return no_nan_this.dot(no_nan_this);
    }

    template <typename Derived>
    inline double nansum(const ArrayBase<Derived> & array)  {
        const auto& no_nan_this = (array == array).select(array,0);
        return no_nan_this.sum();
    }
    template <typename Derived>
    inline ArrayXd nansum_rowwise(const ArrayBase<Derived> & array)  {
        const auto& no_nan_this = (array == array).select(array,0);
        return no_nan_this.rowwise().sum();
    }

    template <typename Derived>
    inline ArrayXd nansum_colwise(const ArrayBase<Derived> & array)  {
        const auto& no_nan_this = (array == array).select(array,0);
        return no_nan_this.colwise().sum();
    }
    template <typename Derived>
    inline double nanmaxCoeff(const ArrayBase<Derived> & array)  {
        const auto& no_nan_this = (array == array).select(array,0);
        return no_nan_this.maxCoeff();
    }

    template <typename Derived>
    inline ArrayXd nanmaxCoeff_rowwise(const ArrayBase<Derived> & array)  {
        const auto& no_nan_this = (array == array).select(array,0);
        return no_nan_this.rowwise().maxCoeff();
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
