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

    extern ArrayXXd mean_depthwise(std::vector<ArrayXXd> &array3D) ;
    extern ArrayXXd err_depthwise(std::vector<ArrayXXd> &array3D, ArrayXXd &avg) ;


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
    typename Derived::Scalar nansquared(const ArrayBase<Derived> & array)  {
        return (array == array).select(array,0).cwiseAbs2();
    }

    template <typename Derived>
    typename Derived::Scalar nansum(const ArrayBase<Derived> & array)  {
        return (array == array).select(array,0).sum();
    }
    template <typename Derived>
    typename Derived::Scalar count_num_elements(const ArrayBase<Derived> & array)  {
        return (array == array && array>0).select(ArrayXd::Ones(array.rows(),array.cols()),0).sum();
    }
    template <typename Derived>
    typename Derived::PlainObject nansum_rowwise(const ArrayBase<Derived> & array)  {
        return (array == array).select(array,0).rowwise().sum();
    }

    template <typename Derived>
    Array<typename Derived::Scalar,Dynamic,Dynamic> nansum_colwise(const ArrayBase<Derived> & array)  {
        return (array == array).select(array,0).colwise().sum();
    }

    template <typename Derived>
    typename Derived::Scalar nanmaxCoeff(const ArrayBase<Derived> & array)  {
        return (array == array).select(array,0).maxCoeff();
    }
    template <typename Derived>
    int nanmaxCoeff_idx(const ArrayBase<Derived> & array)  {
        int idx;
        const auto& temp =  (array == array).select(array,0).maxCoeff(&idx);
        return idx;
    }

    template <typename Derived>
    typename Derived::Scalar nanminCoeff(const ArrayBase<Derived> & array)  {
        return  (array == array).select(array,nanmaxCoeff(array)).minCoeff();
    }

    template <typename Derived>
    typename Derived::PlainObject nanmaxCoeff_rowwise(const ArrayBase<Derived> & array)  {
        return (array == array).select(array,0).rowwise().maxCoeff();
    }

    template <typename Derived>
    typename Derived::Scalar nanmean(const ArrayBase<Derived> & array)  {
        return (array == array).select(array,0).mean();
    }

    template <typename Derived>
    typename Derived::PlainObject Zero_to_NaN(const ArrayBase<Derived> &array){
        return (array > 0).select(array, std::numeric_limits<typename Derived::Scalar>::quiet_NaN());
    }

    template <typename Derived>
    typename Derived::PlainObject NaN_to_Zero(const ArrayBase<Derived> &array){
        return (array == array).select(array,0);

    }

    template <typename Derived>
    typename Derived::Scalar find_min_positive(const ArrayBase<Derived> &array){
        return nanminCoeff((array > 0 && array == array ).select(array, nanmaxCoeff(array)));
    }

    template <typename Derived>
    void subtract_min_nonnan(ArrayBase<Derived> &array){
        auto min_nonnan = nanminCoeff(array);
        array = (array > 0 && array == array).select(array-min_nonnan, array);
    }

    template <typename Derived>
    void subtract_min_nonzero(ArrayBase<Derived> &array){
        auto min_positive = find_min_positive(array);
        array = (array > 0 && array == array).select(array-min_positive, array);
    }


    template <typename Derived>
    void subtract_min_nonzero_nan(ArrayBase<Derived> &array){
        array = Zero_to_NaN(array);
        auto min_positive = find_min_positive(array);
        array = (array == array).select(array-min_positive, array);
    }

    template <typename Derived>
    void subtract_min_nonzero_one(ArrayBase<Derived> &array){
        auto min_positive = find_min_positive(array);
        array = (array > 0 && array == array).select(array-min_positive+1, array);
    }



    template <typename Derived, typename T>
    void add_to_nonzero_nonnan(ArrayBase<Derived> &array, const T &x){
        array = (array > 0 && array == array).select(array+x, array);
    }


    template <typename Derived1, typename Derived2>
    void remove_nan_rows(ArrayBase<Derived1> &mat, ArrayBase<Derived2> &vec){
        Array<typename Derived1::Scalar,Dynamic, Dynamic>  mat_temp;
        Array<typename Derived2::Scalar,Dynamic, Dynamic>  vec_temp;
        mat_temp = mat;
        vec_temp = vec;

        int k = 0;
        for (int j = 0; j < mat.rows(); j++){
            if ((mat.row(j) == mat.row(j)).any()){
                mat_temp.row(k) = mat.row(j);
                vec_temp.row(k) = vec.row(j);
                k++;
            }
        }
        mat = mat_temp.topRows(k);
        vec = vec_temp.topRows(k);
    }

//
//
//    template <typename Derived, typename T>
//    inline int binary_search(const ArrayBase<Derived> &list , const T& x){
//        //Now find the point in list closest to x
//        auto low  = std::lower_bound(list.derived().data(), list.derived().data() + list.size(), x);
//        if (low-list.derived().data() >= list.size() ){
//            low--;
//        }
//        return  low-list.derived().data();
//
//    }
//
//    template <typename Derived, typename T, typename T_idx>
//    inline int binary_search(const ArrayBase<Derived> &list , const T& x, const T &y, const T_idx &y_idx){
//        //Now find the point in list closest to x, from below
//        if (x == y){
//            return y_idx;
//        }
//        if (x > y){
//            auto low = std::lower_bound(list.derived().data() + y_idx, list.derived().data() + list.size(), x);
//            while (low-list.derived().data() >= list.size() ){
//                low--;
//            }
//            return  low-list.derived().data();
//
//        }
//        else if(x < y){
//            auto low =  std::lower_bound(list.derived().data(), list.derived().data() + y_idx, x) ;
//            while (low-list.derived().data() >= list.size() ){
//                low--;
//            }
//            return  low-list.derived().data();
//
//        }
//
//    }

//Finds the element nearest x in an Eigen array
    template <typename Derived, typename T>
    inline int binary_search(const ArrayBase<Derived> &list , const T& x){
        //Now find the point in list closest to x

        //CPP REFERENCE lower_bound: Iterator pointing to the first element that is not less than value,
        // or last if no such element is found.
        auto idx  = std::lower_bound(list.derived().data(), list.derived().data() + list.size(), x) - list.derived().data() ;
        //This number idx is potentially out of bounds, one past last element!!
        idx = idx >= list.size() ? idx-1: idx;

        if (list(idx) == x){
            return idx;
        }
        if(idx > 0) {
            if (fabs(list(idx - 1) - x) < fabs(list(idx) - x)) { idx--; }
        }

        return idx;

    }

//Finds the element nearest x in an Eigen array if we already know the index of the current value
    template <typename Derived, typename T, typename T_idx>
    inline int binary_search(const ArrayBase<Derived> &list , const T& x, const T &y, const T_idx &y_idx) {
        //Now find the point in list closest to x, from below

        //CPP REFERENCE lower_bound: Iterator pointing to the first element that is not less than value,
        // or last if no such element is found.
//        std::cout << "x = "<< x << " size = " << list.size() << " y = " << y <<" y_idx = " << y_idx << std::endl;
        if (x == y) {
            return  y_idx >= list.size() ? y_idx -1 : y_idx;
        }

        if (x > y) {
            auto idx = std::lower_bound(list.derived().data() + y_idx, list.derived().data() + list.size(), x) -
                       list.derived().data();
            idx = idx >= list.size() ? idx-1: idx;
            if (list(idx) == x) {
                return idx;
            }
            if (idx > 0) {
                if (fabs(list(idx - 1) - x) < fabs(list(idx) - x)) { idx--; }
            }
            return idx;

        } else {
            auto idx =  std::lower_bound(list.derived().data(), list.derived().data() + y_idx, x) - list.derived().data();
            idx = idx >= list.size() ? idx-1: idx;
            if (list(idx) == x) {
                return idx;
            }
            if (idx > 0) {
                if (fabs(list(idx - 1) - x) < fabs(list(idx) - x)) { idx--; }
            }
            return idx;
        }
    }
}

#endif //WL_NMSPC_MATH_ALGORITHMS_H
