//
// Created by david on 2016-08-05.
//

#ifndef NMSPC_MATH_ALGORITHMS_H
#define NMSPC_MATH_ALGORITHMS_H
#include <Eigen/Core>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace math {
    // Find index of maximum element in an Eigen-type array
    inline int __attribute__((always_inline)) mod(const int x, const int y) {
        //        return x >= 0 ? x%y : x%y + y;
        return (x % y + y) % y;
    }
    extern double volume(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M);
    int           volume_idx(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M, const double &vol_limit);
    extern double area(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M);
    int           area_idx(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M, const double &area_limit);

    extern Eigen::Vector3d gradient_vector(const Eigen::ArrayXXd &dos, const Eigen::ArrayXd &E, const Eigen::ArrayXd &M, const long &i, const long &j);
    extern int             find_matching_slope(const Eigen::ArrayXXd &dos1, const Eigen::ArrayXXd &dos2, const Eigen::ArrayXd &E1, const Eigen::ArrayXd &E2,
                                               const Eigen::ArrayXd &M1, const Eigen::ArrayXd &M2);

    extern Eigen::ArrayXXd mean_depthwise(std::vector<Eigen::ArrayXXd> &array3D);
    extern Eigen::ArrayXXd err_depthwise(std::vector<Eigen::ArrayXXd> &array3D, Eigen::ArrayXXd &avg);

    template<typename Derived>
    double typical_spacing(const Eigen::ArrayBase<Derived> &array) {
        // Compute the median of the differences to get a typical spacing
        if(array.size() < 2) { return 1; }
        Eigen::ArrayXd arr  = array;
        Eigen::ArrayXd diff = arr.tail(arr.size() - 1) - arr.head(arr.size() - 1);
        std::sort(diff.data(), diff.data() + diff.size());
        int last = (int) diff.size() - 1;
        if(mod(last, 2)) {
            // diff has even number of elements
            return (diff(last / 2) + diff(last / 2 + 1)) / 2.0;

        } else {
            // diff has odd number of elements
            return diff(last / 2);
        }
    }

    template<typename Derived>
    typename Derived::Scalar nansquared(const Eigen::ArrayBase<Derived> &array) {
        return (array == array).select(array, 0).cwiseAbs2();
    }

    template<typename Derived>
    typename Derived::Scalar nansum(const Eigen::ArrayBase<Derived> &array) {
        return (array == array).select(array, 0).sum();
    }

    template<typename Derived>
    typename Derived::PlainObject nansum_rowwise(const Eigen::ArrayBase<Derived> &array) {
        return (array == array).select(array, 0).rowwise().sum();
    }

    template<typename Derived>
    Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> nansum_colwise(const Eigen::ArrayBase<Derived> &array) {
        return (array == array).select(array, 0).colwise().sum();
    }

    template<typename Derived>
    inline typename Derived::Scalar nanmaxCoeff(const Eigen::ArrayBase<Derived> &array) {
        return (array == array).select(array, 0).maxCoeff();
    }
    template<typename Derived>
    int nanmaxCoeff_idx(const Eigen::ArrayBase<Derived> &array) {
        int         idx;
        const auto &temp = (array == array).select(array, 0).maxCoeff(&idx);
        (void) temp;
        return idx;
    }

    template<typename Derived>
    typename Derived::Scalar nanminCoeff(const Eigen::ArrayBase<Derived> &array) {
        return (array == array).select(array, nanmaxCoeff(array)).minCoeff();
    }

    template<typename Derived>
    typename Derived::PlainObject nanmaxCoeff_rowwise(const Eigen::ArrayBase<Derived> &array) {
        return (array == array).select(array, 0).rowwise().maxCoeff();
    }

    template<typename Derived>
    typename Derived::Scalar nanmean(const Eigen::ArrayBase<Derived> &array) {
        return (array == array).select(array, 0).mean();
    }

    template<typename Derived>
    typename Derived::Scalar nanzeromean(const Eigen::ArrayBase<Derived> &array) {
        double sum   = 0;
        int    count = 0;
        for(int j = 0; j < array.cols(); j++) {
            for(int i = 0; i < array.rows(); i++) {
                if(array(i, j) == 0) { continue; }
                if(std::isnan(array(i, j))) { continue; }
                sum += array(i, j);
                count++;
            }
        }
        return sum / std::max(1, count);
    }

    template<typename Derived1, typename Derived2>
    typename Derived1::Scalar nanzerostd(const Eigen::ArrayBase<Derived1> &distance, const Eigen::ArrayBase<Derived2> &array1,
                                         const Eigen::ArrayBase<Derived2> &array2) {
        double sum   = 0;
        double m     = nanzeromean(distance); // mean
        int    count = 0;
        for(int j = 0; j < distance.cols(); j++) {
            for(int i = 0; i < distance.rows(); i++) {
                if(array1(i, j) == 0 || array2(i, j) == 0) { continue; }
                if(std::isnan(distance(i, j))) { continue; }
                sum += (distance(i, j) - m) * (distance(i, j) - m);
                count++;
            }
        }
        return sqrt(sum / std::max(1, count));
    }

    template<typename Derived>
    typename Derived::Scalar dos_distance(const Eigen::ArrayBase<Derived> &array1, const Eigen::ArrayBase<Derived> &array2) {
        if(array2.size() != array1.size()) {
            //            std::cout << "Warning: Array size mismatch in dos_distance" << std::endl;
            return nanzeromean(array1) - nanzeromean(array2);
        }
        Eigen::ArrayXXd distance = array1 - array2;
        auto            stDev    = nanzerostd(distance, array1, array2);
        double          m        = nanzeromean(distance); // mean
        double          sum      = 0;
        int             count    = 0;
        for(int j = 0; j < array1.cols(); j++) {
            for(int i = 0; i < array1.rows(); i++) {
                if(array1(i, j) == 0 || array2(i, j) == 0) { continue; }
                if(std::isnan(distance(i, j))) { continue; }
                if(std::fabs(distance(i, j) - m) < 0.5 * std::fabs(stDev)) {
                    sum += distance(i, j);
                    count++;
                }
            }
        }
        if(count == 0) {
            return nanzeromean(array1) - nanzeromean(array2);
        } else {
            return sum / std::max(1, count);
        }
    }

    template<typename Derived>
    typename Derived::PlainObject Zero_to_NaN(const Eigen::ArrayBase<Derived> &array) {
        return (array > 0).select(array, std::numeric_limits<typename Derived::Scalar>::quiet_NaN());
    }

    template<typename Derived>
    typename Derived::PlainObject NaN_to_Zero(const Eigen::ArrayBase<Derived> &array) {
        return (array == array).select(array, 0);
    }

    template<typename Derived>
    int count_num_elements(const Eigen::ArrayBase<Derived> &array) {
        auto temp = Zero_to_NaN(array);
        return (int) nansum((temp == temp).select(temp / temp, 0));
    }

    //    template <typename Derived>
    //    inline typename Derived::Scalar find_min_positive(const Eigen::ArrayBase<Derived> &array){
    //        return nanminCoeff((array > 0 && array == array ).select(array, nanmaxCoeff(array)));
    //    }

    template<typename Derived>
    inline typename Derived::Scalar find_min_positive(const Eigen::ArrayBase<Derived> &array) {
        return ((array > 0 && array == array).select(array, array.maxCoeff())).minCoeff();
    }

    template<typename Derived>
    void subtract_min_nonnan(Eigen::ArrayBase<Derived> &array) {
        auto min_nonnan = nanminCoeff(array);
        array           = (array > 0 && array == array).select(array - min_nonnan, array);
    }

    template<typename Derived, typename T>
    void add_to_nonzero_nonnan(Eigen::ArrayBase<Derived> &array, const T &x) {
        array = (array > 0 && array == array).select(array + x, array);
        array = (array < 0).select(0, array);
    }

    template<typename Derived>
    void subtract_min_nonzero(Eigen::ArrayBase<Derived> &array) {
        auto min_positive = find_min_positive(array);
        array             = (array > 0 && array == array).select(array - min_positive, array);
    }

    template<typename Derived>
    void subtract_min_nonzero_nan(Eigen::ArrayBase<Derived> &array) {
        array             = Zero_to_NaN(array);
        auto min_positive = find_min_positive(array);
        array             = (array == array).select(array - min_positive, array);
    }

    template<typename Derived>
    inline void subtract_min_nonzero_one(Eigen::ArrayBase<Derived> &array) {
        auto min_positive = find_min_positive(array);
        array             = (array > 0 && array == array).select(array - min_positive + 1, array);
    }

    template<typename Derived1, typename Derived2>
    void remove_nan_rows(Eigen::ArrayBase<Derived1> &mat, Eigen::ArrayBase<Derived2> &vec) {
        Eigen::Array<typename Derived1::Scalar, Eigen::Dynamic, Eigen::Dynamic> mat_temp;
        Eigen::Array<typename Derived2::Scalar, Eigen::Dynamic, Eigen::Dynamic> vec_temp;
        mat_temp = mat;
        vec_temp = vec;

        int k    = 0;
        for(int j = 0; j < mat.rows(); j++) {
            if((mat.row(j) == mat.row(j)).any()) {
                mat_temp.row(k) = mat.row(j);
                vec_temp.row(k) = vec.row(j);
                k++;
            }
        }
        mat = mat_temp.topRows(k);
        vec = vec_temp.topRows(k);
    }

    // Finds the element nearest x in an Eigen array
    template<typename Derived, typename T>
    inline long binary_search_nearest(const Eigen::ArrayBase<Derived> &list, const T x) {
        // Now find the point in list closest to x

        // CPP REFERENCE lower_bound: Iterator pointing to the first element that is not less than value,
        //  or last if no such element is found.
        auto idx = std::lower_bound(list.derived().data(), list.derived().data() + list.size(), x) - list.derived().data();
        // This number idx is potentially out of bounds, one past last element!!
        idx = idx >= list.size() ? idx - 1 : idx;

        if(list(idx) == x) { return idx; }
        if(idx > 0) {
            if(fabs(list(idx - 1) - x) < fabs(list(idx) - x)) { idx--; }
        }

        return idx;
    }

    // Finds the element nearest x in an Eigen array if we already know the index of the current value
    template<typename Derived, typename T, typename T_idx>
    inline int binary_search_nearest(const Eigen::ArrayBase<Derived> &list, const T x, const T y, const T_idx y_idx) {
        // Now find the point in list closest to x, from below

        // CPP REFERENCE lower_bound: Iterator pointing to the first element that is not less than value,
        //  or last if no such element is found.
        //        std::cout << "x = "<< x << " size = " << list.size() << " y = " << y <<" y_idx = " << y_idx << std::endl;
        if(x == y) { return y_idx >= list.size() ? y_idx - 1 : y_idx; }

        if(x > y) {
            auto idx = std::lower_bound(list.derived().data() + y_idx, list.derived().data() + list.size(), x) - list.derived().data();
            idx      = idx >= list.size() ? idx - 1 : idx;
            if(list(idx) == x) { return idx; }
            if(idx > 0) {
                if(fabs(list(idx - 1) - x) < fabs(list(idx) - x)) { idx--; }
            }
            return idx;

        } else {
            auto idx = std::lower_bound(list.derived().data(), list.derived().data() + y_idx, x) - list.derived().data();
            idx      = idx >= list.size() ? idx - 1 : idx;
            if(list(idx) == x) { return idx; }
            if(idx > 0) {
                if(fabs(list(idx - 1) - x) < fabs(list(idx) - x)) { idx--; }
            }
            return idx;
        }
    }

    // Finds the element exactly x in an Eigen array
    template<typename Derived, typename T>
    long binary_search_exact(const Eigen::ArrayBase<Derived> &list, const T x) {
        // Now find the point in list that exactly matches x.
        // If none is found, return -1;
        // CPP REFERENCE lower_bound: Iterator pointing to the first element that is not less than value,
        //  or last if no such element is found.
        auto idx = std::lower_bound(list.derived().data(), list.derived().data() + list.size(), x) - list.derived().data();
        // This number idx is potentially out of bounds, one past last element!!
        if(idx >= list.size()) { return -1; }
        if(list(idx) == x) {
            return idx;
        } else {
            return -1;
        }
    }

    // Finds the element nearest x in an Eigen array if we already know the index of the current value
    template<typename Derived, typename T, typename T_idx>
    int binary_search_exact(const Eigen::ArrayBase<Derived> &list, const T x, const T y, const T_idx y_idx) {
        // Now find the point in list closest to x, from below

        // CPP REFERENCE lower_bound: Iterator pointing to the first element that is not less than value,
        //  or last if no such element is found.
        //         std::cout << "x = "<< x << " size = " << list.size() << " y = " << y <<" y_idx = " << y_idx << std::endl;
        if(y_idx >= list.size()) { return binary_search_exact(list, x); }
        if(list(y_idx) == x) { return y_idx; }
        if(x > y) {
            auto idx = std::lower_bound(list.derived().data() + y_idx, list.derived().data() + list.size(), x) - list.derived().data();
            if(idx >= list.size()) { return -1; }
            if(list(idx) == x) { return idx; }

        } else {
            auto idx = std::lower_bound(list.derived().data(), list.derived().data() + y_idx, x) - list.derived().data();
            if(idx >= list.size()) { return -1; }

            if(list(idx) == x) { return idx; }
        }
        //        std::cout << "Looking for " << x << " Previous y = " << y << " at " << y_idx << " List: " << list.transpose() << std::endl;

        return binary_search_exact(list, x);
    }
}

#endif // NMSPC_MATH_ALGORITHMS_H
