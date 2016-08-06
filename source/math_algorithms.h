//
// Created by david on 2016-08-05.
//

#ifndef WL_MATH_ALGORITHMS_H
#define WL_MATH_ALGORITHMS_H
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;
namespace math{
    //Find index of maximum element in an Eigen-type array
    //extern int find_max_idx(const Ref<const ArrayXd> &list);
    extern int find_max_idx(const Ref<const VectorXd>  &list);

    //Find the minimum element larger than zero in a matrix
    extern int find_min_positive(MatrixXi &H);

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

    //Finds the element nearest x FROM BELOW in a C-style array
    template <typename List_type, typename T, typename size_type>
    int upper_bound(const List_type &list , const T& x, const size_type &size){
        //Now find the point in list closest to x, from below.
        if (size <= 1){return 0;}
        auto low  = std::upper_bound (list, list + size, x);
        return  low-list-1;
    }

}

#endif //WL_MATH_ALGORITHMS_H
