//
// Created by david on 2016-08-05.
//

#include "math_algorithms.h"

namespace math{
    int find_max_idx(const Ref<const VectorXd>  &list) {
        //int find_max_idx(const Ref<const Array<double,Dynamic,1>> &list){
        int i, idx = 0;
        double peak = list(0);
        for (i = 0; i < list.size(); i++) {
            if (list(i) > peak) {
                peak = list(i);
                idx = i;
            }
        }
        return idx;
    }

    int find_min_positive(MatrixXi &H) {
        int min = 1000000000;
        for (int j = 0; j < H.cols(); j++) {
            for (int i = 0; i < H.rows(); i++) {
                if (H(i, j) < min) {
                    if (H(i, j) > 0) {
                        min = H(i, j);
                    }
                }
            }
        }
        return min;
    }
}