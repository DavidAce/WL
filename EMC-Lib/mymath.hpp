#ifndef MYMATH_H   
#define MYMATH_H  
#include <Eigen/Dense>
#include <Eigen/Core>
#include "randomFunctions.hpp"

using namespace EMC_rnd;
using namespace Eigen;
//extern void linspace(ArrayXd &,double, double, int);
extern ArrayXd LogSpaced(int, double, double);
extern void rndChoice(int [], int , int);

extern int get_sign(double);
extern int tri_randint(int, int );
extern int tri_inv_randint(int, int );
extern int mod(int , int );

template <typename Derived, typename searchType>
bool isvalueinarray(ArrayBase<Derived> &arr, searchType val){
    for (int i = 0; i < arr.size(); i++) {
        if (arr(i) == val)
            return true;
    }
    return false;
};
extern int heaviside(double );
extern int find_index(int [], int, int);
 
extern double mean(double[], int );
extern double var(double [], int);


#endif