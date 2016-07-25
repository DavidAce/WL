//
// Created by david on 2016-07-24.
//

#ifndef WL_CONSTANTS_H
#define WL_CONSTANTS_H
#include <math.h>

namespace constants{

    //Lattice Properties
    static const int d = 2;         //Dimension
    static const int L = 20;        //Linear size
    static const int N = (int)pow(L,d);  //Number of spins/particles

    //DOS and Histogram properties
    static const int rw_dims = 1;      //Dimension of random walks (1D or 2D WL)
    static const int bins = 400;

}


#endif //WL_CONSTANTS_H
