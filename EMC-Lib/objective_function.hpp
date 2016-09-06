//
// Created by david on 2016-09-05.
//

#ifndef OBJECTIVE_FUNCTION_H
#define OBJECTIVE_FUNCTION_H

#include <memory>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "constants.hpp"
using namespace Eigen;
class objective_function{
private:
public:
    template <typename boundaryType, typename... aux_args>
    objective_function(boundaryType lower, boundaryType upper, double tol,aux_args... in)
            : lower_bound(lower),
              upper_bound(upper),
              tolerance(tol),
              aux({in...})
    {
        parameters = (int) lower_bound.size();

    };
    typedef std::function<double(objective_function &, ArrayXd &)> providedType ;
    providedType providedFunction;

    double fitness(ArrayXd &parameters){
        return providedFunction(*this, parameters);
//        return (double)(-1.0 / log(H + EMC_constants::log_param) + EMC_constants::log_const + pow(log( 1/(H + 1)), 2));
    }


    Array<double,Dynamic, Dynamic> lower_bound; //Minimum allowed values for each fitting parameter
    Array<double,Dynamic, Dynamic> upper_bound; //Maximum allowed values for each fitting parameter
    double tolerance;                           //If the fitness saturates within tolerance, the program terminates
    ArrayXd optimum;
    int parameters;
    std::vector<Array<double,Dynamic,Dynamic>>  aux; //Auxiliary data or vectors for use in the fitness function
//    typedef std::function<outputType(objective_function &, inputType &)> fitnessType ;


};


#endif //OBJECTIVE_FUNCTION_H
