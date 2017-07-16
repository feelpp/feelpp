/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-*/

#include <feel/feelopt/nlopt.hpp>


namespace Feel
{
namespace opt
{

double feel_nlopt_objective_function(unsigned n, const double *x, double *grad, void *f_data)
{
    OptimizationNonLinear* theopt = static_cast<OptimizationNonLinear*> ( f_data );
    CHECK ( theopt ) << "bad static_cast of OptimizationNonLinear";
    return theopt->objective_function( n, x, grad, theopt->objective_data );
}


double feel_nlopt_constraint_function(unsigned n, const double *x, double *grad, void *f_data)
{
    std::tuple< int,void*> * constraintConfig = static_cast< std::tuple< int,void*> * >( f_data );
    CHECK( constraintConfig ) << "bad cast of constraintConfig";
    int constraintId = std::get<0>( *constraintConfig );
    OptimizationNonLinear* theopt = static_cast<OptimizationNonLinear*> ( std::get<1>( *constraintConfig ) );
    CHECK ( theopt ) << "bad static_cast of OptimizationNonLinear";
    CHECK( constraintId < theopt->constraints_function_data.size() ) << "invalid constraintId";
    OptimizationNonLinear::nlopt_func_type constraintFunction = std::get<0>( (theopt->constraints_function_data)[constraintId] );
    void* constraintData = std::get<1>( (theopt->constraints_function_data)[constraintId] );
    return constraintFunction( n, x, grad, constraintData );
}

double feel_nlopt_objective_vfunction(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
    OptimizationNonLinear* theopt = static_cast<OptimizationNonLinear*> ( f_data );
    CHECK ( theopt ) << "bad static_cast of OptimizationNonLinear";
    return theopt->objective_vfunction( x, grad, theopt->objective_data );
}


double feel_nlopt_constraint_vfunction(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
    std::tuple< int,void*> * constraintConfig = static_cast< std::tuple< int,void*> * >( f_data );
    CHECK( constraintConfig ) << "bad cast of constraintConfig";
    int constraintId = std::get<0>( *constraintConfig );
    OptimizationNonLinear* theopt = static_cast<OptimizationNonLinear*> ( std::get<1>( *constraintConfig ) );
    CHECK ( theopt ) << "bad static_cast of OptimizationNonLinear";
    CHECK( constraintId < theopt->constraints_vfunction_data.size() ) << "invalid constraintId";
    OptimizationNonLinear::nlopt_vfunc_type constraintFunction = std::get<0>( (theopt->constraints_vfunction_data)[constraintId] );
    void* constraintData = std::get<1>( (theopt->constraints_vfunction_data)[constraintId] );
    return constraintFunction( x, grad, constraintData );
}

} // namespace opt
} // namespace Feel
