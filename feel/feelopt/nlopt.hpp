/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-07-15

  Copyright (C) 2014-2016 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef FEELPP_NLOPT_HPP
#define FEELPP_NLOPT_HPP 1

#include <feel/feelconfig.h>

#if defined(FEELPP_HAS_NLOPT)
#include <nlopt.hpp>
#include <functional>
#include <iostream>
#include <feel/feelcore/feel.hpp>

namespace Feel
{
namespace opt
{
double feel_nlopt_objective_function(unsigned n, const double *x, double *grad, void *f_data);
double feel_nlopt_constraint_function(unsigned n, const double *x, double *grad, void *f_data);
double feel_nlopt_objective_vfunction(const std::vector<double> &x, std::vector<double> &grad, void *f_data);
double feel_nlopt_constraint_vfunction(const std::vector<double> &x, std::vector<double> &grad, void *f_data);

class OptimizationNonLinear : public ::nlopt::opt
{
    typedef ::nlopt::opt super_type;
public :
    typedef std::function<double (unsigned int n, const double *x,
                                  double *gradient, /* NULL if not needed */
                                  void *func_data)> nlopt_func_type;
    typedef std::function<double (const std::vector<double> &x,
                                  std::vector<double> &gradient, /* NULL if not needed */
                                  void *func_data)> nlopt_vfunc_type;

    OptimizationNonLinear( ::nlopt::algorithm a, unsigned n) : super_type(a,n) {}

    void set_min_objective(::nlopt::func f, void *f_data)
        {
            super_type::set_min_objective( f,f_data );
        }
    void set_min_objective( nlopt_func_type f, void *f_data)
        {
            objective_function = f;
            objective_data = f_data;
            super_type::set_min_objective( feel_nlopt_objective_function,this );
        }
    void set_min_objective(::nlopt::vfunc f, void *f_data)
        {
            super_type::set_min_objective( f,f_data );
        }
    void set_min_objective( nlopt_vfunc_type f, void *f_data)
        {
            objective_vfunction = f;
            objective_data = f_data;
            super_type::set_min_objective( feel_nlopt_objective_vfunction,this );
        }
    void set_max_objective(::nlopt::func f, void *f_data)
        {
            super_type::set_max_objective( f,f_data );
        }
    void set_max_objective( nlopt_func_type f, void *f_data)
        {
            objective_function = f;
            objective_data = f_data;
            super_type::set_max_objective( feel_nlopt_objective_function,this );
        }
    void set_max_objective(::nlopt::vfunc f, void *f_data)
        {
            super_type::set_max_objective( f,f_data );
        }
    void set_max_objective( nlopt_vfunc_type f, void *f_data)
        {
            objective_vfunction = f;
            objective_data = f_data;
            super_type::set_max_objective( feel_nlopt_objective_vfunction,this );
        }
    void add_inequality_constraint( ::nlopt::func f, void *f_data, double tol=0)
        {
            super_type::add_inequality_constraint( f,f_data,tol );
        }
    void add_inequality_constraint( nlopt_func_type f, void *f_data, double tol=0)
        {
            int constraintId = constraints_config.size();
            if ( constraintId == 0 )
                constraints_config.reserve(100);
            CHECK( constraintId < 100 ) << "TODO : case were the vector is resized and consequently invalidated previous pointer on vector";

            constraints_config.push_back( std::make_tuple( constraintId,this ) );
            constraints_function_data.push_back( std::make_tuple( f,f_data ) );
            super_type::add_inequality_constraint( feel_nlopt_constraint_function, (&(constraints_config[constraintId])),tol );
        }
    void add_inequality_constraint( ::nlopt::vfunc f, void *f_data, double tol=0)
        {
            super_type::add_inequality_constraint( f,f_data,tol );
        }
    void add_inequality_constraint( nlopt_vfunc_type f, void *f_data, double tol=0)
        {
            int constraintId = vconstraints_config.size();
            if ( constraintId == 0 )
                vconstraints_config.reserve(100);
            CHECK( constraintId < 100 ) << "TODO : case were the vector is resized and consequently invalidated previous pointer on vector";

            vconstraints_config.push_back( std::make_tuple( constraintId,this ) );
            constraints_vfunction_data.push_back( std::make_tuple( f,f_data ) );
            super_type::add_inequality_constraint( feel_nlopt_constraint_vfunction, (&(vconstraints_config[constraintId])),tol );
        }
    void remove_inequality_constraints()
        {
            constraints_config.clear();
            vconstraints_config.clear();
            constraints_function_data.clear();
            constraints_vfunction_data.clear();
            super_type::remove_inequality_constraints();
        }

    nlopt_func_type objective_function;
    nlopt_vfunc_type objective_vfunction;
    void *objective_data;

    std::vector< std::tuple< int,void*> > constraints_config;
    std::vector< std::tuple< int,void*> > vconstraints_config;
    std::vector<std::tuple<nlopt_func_type,void*>> constraints_function_data;
    std::vector<std::tuple<nlopt_vfunc_type,void*>> constraints_vfunction_data;

};

}

namespace nlopt = ::nlopt;
}


#endif // FEELPP_HAS_NLOPT

#endif
