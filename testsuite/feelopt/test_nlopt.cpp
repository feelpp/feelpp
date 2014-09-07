/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-07-15

  Copyright (C) 2014 Feel++ Consortium

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
/**
   \file test_nlopt.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-07-15
 */
#include <iostream>

#include <feel/feelconfig.h>

#if defined( FEELPP_HAS_NLOPT )

#include <nlopt.hpp>

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / std::sqrt(x[1]);
    }
    return std::sqrt(x[1]);
}
typedef struct {
    double a, b;
} my_constraint_data;

double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
 }

int main(int argc, char** argv)
{
    nlopt::opt opt(nlopt::LD_MMA, 2);

    std::vector<double> lb(2);
    lb[0] = -HUGE_VAL; lb[1] = 0;
    opt.set_lower_bounds(lb);

    opt.set_min_objective(myfunc, NULL);

    my_constraint_data data[2] = { {2,0}, {-1,1} };
    opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
    opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);

    opt.set_xtol_rel(1e-4);

    std::vector<double> x(2);
    x[0] = 1.234; x[1] = 5.678;
    double minf;

    nlopt::result result = opt.optimize(x, minf);
    if (result < 0) {
        std::cout << "nlopt failed!\n";
    }
    else {
        std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = " << minf << "\n";
    }
}
#endif
