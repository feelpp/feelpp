/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

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
/**
   \file test_nlopt.cpp
   \author Guillaume Dolle <gdolle@unistra.fr>
   \date 2014-07-15
   Test based on Ipopt tutorial (authors: Carl Laird, Andreas Waechter IBM 2004-11-05)
 */
#include <iostream>
#include <cassert>
#include <feel/feelconfig.h>
#include <feel/feelcore/environment.hpp>

#if defined( FEELPP_HAS_IPOPT )

#include <coin/IpTNLP.hpp>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>

using namespace Ipopt;

/** C++ Example NLP for interfacing a problem with IPOPT.
 *
 * min_x f(x) = -(x2-2)^2
 *  s.t.
 *       0 = x1^2 + x2 - 1
 *       -1 <= x1 <= 1
 *
 */
class MyNLP : public TNLP
{
    public:
        /** default constructor */
        MyNLP(){};

        /** default destructor */
        virtual ~MyNLP(){};

        /**@name Overloaded from TNLP */
        //@{
        /** Method to return some info about the nlp */
        virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                Index& nnz_h_lag, IndexStyleEnum& index_style)
        {
            n = 2;
            m = 1;
            nnz_jac_g = 2;
            nnz_h_lag = 2;
            index_style = FORTRAN_STYLE;
            return true;
        }

        /** Method to return the bounds for my problem */
        virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                Index m, Number* g_l, Number* g_u)
        {
            assert(n == 2);
            assert(m == 1);
            // upper and lower bounds
            x_l[0] = -1.0;
            x_u[0] = 1.0;
            x_l[1] = -1.0e19;
            x_u[1] = +1.0e19;
            g_l[0] = g_u[0] = 0.0;

            return true;
        }

        /** Method to return the starting point for the algorithm */
        virtual bool get_starting_point(Index n, bool init_x, Number* x,
                bool init_z, Number* z_L, Number* z_U,
                Index m, bool init_lambda,
                Number* lambda)
        {
            assert(init_x == true);
            assert(init_z == false);
            assert(init_lambda == false);

            // we initialize x in bounds, in the upper right quadrant
            x[0] = 0.5;
            x[1] = 1.5;

            return true;
        }

        /** Method to return the objective value */
        virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
        {
            Number x2 = x[1];
            obj_value = -(x2 - 2.0) * (x2 - 2.0);
            return true;
        }

        /** Method to return the gradient of the objective */
        virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
        {
            grad_f[0] = 0.0;
            Number x2 = x[1];
            grad_f[1] = -2.0*(x2 - 2.0);
            return true;
        }

        /** Method to return the constraint residuals */
        virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
        {
            Number x1 = x[0];
            Number x2 = x[1];
            g[0] = -(x1*x1 + x2 - 1.0);
            return true;
        }

        /** Method to return:
         *   1) The structure of the jacobian (if "values" is NULL)
         *   2) The values of the jacobian (if "values" is not NULL)
         */
        virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                Index m, Index nele_jac, Index* iRow, Index *jCol,
                Number* values)
        {
            if (values == NULL) {
                iRow[0] = 1;
                jCol[0] = 1;
                iRow[1] = 1;
                jCol[1] = 2;
            }
            else {
                Number x1 = x[0];
                values[0] = -2.0 * x1;
                values[1] = -1.0;
            }

            return true;
        }

        /** Method to return:
         *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
         *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
         */
        virtual bool eval_h(Index n, const Number* x, bool new_x,
                Number obj_factor, Index m, const Number* lambda,
                bool new_lambda, Index nele_hess, Index* iRow,
                Index* jCol, Number* values)
        {
            if (values == NULL)
            {
                iRow[0] = 1;
                jCol[0] = 1;
                iRow[1] = 2;
                jCol[1] = 2;
            }
            else
            {
                values[0] = -2.0 * lambda[0];
                values[1] = -2.0 * obj_factor;
            }
            return true;
        }

        virtual void finalize_solution(SolverReturn status,
                Index n, const Number* x, const Number* z_L, const Number* z_U,
                Index m, const Number* g, const Number* lambda,
                Number obj_value,
                const IpoptData* ip_data,
                IpoptCalculatedQuantities* ip_cq)
        {
        }

    private:
        MyNLP(const MyNLP&);
        MyNLP& operator=(const MyNLP&);
};


int main(int argc, char* argv[])
{
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="test_ipopt",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    SmartPtr<TNLP> mynlp = new MyNLP();
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    auto status = app->Initialize();
    if (status != Solve_Succeeded)
    {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int) status;
    }

    status = app->OptimizeTNLP(mynlp);

    if (status == Solve_Succeeded)
    {
        // Retrieve some statistics about the solve
        Index iter_count = app->Statistics()->IterationCount();
        std::cout << std::endl << std::endl << "*** The problem solved in " << iter_count << " iterations!" << std::endl;

        Number final_obj = app->Statistics()->FinalObjective();
        std::cout << std::endl << std::endl << "*** The final value of the objective function is " << final_obj << '.' << std::endl;
    }
    return 0;
}

#endif
