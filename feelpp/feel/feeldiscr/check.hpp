/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2020-03-18

  Copyright (C) 2020 Feel++ Consortium

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
#pragma once

#include <feel/feeldiscr/traits.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel {

// solution enum type
enum class solution_t {
    unique,
    up_to_a_constant
};

/*
 * check L2, H1 and semiH1 norms of an function u
 * the check is triggerred only if @p thechecker @c check() member function returns true
 * @param thechecker is a data structure that holds information about the exact solution and its gradient
 * @param u is the function to be checked, it can be scalar or vectorial
 * 
 * @return 0 if ok, 1 otherwise
 */
template<typename CheckerT, typename ElementT, 
         typename = std::enable_if_t<is_scalar_field_v<ElementT> || is_vector_field_v<ElementT> || is_matrix_field_v<ElementT> >>
int check( CheckerT&& thechecker, ElementT const& u, solution_t s = solution_t::unique )
{
    int status = 0;
    // tag::check[]
    if ( thechecker.check() )
    {
        auto mesh = u.functionSpace()->mesh();
        auto Vh = u.functionSpace();
        
        // compute l2 and h1 norm of u-u_h where u=solution
        auto norms = [=]( std::string const& solution ) -> std::map<std::string, double> {
             
            constexpr int dim = dimension_v<decltype(support(Vh))>;
            auto get_sol_ex = [=]( auto const & c ) {
                if constexpr ( is_scalar_field_v<ElementT> ) 
                {
                    return expr( solution );
                }
                else if constexpr ( is_vector_field_v<ElementT> ) 
                {
                    return expr<dim,1>( solution );
                }
                else if constexpr ( is_matrix_field_v<ElementT> ) 
                {
                    return expr<dim,dim>( solution );
                }
            };
            auto sol_ex = get_sol_ex( thechecker ); 
            sol_ex.setParameterValues( thechecker.parameterValues() );

            tic(); 
            double l2_p=1, l2=0;
    
            if ( s == solution_t::unique )
            {
                l2_p = normL2(_range=elements(support(Vh)), _expr=sol_ex );
                l2 = normL2(_range=elements(support(Vh)), _expr=idv(u)-(sol_ex) );
            }
            else if ( s == solution_t::up_to_a_constant )
            {
                if constexpr ( is_scalar_field_v<ElementT> )
                {
                    auto mean_p_exact = mean( _range=elements(support(Vh)), _expr=sol_ex )(0,0);
                    auto mean_p = mean( _range=elements(support(Vh)), _expr=idv(u) )(0,0);
                    l2 = normL2( _range=elements(support(Vh)),
                                 _expr=(sol_ex - cst(mean_p_exact)) - (idv(u) - cst(mean_p)) );
                    l2_p = normL2( _range=elements(support(Vh)), _expr=sol_ex- cst(mean_p_exact) );
                }
            }

            l2_p = (l2_p<1e-10)?1:l2_p; 
            toc("L2 error norm");
            if ( !thechecker.hasGradient() )
                return { { "L2", l2/l2_p } };
            if ( is_matrix_field_v<ElementT> )
                throw std::invalid_argument( "H1 and SemiH1 norms not supported for Matrix Fields" );
            
            auto get_grad_ex = [=]( auto const & c ) {
                if constexpr ( is_scalar_field_v<ElementT> ) 
                {
                    return expr<1,dim>( thechecker.gradient().value() );
                }
                else if constexpr ( is_vector_field_v<ElementT> ) 
                {
                    return expr<dim,dim>( thechecker.gradient().value() );
                }
            };
            auto grad_ex = get_grad_ex( thechecker );
            grad_ex.setParameterValues( thechecker.parameterValues() );
            tic();
            double h1_p = normH1(_range=elements(support(Vh)), _expr=(sol_ex), _grad_expr=grad_ex );
            double h1 = normH1(_range=elements(support(Vh)), _expr=idv(u)-(sol_ex), _grad_expr=gradv(u)-grad_ex );
            toc("H1 error norm");   
            tic();
            double semih1 = normL2(_range=elements(support(Vh)), _expr=gradv(u)-grad_ex );
            double semih1_p = normL2(_range=elements(support(Vh )), _expr=grad_ex );
            toc("semi H1 error norm");
            if ( semih1_p < 1e-10 )
                return { { "L2", l2/l2_p }, {  "H1", h1/h1_p }, {"semih1",semih1 } }; 
            return { { "L2", l2/l2_p }, {  "H1", h1/h1_p }, {"semih1",semih1/semih1_p} }; 
        };

        status = !thechecker.runOnce( norms, rate::hp( mesh->hMax(), Vh->fe()->order() ) );
    }
    // end::check[]
    return status;
}



} // Feel
