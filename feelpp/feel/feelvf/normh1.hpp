/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2014-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file normh1.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-01-04
 */
#ifndef FEELPP_VF_NORMH1_H
#define FEELPP_VF_NORMH1_H

#include <feel/feelvf/field.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/integrator.hpp>

namespace Feel  {

template <typename ... Ts>
double normH1( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    auto expr = Feel::vf::toValueExpression( args.get(_expr) );
    auto && grad_expr = args.get(_grad_expr);
    bool parallel = args.get_else( _parallel, true );
    auto && quad = args.get_else( _quad, quad_order_from_expression );
    auto && quad1 = args.get_else( _quad1, quad_order_from_expression );
    GeomapStrategyType geomap = args.get_else( _geomap,GeomapStrategyType::GEOMAP_OPT );
    bool use_harts = args.get_else( _use_harts, false );
    bool verbose = args.get_else( _verbose, false );

    double a = integrate( _range=range, _expr=inner(expr), _quad=quad, _geomap=geomap,
                          _quad1=quad1, _use_harts=use_harts, _verbose=verbose ).evaluate(parallel)( 0, 0 );
    double b = integrate( _range=range, _expr=inner(grad_expr), _quad=quad, _geomap=geomap,
                          _quad1=quad1, _use_harts=use_harts, _verbose=verbose ).evaluate(parallel)( 0, 0 );
    return math::sqrt( a + b );
}

}

#endif
