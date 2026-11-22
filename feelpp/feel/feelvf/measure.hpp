/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2014-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_VF_MEASURE_HPP
#define FEELPP_VF_MEASURE_HPP 1

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/integrate.hpp>

namespace Feel {

template <typename ... Ts>
double measure( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    auto && expr = args.get_else(_expr,cst(1.));
    bool parallel = args.get_else( _parallel, true );
    auto && quad = args.get_else( _quad, quad_order_from_expression );
    auto && quad1 = args.get_else( _quad1, quad_order_from_expression );
    GeomapStrategyType geomap = args.get_else( _geomap,GeomapStrategyType::GEOMAP_OPT );
    bool use_harts = args.get_else( _use_harts, false );
    bool verbose = args.get_else( _verbose, false );

    double meas = integrate( _range=range, _expr=expr, _quad=quad, _quad1=quad1, _geomap=geomap,
                             _use_harts=use_harts, _verbose=verbose ).evaluate(parallel)( 0, 0 );
    DLOG(INFO) << "[mean] measure = " << meas << "\n";
    if ( math::abs(meas) < 1e-13 ) 
      LOG(WARNING) << "Invalid domain measure : " << meas << ", domain range: " << nelements( range ) << "\n";
    return meas;
}

}

#endif /* FEELPP_VF_MEASURE_HPP */
