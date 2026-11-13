/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2014-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_VF_MEAN_H
#define FEELPP_VF_MEAN_H

#include <feel/feelvf/field.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/integrator.hpp>
#include <feel/feelvf/form.hpp>

namespace Feel {

template <typename ... Ts>
auto mean( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    auto expr = Feel::vf::toValueExpression( args.get(_expr) );
    bool parallel = args.get_else( _parallel, parallelEvaluation );
    auto && quad = args.get_else( _quad, quad_order_from_expression );
    auto && quad1 = args.get_else( _quad1, quad_order_from_expression );
    GeomapStrategyType geomap = args.get_else( _geomap,GeomapStrategyType::GEOMAP_OPT );
    bool use_harts = args.get_else( _use_harts, false );
    bool verbose = args.get_else( _verbose, false );
    worldcomm_ptr_t worldcomm = args.get_else( _worldcomm, Environment::worldCommPtr() );

    double meas = integrate( _range=range, _expr=cst(1.0), _quad=quad, _quad1=quad1, _geomap=geomap,
                             _use_harts=use_harts, _verbose=verbose ).evaluate( parallel )( 0, 0 );
    DLOG(INFO) << "[mean] nelements = " << nelements(range) << "\n";
    DLOG(INFO) << "[mean] measure = " << meas << "\n";
    CHECK( math::abs(meas) > 1e-16 ) << "Invalid domain measure : " << meas << ", domain range: " << nelements( range ) << "\n";
    auto eint = integrate( _range=range, _expr=expr, _quad=quad, _geomap=geomap,
                           _quad1=quad1, _use_harts=use_harts, _verbose=verbose ).evaluate( parallel );
    DLOG(INFO) << "[mean] integral = " << eint << "\n";
    auto res = (eint/meas).eval();
    DLOG(INFO) << "[mean] mean = " << res << "\n";
    return res;
}

template <typename ... Ts>
auto form1_mean( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && test = args.get(_test);
    auto && range = args.get(_range);
    auto expr = Feel::vf::toValueExpression( args.get(_expr) );
    bool parallel = args.get_else( _parallel, true );
    auto && quad = args.get_else( _quad, quad_order_from_expression );
    auto && quad1 = args.get_else( _quad1, quad_order_from_expression );
    GeomapStrategyType geomap = args.get_else( _geomap,GeomapStrategyType::GEOMAP_OPT );
    bool use_harts = args.get_else( _use_harts, false );
    bool verbose = args.get_else( _verbose, false );
    worldcomm_ptr_t worldcomm = args.get_else( _worldcomm, Environment::worldCommPtr() );

    double meas = integrate( _range=range, _expr=cst(1.0), _quad=quad, _quad1=quad1, _geomap=geomap,
                             _use_harts=use_harts, _verbose=verbose ).evaluate( parallel )( 0, 0 );
    DLOG(INFO) << "[mean] nelements = " << nelements(range) << "\n";
    DLOG(INFO) << "[mean] measure = " << meas << "\n";
    CHECK( math::abs(meas) > 1e-16 ) << "Invalid domain measure : " << meas << ", domain range: " << nelements( range ) << "\n";
    auto l = form1( _test=test );
    l = integrate( _range=range, _expr=expr/meas, _quad=quad, _geomap=geomap,
                   _quad1=quad1, _use_harts=use_harts, _verbose=verbose );
    return l;
}

}
#endif // FEELPP_VF_MEAN_HPP
