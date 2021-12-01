/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-04

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
#ifndef FEELPP_VF_MEAN_H
#define FEELPP_VF_MEAN_H

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/integrator.hpp>
#include <feel/feelvf/form.hpp>

namespace Feel {

template <typename ... Ts>
auto mean( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    auto && expr = args.get(_expr);
    bool parallel = args.get_else( _parallel, true );
    auto && quad = args.get_else( _quad, quad_order_from_expression );
    auto && quad1 = args.get_else( _quad1, quad_order_from_expression );
    GeomapStrategyType geomap = args.get_else( _geomap,GeomapStrategyType::GEOMAP_OPT );
    bool use_tbb = args.get_else( _use_tbb, false );
    bool use_harts = args.get_else( _use_harts, false );
    int grainsize = args.get_else( _grainsize, 100 );
    std::string const& partitioner = args.get_else( _partitioner, "auto");
    bool verbose = args.get_else( _verbose, false );
    worldcomm_ptr_t worldcomm = args.get_else( _worldcomm, Environment::worldCommPtr() );

    double meas = integrate( _range=range, _expr=cst(1.0), _quad=quad, _quad1=quad1, _geomap=geomap,
                             _use_tbb=use_tbb, _use_harts=use_harts, _grainsize=grainsize,
                             _partitioner=partitioner, _verbose=verbose ).evaluate( parallel,worldcomm )( 0, 0 );
    DLOG(INFO) << "[mean] nelements = " << nelements(range) << "\n";
    DLOG(INFO) << "[mean] measure = " << meas << "\n";
    CHECK( math::abs(meas) > 1e-16 ) << "Invalid domain measure : " << meas << ", domain range: " << nelements( range ) << "\n";
    auto eint = integrate( _range=range, _expr=expr, _quad=quad, _geomap=geomap,
                           _quad1=quad1, _use_tbb=use_tbb, _use_harts=use_harts, _grainsize=grainsize,
                           _partitioner=partitioner, _verbose=verbose ).evaluate( parallel,worldcomm );
    DLOG(INFO) << "[mean] integral = " << eint << "\n";
    DLOG(INFO) << "[mean] mean = " << eint/meas << "\n";
    return eint/meas;
}

template <typename ... Ts>
auto form1_mean( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && test = args.get(_test);
    auto && range = args.get(_range);
    auto && expr = args.get(_expr);
    bool parallel = args.get_else( _parallel, true );
    auto && quad = args.get_else( _quad, quad_order_from_expression );
    auto && quad1 = args.get_else( _quad1, quad_order_from_expression );
    GeomapStrategyType geomap = args.get_else( _geomap,GeomapStrategyType::GEOMAP_OPT );
    bool use_tbb = args.get_else( _use_tbb, false );
    bool use_harts = args.get_else( _use_harts, false );
    int grainsize = args.get_else( _grainsize, 100 );
    std::string const& partitioner = args.get_else( _partitioner, "auto");
    bool verbose = args.get_else( _verbose, false );
    worldcomm_ptr_t worldcomm = args.get_else( _worldcomm, Environment::worldCommPtr() );

    double meas = integrate( _range=range, _expr=cst(1.0), _quad=quad, _quad1=quad1, _geomap=geomap,
                             _use_tbb=use_tbb, _use_harts=use_harts, _grainsize=grainsize,
                             _partitioner=partitioner, _verbose=verbose ).evaluate( parallel,worldcomm )( 0, 0 );
    DLOG(INFO) << "[mean] nelements = " << nelements(range) << "\n";
    DLOG(INFO) << "[mean] measure = " << meas << "\n";
    CHECK( math::abs(meas) > 1e-16 ) << "Invalid domain measure : " << meas << ", domain range: " << nelements( range ) << "\n";
    auto l = form1( _test=test );
    l = integrate( _range=range, _expr=expr/meas, _quad=quad, _geomap=geomap,
                   _quad1=quad1, _use_tbb=use_tbb, _use_harts=use_harts, _grainsize=grainsize,
                   _partitioner=partitioner, _verbose=verbose );
    return l;
}

}
#endif // FEELPP_VF_MEAN_HPP
