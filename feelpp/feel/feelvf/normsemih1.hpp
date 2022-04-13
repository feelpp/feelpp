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
#ifndef FEELPP_VF_NORMSEMIH1_H
#define FEELPP_VF_NORMSEMIH1_H

namespace Feel {

template <typename ... Ts>
double normSemiH1( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    auto && grad_expr = args.get(_grad_expr);
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
    double a = integrate( _range=range, _expr=inner(grad_expr), _quad=quad, _geomap=geomap,
                          _quad1=quad1, _use_tbb=use_tbb, _use_harts=use_harts, _grainsize=grainsize,
                          _partitioner=partitioner, _verbose=verbose ).evaluate(parallel,worldcomm)( 0, 0 );
    return math::sqrt( a );
}


}

#endif
