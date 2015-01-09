/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-04

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
#ifndef FEELPP_VF_MEAN_HPP
#define FEELPP_VF_MEAN_HPP 1

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/integrator.hpp>


namespace Feel {

BOOST_PARAMETER_FUNCTION(
    ( typename vf::detail::integrate_type<Args>::expr_type::expression_type::matrix_type ), // return type
    mean,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( range, *  )
      ( expr,   * )
    ) // 4. one required parameter, and

    ( optional
      ( quad,   *, typename vf::detail::integrate_type<Args>::_quad_type() )
      ( geomap, *, GeomapStrategyType::GEOMAP_OPT )
      ( quad1,   *, typename vf::detail::integrate_type<Args>::_quad1_type() )
      ( use_tbb,   ( bool ), false )
      ( use_harts,   ( bool ), false )
      ( grainsize,   ( int ), 100 )
      ( partitioner,   *, "auto" )
      ( verbose,   ( bool ), false )
      ( parallel,   *( boost::is_integral<mpl::_> ), 1 )
      ( worldcomm,       (WorldComm), Environment::worldComm() )
    )
)
{
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

}
#endif // FEELPP_VF_MEAN_HPP
