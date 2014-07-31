/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
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
/**
   \file integrate.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2014-01-04
 */
#ifndef FEELPP_VF_INTEGRATE_HPP
#define FEELPP_VF_INTEGRATE_HPP 1

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/integrator.hpp>

namespace Feel {

BOOST_PARAMETER_FUNCTION(
    ( typename vf::detail::integrate_type<Args>::expr_type ), // return type
    integrate,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( range, *  )
      ( expr,   * )
    ) // 4. one required parameter, and

    ( optional
      ( quad,   *, typename vf::detail::integrate_type<Args>::_quad_type() )
      ( geomap, *, (vf::detail::integrate_type<Args>::geoOrder > 1 )?GeomapStrategyType::GEOMAP_OPT:GeomapStrategyType::GEOMAP_HO )
      ( quad1,   *, typename vf::detail::integrate_type<Args>::_quad1_type() )
      ( use_tbb,   ( bool ), false )
      ( use_harts,   ( bool ), false )
      ( grainsize,   ( int ), 100 )
      ( partitioner,   *, "auto" )
      ( verbose,   ( bool ), false )
      ( quadptloc, *, typename vf::detail::integrate_type<Args>::_quadptloc_ptrtype() )
    )
)
{

    auto ret =  integrate_impl( range, quad, expr, geomap, quad1, use_tbb, use_harts, grainsize, partitioner, quadptloc );

    if ( verbose )
    {
        std::cout << " -- integrate: size(range) = " << std::distance( ret.expression().beginElement(),
                  ret.expression().endElement() ) << "\n";
        std::cout << " -- integrate: quad = " << ret.expression().im().nPoints() << "\n";
        std::cout << " -- integrate: quad1 = " << ret.expression().im2().nPoints() << "\n";
        //std::cout << " -- integrate: geomap = " << geomap << "\n";
    }

    return ret;
}

}

#endif
