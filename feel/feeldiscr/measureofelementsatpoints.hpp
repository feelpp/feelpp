/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013 Feel++ Consortium

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
   \file measureofelementsatpoints.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_MEASURE_OF_ELEMENTS_AT_POINTS_HPP)

#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

/**
   Given a function space \p Xh, compute the sum of of the measure of the
   elements shared by a point in the mesh. This is helfpul for example for a
   posteriori error estimation.
 */
template<typename MeshType>
auto
measurePointElements( boost::shared_ptr<FunctionSpace<MeshType,bases<Lagrange<MeshType::nOrder,Scalar> > > >& Xh ) -> decltype( Xh->element() )
{
    auto _fn = Xh->element( "measurePointElements" );
    _fn.setZero();
    std::vector<bool> ptdone( Xh->mesh()->numPoints(), false );
    auto elit = Xh->mesh()->beginElement();
    auto elen = Xh->mesh()->endElement();

    for ( ; elit != elen; ++ elit )
    {
        for ( int p = 0; p < elit->numPoints; ++p )
        {
            if ( ptdone[elit->point( p ).id()] == false )
            {
                BOOST_FOREACH( auto pt, elit->point( p ).elements() )
                {
                    _fn.plus_assign( elit->id(), p, 0, elit->measure() );
                }
                ptdone[elit->point( p ).id()] = true;
            }
        }
    }

    return _fn;
}

}

#endif /* FEELPP_MEASURE_OF_ELEMENTS_AT_POINTS_HPP */
