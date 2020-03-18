/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 12 Mar 2020

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

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feells/reinit_fms.hpp>

namespace Feel {

/**
 * compute the distance field as element of \p Xh  to boundary defined by range \p r
 * \code
 * // distance to boundary faces
 * auto phi = distToEntityRange( Xh, boundaryfaces( Xh->mesh() ) );
 * // distance to faces marked moving
 * auto phi = distToEntityRange( Xh, markedfaces( Xh->mesh(), "moving" ) );
 * \endcode
 * \note Xh must be Lagrange P1
 */
template<typename SpaceType, typename RangeType, typename = std::enable_if_t<is_functionspace_v<SpaceType>>>
element_t<SpaceType> distToEntityRange( std::shared_ptr<SpaceType> const& Xh, RangeType const& r )
{
    auto thefms = fms( Xh );
    auto phio = Xh->element();
    phio.on( _range=elements(Xh->mesh()), _expr=h());
    //phio.on( _range=markedfaces(Xh->mesh(),r), _expr= -h()/100. );
    phio.on( _range=r, _expr= -h()/100. );
    return thefms->march(phio);
}

}
