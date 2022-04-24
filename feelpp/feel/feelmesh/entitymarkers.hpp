/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 26 Nov 2019

 Copyright (C) 2019 Feel++ Consortium

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
#ifndef FEELPP_ENTITYMARKERS_HPP
#define FEELPP_ENTITYMARKERS_HPP 1

#include <optional>
#include <map>
#include <feel/feelcore/feel.hpp>


namespace Feel {

/**
 * get the map of markers of a entity of a geoelement
 * \param e element id to get face markers from
 * \param f entity id in the element
 * \return an optional map of markers for the entity \p f of the element \p e
 */
template<int EntityCoDim, typename ElementT>
inline
std::optional<std::map<uint16_type,typename ElementT::marker_type/*Marker1*/>> entityMarkers( ElementT const& e, uint16_type f )
{
    if constexpr ( dimension_v<ElementT> >= 1  && is_geoelement_v<ElementT> )
    {
        if constexpr ( EntityCoDim == 1 )
        {
            if ( f != invalid_v<uint16_type> && e.hasFace(f) )
                return e.face(f).markers();
        }
        else if  constexpr ( dimension_v<ElementT> == 3 &&  EntityCoDim == 2 )
        {
            if ( f != invalid_v<uint16_type> && e.hasEdge(f) )
                return e.edge(f).markers();
        }
        else if  constexpr ( dimension_v<ElementT> >= 2 &&  EntityCoDim == dimension_v<ElementT>  )
        {
            if ( f != invalid_v<uint16_type> )
                return e.point(f).markers();
        }
    }
    return std::nullopt;
}


}
#endif
