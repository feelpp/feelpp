/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-02-18

  Copyright (C) 2026 Feel++ Consortium

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
   \file neh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-02-18
 */
#ifndef FEELPP_NEH_H
#define FEELPP_NEH_H 1

#include <feel/feeldiscr/ned1h.hpp>
#include <utility>

namespace Feel {

template<typename MeshType, int Order, typename T = double>
using Neh_type = Ned1h_type<MeshType, Order, T>;

template<typename MeshType, int Order, typename T = double>
using Neh_ptrtype = Ned1h_ptrtype<MeshType, Order, T>;

template<int Order, typename MeshType, typename T = double>
inline
Neh_ptrtype<MeshType, Order, T>
Neh( std::shared_ptr<MeshType> const& mesh,
     RuntimeOrder order,
     DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Ned1h<Order, MeshType, T>( mesh, order, dte );
}

template<int Order, typename MeshType, typename T = double>
    requires ( Order >= 0 )
inline
Neh_ptrtype<MeshType, Order, T>
Neh( std::shared_ptr<MeshType> const& mesh,
     DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Ned1h<Order, MeshType, T>( mesh, dte );
}

template<int Order, typename MeshType, typename RangeType, typename T = double>
inline
Neh_ptrtype<MeshType, Order, T>
Neh( std::shared_ptr<MeshType> const& mesh,
     RangeType&& rangeElt,
     RuntimeOrder order,
     DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Ned1h<Order, MeshType, RangeType, T>( mesh,
                                                 std::forward<RangeType>( rangeElt ),
                                                 order,
                                                 dte );
}

template<int Order, typename MeshType, typename RangeType, typename T = double>
    requires ( Order >= 0 )
inline
Neh_ptrtype<MeshType, Order, T>
Neh( std::shared_ptr<MeshType> const& mesh,
     RangeType&& rangeElt,
     DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Ned1h<Order, MeshType, RangeType, T>( mesh,
                                                 std::forward<RangeType>( rangeElt ),
                                                 dte );
}

} // namespace Feel

#endif /* FEELPP_NEH_H */
