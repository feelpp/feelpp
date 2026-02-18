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
   \file bdmh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-02-18
 */
#ifndef FEELPP_BDMH_H
#define FEELPP_BDMH_H 1

#include <feel/feelpoly/brezzidouglasmarini.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/order.hpp>

namespace Feel {

template<int Order, typename MeshType, typename T = double>
using bdmh_type = FunctionSpace<MeshType, bases<BrezziDouglasMarini<Order>>, T, Periodicity<NoPeriodicity>>;

template<int Order, typename MeshType>
using bdmh_ptrtype = std::shared_ptr<bdmh_type<Order, MeshType>>;

template<typename MeshType, int Order, typename T = double>
using BDMh_type = FunctionSpace<MeshType, bases<BrezziDouglasMarini<Order>>, T, Periodicity<NoPeriodicity>>;

template<typename MeshType, int Order>
using BDMh_ptrtype = std::shared_ptr<bdmh_type<Order, MeshType>>;

template<int Order, typename MeshType>
inline
bdmh_ptrtype<Order, MeshType>
BDMh( std::shared_ptr<MeshType> const& mesh,
      RuntimeOrder order,
      DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return bdmh_type<Order, MeshType>::New( _mesh = mesh,
                                            _worldscomm = makeWorldsComm( 1, mesh->worldComm() ),
                                            _extended_doftable = dte,
                                            _runtime_order = order );
}

template<int Order, typename MeshType>
    requires ( Order >= 0 )
inline
bdmh_ptrtype<Order, MeshType>
BDMh( std::shared_ptr<MeshType> const& mesh,
      DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return BDMh<Order>( mesh, RuntimeOrder{ static_cast<uint16_type>( Order ) }, dte );
}

} // namespace Feel

#endif /* FEELPP_BDMH_H */
