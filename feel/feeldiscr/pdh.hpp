/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013-2016 Feel++ Consortium

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
   \file pdh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined( FEELPP_PDH_HPP )
#define FEELPP_PDH_HPP

#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{

template <typename MeshType, int Order, template <class, uint16_type, class> class Pts = PointSetEquiSpaced>
using Pdh_type = FunctionSpace<MeshType, bases<Lagrange<Order, Scalar, Discontinuous, Pts>>>;
template <typename MeshType, int Order, template <class, uint16_type, class> class Pts = PointSetEquiSpaced>
using Pdh_ptrtype = boost::shared_ptr<Pdh_type<MeshType, Order, Pts>>;

/**
   Given a \p mesh, build a function space of discontinuous function which are
   piecewise polynomial of degree (total or in each variable) less than k.
*/
template <int Order, template <class, uint16_type, class> class Pts = PointSetEquiSpaced, typename MeshType>
inline Pdh_ptrtype<MeshType, Order, Pts>
Pdh( boost::shared_ptr<MeshType> mesh, bool buildExtendedDofTable = false )
{
    return Pdh_type<MeshType, Order, Pts>::New( _mesh = mesh,
                                                _worldscomm = std::vector<WorldComm>( 1, mesh->worldComm() ),
                                                _extended_doftable = std::vector<bool>( 1, buildExtendedDofTable ) );
}
}

#endif /* FEELPP_PDH_HPP */
