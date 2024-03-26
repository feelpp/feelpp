/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-30

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
/**
   \file dh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-01-30
 */
#ifndef FEELPP_DH_H
#define FEELPP_DH_H 1

#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

template<int Order,typename MeshType,typename T = double>
using dh_type = FunctionSpace<MeshType,bases<RaviartThomas<Order>>,T,Periodicity <NoPeriodicity>>;

template<int Order,typename MeshType>
using dh_ptrtype = std::shared_ptr<dh_type<Order,MeshType>>;


template<typename MeshType, int Order,typename T = double>
using Dh_type = FunctionSpace<MeshType,bases<RaviartThomas<Order>>,T,Periodicity <NoPeriodicity>>;

template<typename MeshType, int Order>
using Dh_ptrtype = std::shared_ptr<dh_type<Order,MeshType>>;

/**
 * \fn Dh<k,MeshType>
 *
 * build a function space of continuous function which are piecewise polynomial
 * of degree (total or in each variable) less than k.
 */
template<int Order,typename MeshType>
inline
dh_ptrtype<Order,MeshType>
Dh( std::shared_ptr<MeshType> mesh, bool buildExtendedDofTable=false )
{
    return dh_type<Order,MeshType>::New( _mesh=mesh,
                                         _worldscomm=makeWorldsComm( 1, mesh->worldComm() ),
                                         _extended_doftable=std::vector<bool>( 1,buildExtendedDofTable ) );
}


} // Feel

#endif /* FEELPP_DH_H */
