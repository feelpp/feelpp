/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-30

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
   \file ned1h.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-01-30
 */
#ifndef FEELPP_NED1H_H
#define FEELPP_NED1H_H 1

#include <feel/feelpoly/nedelec.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

namespace meta {
template<typename MeshType,int Order>
struct Ned1h
{
    typedef FunctionSpace<MeshType,
                          bases<Nedelec<Order,NedelecKind::NED1>>,
                          double,
                          Periodicity <NoPeriodicity>,
                          mortars<NoMortar>> type;
    typedef boost::shared_ptr<type> ptrtype;
};

} // meta

/**
 * \fn Ned1h<k,MeshType>
 *
 */
template<int Order,typename MeshType>
inline
typename meta::Ned1h<MeshType,Order>::ptrtype
Ned1h( boost::shared_ptr<MeshType> mesh, bool buildExtendedDofTable=false )
{
    typedef typename meta::Ned1h<MeshType,Order>::type space_type;
    return space_type::New( _mesh=mesh,
                            _worldscomm=worldsComm( mesh->worldComm() ),
                            _extended_doftable=std::vector<bool>( 1,buildExtendedDofTable ) );
}


} // Feel

#endif /* FEELPP_NED1H_H */
