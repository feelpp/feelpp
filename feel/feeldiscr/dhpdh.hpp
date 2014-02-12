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
#if !defined(FEELPP_DHPDH_HPP)
#define FEELPP_DHPDH_HPP 1

#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

/**

   \code
   #include <feel/feeldiscr/dhpdh.hpp>
   auto Xh = DhPdh<2>( mesh );
   \endcode
 */
template<int Order,typename MeshType>
inline
boost::shared_ptr<FunctionSpace<MeshType,bases<RaviartThomas<Order>,Lagrange<Order,Scalar,Discontinuous> > > >
DhPdh( boost::shared_ptr<MeshType> mesh,
       std::vector<bool> buildExtendedDofTable = std::vector<bool>( 2,false ) )
{
    CHECK( buildExtendedDofTable.size() == 2 ) << " vector activation for extended dof table must be equal to 2 but here " << buildExtendedDofTable.size() << "\n";
    return FunctionSpace<MeshType,bases<RaviartThomas<Order>,Lagrange<Order,Scalar,Discontinuous>>>::New( _mesh=mesh,
                                                                                                          _worldscomm=std::vector<WorldComm>( 2,mesh->worldComm() ),
                                                                                                          _extended_doftable=buildExtendedDofTable );
}


}
#endif /* FEELPP_DHPDH_HPP */
