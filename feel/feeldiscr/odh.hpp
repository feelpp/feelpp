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
   \file odh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_ODH_HPP)
#define FEELPP_ODH_HPP 1

#include <feel/feelpoly/orthonormalpolynomialset.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

namespace meta {
template<typename MeshType,
         int Order,
         int Tag = 0>
struct Odh
{
    typedef FunctionSpace<MeshType,
                          bases<OrthonormalPolynomialSet<Order,Scalar>>,
                          double,
                          Periodicity <NoPeriodicity>,
                          mortars<NoMortar>> type;
    typedef boost::shared_ptr<type> ptrtype;
};

} // meta

/**
   Given a \p mesh, build a function space of discontinuous function which are
  piecewise polynomial of degree and \f$L_2\f$ Orthonormal
 */
template<int Order,typename MeshType, int Tag = 0>
inline
typename meta::Odh<MeshType,Order,Tag>::ptrtype
Odh( boost::shared_ptr<MeshType> mesh, bool buildExtendedDofTable=false )
{
    typedef typename meta::Odh<MeshType,Order,Tag>::type space_type;
    return space_type::New( _mesh=mesh,
                            _worldscomm=std::vector<WorldComm>( 1,mesh->worldComm() ),
                            _extended_doftable=std::vector<bool>( 1,buildExtendedDofTable ) );
}

}
#endif /* FEELPP_ODH_HPP */
