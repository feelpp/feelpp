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
   \file pch.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#ifndef FEELPP_PCH_H
#define FEELPP_PCH_H 1

#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

namespace meta {
template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
struct Pch
{
    typedef FunctionSpace<MeshType,
                          bases<Lagrange<Order,Scalar,Continuous,Pts,Tag>>,
                          T,
                          Periodicity <NoPeriodicity>,
                          mortars<NoMortar>> type;
    typedef std::shared_ptr<type> ptrtype;
};

} // meta

template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pch_type = typename meta::Pch<MeshType,Order,T,Pts,Tag>::type;
template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pch_ptrtype = typename meta::Pch<MeshType,Order,T,Pts,Tag>::ptrtype;

template<typename MeshType,int Order,typename T = double, template<class, uint16_type, class> class Pts = PointSetFekete, int Tag = 0>
using Pch_element_t=typename Pch_type<MeshType,Order, T,Pts, Tag>::element_type;

template<typename MeshType,int Order,typename T = double,template<class, uint16_type, class> class Pts = PointSetFekete, int Tag = 0>
using Pch_element_type=Pch_element_t<MeshType,Order,T,Pts, Tag>;


/**
 * \fn Pch<k,MeshType>
 *
 * build a function space of continuous function which are piecewise polynomial
 * of degree (total or in each variable) less than k.
 */
template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,
         int Tag = 0>
inline
Pch_ptrtype<MeshType,Order,T,Pts,Tag>
Pch( std::shared_ptr<MeshType> mesh, bool buildExtendedDofTable=false )
{
    return Pch_type<MeshType,Order,T,Pts,Tag>::New( _mesh=mesh,
                                                    _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                                    _extended_doftable=buildExtendedDofTable );
}

/**
 * \fn Pch<k,MeshType>
 *
 * build a function space of continuous function which are piecewise polynomial
 * of degree (total or in each variable) less than k.
 */
template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,
         int Tag = 0>
inline
Pch_ptrtype<MeshType,Order,T,Pts,Tag>
Pch( std::shared_ptr<MeshType> mesh, elements_reference_wrapper_t<MeshType> const& rangeElt, bool buildExtendedDofTable=false, size_type components = 0 )
{
    return Pch_type<MeshType,Order,T,Pts,Tag>::New( _mesh=mesh,
                                                    _range=rangeElt,
                                                    _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                                    _extended_doftable=buildExtendedDofTable,
                                                    _components=components );
}

#if !defined( FEELPP_INSTANTIATE )
extern template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<0,Scalar>>>;
extern template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<1,Scalar>>>;
extern template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<2,Scalar>>>;
extern template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<3,Scalar>>>;
extern template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<0,Scalar>>>;
extern template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<1,Scalar>>>;
extern template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<2,Scalar>>>;
extern template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<3,Scalar>>>;

extern template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<0,Scalar>>,double, Periodicity <NoPeriodicity>,mortars<NoMortar>>;
extern template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<1,Scalar>>,double, Periodicity <NoPeriodicity>,mortars<NoMortar>>;
extern template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<2,Scalar>>,double, Periodicity <NoPeriodicity>,mortars<NoMortar>>;
extern template class FunctionSpace<Mesh<Simplex<2>>,bases<Lagrange<3,Scalar>>,double, Periodicity <NoPeriodicity>,mortars<NoMortar>>;
extern template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<0,Scalar>>,double, Periodicity <NoPeriodicity>,mortars<NoMortar>>;
extern template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<1,Scalar>>,double, Periodicity <NoPeriodicity>,mortars<NoMortar>>;
extern template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<2,Scalar>>,double, Periodicity <NoPeriodicity>,mortars<NoMortar>>;
extern template class FunctionSpace<Mesh<Simplex<3>>,bases<Lagrange<3,Scalar>>,double, Periodicity <NoPeriodicity>,mortars<NoMortar>>;
#endif

} // Feel

#endif /* FEELPP_PCH_H */
