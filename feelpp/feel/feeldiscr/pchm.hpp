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
   \file pchm.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#ifndef FEELPP_PCHM_H
#define FEELPP_PCHM_H 1

#include <boost/mp11/utility.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/functionspacemanager.hpp>

namespace Feel {

namespace meta {
template<typename MeshType,
         int Order,
         template <uint16_type> class Pset = Tensor2,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
struct Pchmg
{
    using type = boost::mp11::mp_if_c<Tag == 0 && std::is_same_v<T, double>,
                                      FunctionSpace<MeshType, bases<Lagrange<Order, Pset, Continuous, Pts>>>,
                                      FunctionSpace<MeshType, bases<Lagrange<Order, Pset, Continuous, Pts, Tag>>, T>>;
    typedef std::shared_ptr<type> ptrtype;
};

template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pchm = Pchmg<MeshType,Order,Tensor2,T,Pts,Tag>;

template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pchms = Pchmg<MeshType,Order,Tensor2Symm,T,Pts,Tag>;

} // meta

template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pchm_type = typename meta::Pchm<MeshType,Order,T,Pts,Tag>::type;
template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pchm_ptrtype = typename meta::Pchm<MeshType,Order,T,Pts,Tag>::ptrtype;

template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pchms_type = typename meta::Pchms<MeshType,Order,T,Pts,Tag>::type;
template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pchms_ptrtype = typename meta::Pchms<MeshType,Order,T,Pts,Tag>::ptrtype;

/**
 * \fn Pchm<k,MeshType>
 *
 * build a function space of continuous matrix fields which are piecewise polynomial
 * of degree (total or in each variable) less than k.
 */
template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,
         int Tag = 0>
inline
Pchm_ptrtype<MeshType,Order,T,Pts,Tag>
Pchm( std::shared_ptr<MeshType> mesh, DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    using space_type = Pchm_type<MeshType,Order,T,Pts,Tag>;
    return getOrCreateWholeMeshFunctionSpace<space_type>(
        mesh, dte,
        FunctionSpaceReusePolicy::automatic,
        [&]()
        {
            return space_type::New( _mesh=mesh,
                                    _worldscomm=makeWorldsComm( 1,mesh->worldCommPtr() ),
                                    _extended_doftable=dte );
        } );
}

template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0,
         typename... Ts>
    requires ( sizeof...( Ts ) != 0 ) && ( NA::is_named_argument_v<Ts> && ... )
inline auto
Pchm( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto mesh = args.get( _mesh );
    auto dte = args.get_else( _extended_doftable, DofTableExtendedType::DEFAULT );
    auto policy = args.get_else( _fspace_reuse_policy,
                                 FunctionSpaceReusePolicy::automatic );
    using mesh_type = typename std::decay_t<decltype( mesh )>::element_type;
    using space_type = Pchm_type<mesh_type,Order,T,Pts,Tag>;
    return getOrCreateWholeMeshFunctionSpace<space_type>(
        mesh, dte,
        policy,
        [&]()
        {
            return space_type::New( _mesh=mesh,
                                    _worldscomm=makeWorldsComm( 1,mesh->worldCommPtr() ),
                                    _extended_doftable=dte );
        } );
}

/**
 * build a function space of continuous symmetric matrix fields which are piecewise polynomial
 * of degree (total or in each variable) less than k.
 */
template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,
         int Tag = 0>
inline
Pchms_ptrtype<MeshType,Order,T,Pts,Tag>
Pchms( std::shared_ptr<MeshType> const& mesh, DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    using space_type = Pchms_type<MeshType,Order,T,Pts,Tag>;
    return getOrCreateWholeMeshFunctionSpace<space_type>(
        mesh, dte,
        FunctionSpaceReusePolicy::automatic,
        [&]()
        {
            return space_type::New( _mesh=mesh,
                                    _worldscomm=makeWorldsComm( 1,mesh->worldCommPtr() ),
                                    _extended_doftable=dte );
        } );
}

template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0,
         typename... Ts>
    requires ( sizeof...( Ts ) != 0 ) && ( NA::is_named_argument_v<Ts> && ... )
inline auto
Pchms( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto mesh = args.get( _mesh );
    auto dte = args.get_else( _extended_doftable, DofTableExtendedType::DEFAULT );
    auto policy = args.get_else( _fspace_reuse_policy,
                                 FunctionSpaceReusePolicy::automatic );
    using mesh_type = typename std::decay_t<decltype( mesh )>::element_type;
    using space_type = Pchms_type<mesh_type,Order,T,Pts,Tag>;
    return getOrCreateWholeMeshFunctionSpace<space_type>(
        mesh, dte,
        policy,
        [&]()
        {
            return space_type::New( _mesh=mesh,
                                    _worldscomm=makeWorldsComm( 1,mesh->worldCommPtr() ),
                                    _extended_doftable=dte );
        } );
}


} // Feel

#endif /* FEELPP_PCHM_H */
