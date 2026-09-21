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
   \file pdhm.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#ifndef FEELPP_PDHM_H
#define FEELPP_PDHM_H 1

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
struct Pdhmg
{
    typedef FunctionSpace<MeshType,
                          bases<Lagrange<Order,Pset,Discontinuous,Pts,Tag>>,
                          T,
                          Periodicity <NoPeriodicity>,
                          mortars<NoMortar>> type;
    typedef std::shared_ptr<type> ptrtype;
};

template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pdhm = Pdhmg<MeshType,Order,Tensor2,T,Pts,Tag>;

template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pdhms = Pdhmg<MeshType,Order,Tensor2Symm,T,Pts,Tag>;

} // meta

template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pdhm_type = typename meta::Pdhm<MeshType,Order,T,Pts,Tag>::type;
template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pdhm_ptrtype = typename meta::Pdhm<MeshType,Order,T,Pts,Tag>::ptrtype;

template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pdhms_type = typename meta::Pdhms<MeshType,Order,T,Pts,Tag>::type;
template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pdhms_ptrtype = typename meta::Pdhms<MeshType,Order,T,Pts,Tag>::ptrtype;

/**
 * \fn Pdhm<k,MeshType>
 *
 * build a function space of discontinuous matrix fields which are piecewise polynomial
 * of degree (total or in each variable) less than k.
 */
template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,
         int Tag = 0>
inline
Pdhm_ptrtype<MeshType,Order,T,Pts,Tag>
Pdhm( std::shared_ptr<MeshType> mesh, DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    using space_type = Pdhm_type<MeshType,Order,T,Pts,Tag>;
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
Pdhm( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto mesh = args.get( _mesh );
    auto dte = args.get_else( _extended_doftable, DofTableExtendedType::DEFAULT );
    auto policy = args.get_else( _fspace_reuse_policy,
                                 FunctionSpaceReusePolicy::automatic );
    using mesh_type = typename std::decay_t<decltype( mesh )>::element_type;
    using space_type = Pdhm_type<mesh_type,Order,T,Pts,Tag>;
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
 * build a function space of discontinuous symmetric matrix fields which are piecewise polynomial
 * of degree (total or in each variable) less than k.
 */
template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,
         int Tag = 0>
inline
Pdhms_ptrtype<MeshType,Order,T,Pts,Tag>
Pdhms( std::shared_ptr<MeshType> const& mesh, DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    using space_type = Pdhms_type<MeshType,Order,T,Pts,Tag>;
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
Pdhms( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto mesh = args.get( _mesh );
    auto dte = args.get_else( _extended_doftable, DofTableExtendedType::DEFAULT );
    auto policy = args.get_else( _fspace_reuse_policy,
                                 FunctionSpaceReusePolicy::automatic );
    using mesh_type = typename std::decay_t<decltype( mesh )>::element_type;
    using space_type = Pdhms_type<mesh_type,Order,T,Pts,Tag>;
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

#endif /* FEELPP_PDHM_H */
