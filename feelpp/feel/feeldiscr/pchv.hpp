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
   \file pchv.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_PCHV_HPP)
#define FEELPP_PCHV_HPP 1

#include <boost/mp11/utility.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/functionspacemanager.hpp>

namespace Feel
{

namespace meta
{

template<typename MeshType,
         int Order,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename T = double,
         int Tag = 0>
struct Pchv
{
    using type = FunctionSpace<MeshType, bases<Lagrange<Order, Vectorial, Continuous, Pts, Tag>>,T>;
    typedef std::shared_ptr<type> ptrtype;
};

} // meta

template<typename MeshType,
         int Order,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename T = double,
         int Tag = 0>
using Pchv_type = typename meta::Pchv<MeshType,Order,Pts,T,Tag>::type;
template<typename MeshType,
         int Order,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename T = double,
         int Tag = 0>
using Pchv_ptrtype = typename meta::Pchv<MeshType,Order,Pts,T,Tag>::ptrtype;

template<typename MeshType,int Order,template<class, uint16_type, class> class Pts = PointSetFekete,typename T = double,int Tag = 0>
using Pchv_element_t=typename Pchv_type<MeshType,Order,Pts,T,Tag>::element_type;

template<typename MeshType,int Order,template<class, uint16_type, class> class Pts = PointSetFekete,typename T = double,int Tag = 0>
using Pchv_element_type=Pchv_element_t<MeshType,Order,Pts,T,Tag>;


/**
 * @brief Build a whole-mesh continuous vector Lagrange function space.
 *
 * The request uses @ref FunctionSpaceReusePolicy::automatic. It therefore
 * reuses a managed space when global reuse is enabled and otherwise preserves
 * the historical always-new behavior.
 *
 * @tparam Order polynomial order
 * @tparam Pts interpolation point-set family
 * @tparam MeshType concrete mesh type
 * @tparam T coefficient value type
 * @tparam Tag basis tag used to distinguish otherwise identical spaces
 * @param mesh mesh on which the function space is defined
 * @param dte extended DOF-table mode
 * @return continuous vector function space
 */
template<int Order,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,
         typename T = double,
         int Tag = 0>
inline
Pchv_ptrtype<MeshType,Order,Pts,T,Tag>
Pchv( std::shared_ptr<MeshType> const& mesh, DofTableExtendedType dte = DofTableExtendedType::DEFAULT  )
{
    using space_type = Pchv_type<MeshType,Order,Pts,T,Tag>;
    return getOrCreateWholeMeshFunctionSpace<space_type>(
        mesh, dte,
        FunctionSpaceReusePolicy::automatic,
        [&]()
        {
            return space_type::New( _mesh=mesh,
                                    _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                    _extended_doftable=dte );
        } );
}

/**
 * @brief Build a whole-mesh continuous vector Lagrange space using named arguments.
 *
 * Supported arguments are the required @c _mesh and the optional
 * @c _extended_doftable and @c _fspace_reuse_policy keywords. The reuse policy
 * defaults to @ref FunctionSpaceReusePolicy::automatic.
 *
 * @tparam Order polynomial order
 * @tparam Pts interpolation point-set family
 * @tparam T coefficient value type
 * @tparam Tag basis tag used to distinguish otherwise identical spaces
 * @tparam Ts named-argument types
 * @param v named arguments controlling mesh, DOF table, and reuse policy
 * @return continuous vector function space
 */
template<int Order,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename T = double,
         int Tag = 0,
         typename... Ts>
    requires ( sizeof...( Ts ) != 0 ) && ( NA::is_named_argument_v<Ts> && ... )
inline auto
Pchv( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto mesh = args.get( _mesh );
    auto dte = args.get_else( _extended_doftable, DofTableExtendedType::DEFAULT );
    auto policy = args.get_else( _fspace_reuse_policy,
                                 FunctionSpaceReusePolicy::automatic );
    using mesh_type = typename std::decay_t<decltype( mesh )>::element_type;
    using space_type = Pchv_type<mesh_type,Order,Pts,T,Tag>;
    return getOrCreateWholeMeshFunctionSpace<space_type>(
        mesh, dte,
        policy,
        [&]()
        {
            return space_type::New( _mesh=mesh,
                                    _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                    _extended_doftable=dte );
        } );
}

/**
 Given a \p mesh, build a function space of vectorial continuous function
 which are piecewise polynomial of degree (total or in each variable) less
 than k using Lagrange basis functions
 */
template<int Order,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,typename RangeType,
         typename T=double,
         int Tag = 0>
inline
Pchv_ptrtype<MeshType,Order,Pts,T,Tag>
Pchv( std::shared_ptr<MeshType> const& mesh, RangeType && rangeElt, DofTableExtendedType dte = DofTableExtendedType::DEFAULT  )
{
    return Pchv_type<MeshType,Order,Pts,T,Tag>::New( _mesh=mesh,
                                                   _range=std::forward<RangeType>(rangeElt),
                                                   _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                                   _extended_doftable=dte );
}


}

#endif /* FEELPP_PCHV_HPP */
