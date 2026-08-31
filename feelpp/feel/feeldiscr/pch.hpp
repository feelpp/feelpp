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
#include <boost/mp11/utility.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/functionspacemanager.hpp>

namespace Feel {

namespace meta {
template<typename MeshType,
         int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0>
struct Pch
{
    using type = boost::mp11::mp_if_c<Tag==0 && std::is_same_v<T,double>,
                                      FunctionSpace<MeshType,bases<Lagrange<Order,Scalar,Continuous,Pts>>>,
                                      FunctionSpace<MeshType,bases<Lagrange<Order,Scalar,Continuous,Pts,Tag>>,T> >;
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
 * @brief Build a whole-mesh continuous scalar Lagrange function space.
 *
 * The request uses @ref FunctionSpaceReusePolicy::automatic. It therefore
 * reuses a managed space when global reuse is enabled and otherwise preserves
 * the historical always-new behavior.
 *
 * @tparam Order polynomial order
 * @tparam T coefficient value type
 * @tparam Pts interpolation point-set family
 * @tparam MeshType concrete mesh type
 * @tparam Tag basis tag used to distinguish otherwise identical spaces
 * @param mesh mesh on which the function space is defined
 * @param dte extended DOF-table mode
 * @return continuous scalar function space
 */
template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,
         int Tag = 0>
inline
Pch_ptrtype<MeshType,Order,T,Pts,Tag>
Pch( std::shared_ptr<MeshType> const& mesh, DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    using space_type = Pch_type<MeshType,Order,T,Pts,Tag>;
    return getOrCreateFunctionSpace<space_type>(
        mesh,
        FunctionSpaceManagerOptions{ normalizeFunctionSpaceDofTable( dte ),
                                        MESH_RENUMBER | MESH_CHECK },
        FunctionSpaceReusePolicy::automatic,
        [&]()
        {
            return space_type::New( _mesh=mesh,
                                    _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                    _extended_doftable=dte );
        } );
}

/**
 * @brief Build a whole-mesh continuous scalar Lagrange space using named arguments.
 *
 * Supported arguments are the required @c _mesh and the optional
 * @c _extended_doftable and @c _fspace_reuse_policy keywords. The reuse policy
 * defaults to @ref FunctionSpaceReusePolicy::automatic.
 *
 * @tparam Order polynomial order
 * @tparam T coefficient value type
 * @tparam Pts interpolation point-set family
 * @tparam Tag basis tag used to distinguish otherwise identical spaces
 * @tparam Ts named-argument types
 * @param v named arguments controlling mesh, DOF table, and reuse policy
 * @return continuous scalar function space
 */
template<int Order,
         typename T = double,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         int Tag = 0,
         typename... Ts>
    requires ( sizeof...( Ts ) != 0 ) && ( NA::is_named_argument_v<Ts> && ... )
inline auto
Pch( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto mesh = args.get( _mesh );
    auto dte = args.get_else( _extended_doftable, DofTableExtendedType::DEFAULT );
    auto policy = args.get_else( _fspace_reuse_policy,
                                 FunctionSpaceReusePolicy::automatic );
    using mesh_type = typename std::decay_t<decltype( mesh )>::element_type;
    using space_type = Pch_type<mesh_type,Order,T,Pts,Tag>;
    return getOrCreateFunctionSpace<space_type>(
        mesh,
        FunctionSpaceManagerOptions{ normalizeFunctionSpaceDofTable( dte ),
                                        MESH_RENUMBER | MESH_CHECK },
        policy,
        [&]()
        {
            return space_type::New( _mesh=mesh,
                                    _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                    _extended_doftable=dte );
        } );
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
         typename MeshType, typename RangeType,
         int Tag = 0>
inline
Pch_ptrtype<MeshType,Order,T,Pts,Tag>
Pch( std::shared_ptr<MeshType> const& mesh, RangeType&& rangeElt, DofTableExtendedType dte = DofTableExtendedType::DEFAULT, size_type components = 0 )
{
    return Pch_type<MeshType,Order,T,Pts,Tag>::New( _mesh=mesh,
                                                    _range=std::forward<RangeType>(rangeElt),
                                                    _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                                    _extended_doftable=dte,
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

#endif

} // Feel

#endif /* FEELPP_PCH_H */
