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
#include <feel/feelpoly/order.hpp>

namespace Feel {

namespace meta {

/**
 * @brief Metafunction to create Pch function space type
 *
 * For static order (Order >= 0), creates a standard FunctionSpace.
 * For dynamic order (Order == Dynamic), creates a FunctionSpace with
 * Lagrange<Dynamic> basis that stores order at runtime.
 *
 * @tparam MeshType Mesh type
 * @tparam Order Polynomial order (can be Dynamic = -1 for runtime order)
 * @tparam T Value type
 * @tparam Pts Point set type
 * @tparam Tag Space tag
 */
template<typename MeshType,
         int Order = Dynamic,
         typename T = double,
         template<class, int, class> class Pts = PointSetFekete,
         int Tag = 0>
struct Pch
{
    using type = boost::mp11::mp_if_c<Tag==0 && std::is_same_v<T,double>,
                                      FunctionSpace<MeshType,bases<Lagrange<Order,Scalar,Continuous,Pts>>>,
                                      FunctionSpace<MeshType,bases<Lagrange<Order,Scalar,Continuous,Pts,Tag>>,T> >;
    typedef std::shared_ptr<type> ptrtype;

    //! @brief True if order is determined at runtime
    static constexpr bool is_dynamic = (Order == Dynamic);
};

} // meta

template<typename MeshType,
         int Order = Dynamic,
         typename T = double,
         template<class, int, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pch_type = typename meta::Pch<MeshType,Order,T,Pts,Tag>::type;
template<typename MeshType,
         int Order = Dynamic,
         typename T = double,
         template<class, int, class> class Pts = PointSetFekete,
         int Tag = 0>
using Pch_ptrtype = typename meta::Pch<MeshType,Order,T,Pts,Tag>::ptrtype;

template<typename MeshType,int Order = Dynamic,typename T = double, template<class, int, class> class Pts = PointSetFekete, int Tag = 0>
using Pch_element_t=typename Pch_type<MeshType,Order, T,Pts, Tag>::element_type;

template<typename MeshType,int Order = Dynamic,typename T = double,template<class, int, class> class Pts = PointSetFekete, int Tag = 0>
using Pch_element_type=Pch_element_t<MeshType,Order,T,Pts, Tag>;


/**
 * @brief Create Pch function space with explicit RuntimeOrder
 *
 * Unified interface for both static and dynamic order function spaces.
 * For static orders, the RuntimeOrder is passed through but the compile-time
 * order takes precedence. For dynamic orders (Order == Dynamic), the
 * RuntimeOrder specifies the polynomial degree.
 *
 * @code
 * // Static order P2 (RuntimeOrder ignored)
 * auto Vh = Pch<2>(mesh, RuntimeOrder(2));
 *
 * // Dynamic order P2
 * auto Vh = Pch<Dynamic>(mesh, RuntimeOrder(2));
 *
 * // Both produce equivalent function spaces
 * @endcode
 *
 * @param mesh The mesh
 * @param order Runtime order specification
 * @param dte Extended doftable type
 * @return Shared pointer to the function space
 */
template<int Order = Dynamic,
         typename T = double,
         template<class, int, class> class Pts = PointSetFekete,
         typename MeshType,
         int Tag = 0>
inline
Pch_ptrtype<MeshType, Order, T, Pts, Tag>
Pch( std::shared_ptr<MeshType> const& mesh,
     RuntimeOrder order,
     DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Pch_type<MeshType, Order, T, Pts, Tag>::New(
        _mesh = mesh,
        _worldscomm = makeWorldsComm( 1, mesh->worldComm() ),
        _extended_doftable = dte,
        _runtime_order = order );
}

/**
 * @brief Create Pch function space (static order, backward compatible)
 *
 * Build a function space of continuous functions which are piecewise polynomial
 * of degree (total or in each variable) less than k.
 *
 * @note This overload is only available for static orders (Order >= 0).
 *       For dynamic orders, use Pch<Dynamic>(mesh, RuntimeOrder(k)).
 *
 * @param mesh The mesh
 * @param dte Extended doftable type
 * @return Shared pointer to the function space
 */
template<int Order,
         typename T = double,
         template<class, int, class> class Pts = PointSetFekete,
         typename MeshType,
         int Tag = 0>
    requires ( Order >= 0 )
inline
Pch_ptrtype<MeshType, Order, T, Pts, Tag>
Pch( std::shared_ptr<MeshType> const& mesh,
     DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Pch<Order, T, Pts, MeshType, Tag>( mesh, RuntimeOrder{ static_cast<uint16_type>( Order ) }, dte );
}

/**
 * @brief Create Pch function space with range (static order, backward compatible)
 *
 * Build a function space of continuous functions which are piecewise polynomial
 * of degree (total or in each variable) less than k, restricted to a range of elements.
 *
 * @note This overload is only available for static orders (Order >= 0).
 *
 * @param mesh The mesh
 * @param rangeElt Range of elements
 * @param dte Extended doftable type
 * @param components Mesh components
 * @return Shared pointer to the function space
 */
template<int Order = Dynamic,
         typename T = double,
         template<class, int, class> class Pts = PointSetFekete,
         typename MeshType, typename RangeType,
         int Tag = 0>
    requires ( Order >= 0 )
inline
Pch_ptrtype<MeshType, Order, T, Pts, Tag>
Pch( std::shared_ptr<MeshType> const& mesh,
     RangeType&& rangeElt,
     DofTableExtendedType dte = DofTableExtendedType::DEFAULT,
     size_type components = 0 )
{
    return Pch_type<MeshType, Order, T, Pts, Tag>::New(
        _mesh = mesh,
        _range = std::forward<RangeType>( rangeElt ),
        _worldscomm = makeWorldsComm( 1, mesh->worldComm() ),
        _extended_doftable = dte,
        _components = components,
        _runtime_order = RuntimeOrder{ static_cast<uint16_type>( Order ) } );
}

/**
 * @brief Create Pch function space with range and RuntimeOrder
 *
 * Unified interface for both static and dynamic order function spaces,
 * restricted to a range of elements.
 *
 * @param mesh The mesh
 * @param rangeElt Range of elements
 * @param order Runtime order specification
 * @param dte Extended doftable type
 * @param components Mesh components
 * @return Shared pointer to the function space
 */
template<int Order,
         typename T = double,
         template<class, int, class> class Pts = PointSetFekete,
         typename MeshType, typename RangeType,
         int Tag = 0>
inline
Pch_ptrtype<MeshType, Order, T, Pts, Tag>
Pch( std::shared_ptr<MeshType> const& mesh,
     RangeType&& rangeElt,
     RuntimeOrder order,
     DofTableExtendedType dte = DofTableExtendedType::DEFAULT,
     size_type components = 0 )
{
    return Pch_type<MeshType, Order, T, Pts, Tag>::New(
        _mesh = mesh,
        _range = std::forward<RangeType>( rangeElt ),
        _worldscomm = makeWorldsComm( 1, mesh->worldComm() ),
        _extended_doftable = dte,
        _components = components,
        _runtime_order = order );
}

} // Feel

#endif /* FEELPP_PCH_H */
