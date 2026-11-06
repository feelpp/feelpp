/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2013-2016 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <feel/feelpoly/fecapabilities.hpp>
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
 * // Static order P2 (RuntimeOrder must agree with the static type)
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
    using space_type = Pch_type<MeshType, Order, T, Pts, Tag>;
    order = validatedOrder<typename space_type::reference_element_type>( order );
    fe::requireCapability<fe::Family::Lagrange>( *mesh, order );
    return getOrCreateWholeMeshFunctionSpace<space_type>(
        mesh, dte, FunctionSpaceReusePolicy::automatic,
        [&]()
        {
            return space_type::New( _mesh=mesh,
                                   _worldscomm=makeWorldsComm( 1, mesh->worldComm() ),
                                   _extended_doftable=dte,
                                   _runtime_order=order );
        }, Order == Dynamic ? static_cast<int>( order.value ) : -1 );
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
    fe::requireCapability<fe::Family::Lagrange>(
        *mesh, RuntimeOrder{ static_cast<uint16_type>( Order ) } );
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
    fe::requireCapability<fe::Family::Lagrange>( *mesh, order );
    return Pch_type<MeshType, Order, T, Pts, Tag>::New(
        _mesh = mesh,
        _range = std::forward<RangeType>( rangeElt ),
        _worldscomm = makeWorldsComm( 1, mesh->worldComm() ),
        _extended_doftable = dte,
        _components = components,
        _runtime_order = order );
}


/** @brief Managed whole-mesh factory; dynamic orders require `_runtime_order`. */
template<int Order,
         typename T = double,
         template<class, int, class> class Pts = PointSetFekete,
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
    auto order = [&]()
    {
        if constexpr ( Order == Dynamic )
            return args.get( _runtime_order );
        else
            return args.get_else( _runtime_order, RuntimeOrder{ Order } );
    }();
    order = validatedOrder<typename space_type::reference_element_type>( order );
    fe::requireCapability<fe::Family::Lagrange>( *mesh, order );
    return getOrCreateWholeMeshFunctionSpace<space_type>(
        mesh, dte,
        policy,
        [&]()
        {
            return space_type::New( _mesh=mesh,
                                   _worldscomm=makeWorldsComm( 1, mesh->worldComm() ),
                                   _extended_doftable=dte,
                                   _runtime_order=order );
        }, Order == Dynamic ? static_cast<int>( order.value ) : -1 );
}

} // Feel

#endif /* FEELPP_PCH_H */
