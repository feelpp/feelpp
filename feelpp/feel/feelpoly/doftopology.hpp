/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2026-02-16

 Copyright (C) 2026 Feel++ Consortium

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
 * @file doftopology.hpp
 * @brief Uniform DOF topology for simplex/hypercube cells.
 */
#ifndef FEELPP_FEELPOLY_DOFTOPOLOGY_HPP
#define FEELPP_FEELPOLY_DOFTOPOLOGY_HPP 1

#include <type_traits>

#include <feel/feelpoly/order.hpp>
#include <feel/feelmesh/traits.hpp>

namespace Feel
{

struct SimplexTag {};
struct HypercubeTag {};

namespace detail
{
template <typename... T>
inline constexpr bool alwaysFalseV = false;
} // namespace detail

/**
 * @brief Map a cell type to a DOF topology tag.
 */
template <typename CellShape>
struct CellShapeTag
{
  private:
    using cell_shape_type = std::remove_cvref_t<CellShape>;

  public:
    using type = std::conditional_t<
        is_simplex_v<cell_shape_type>,
        SimplexTag,
        std::conditional_t<
            is_hypercube_v<cell_shape_type>,
            HypercubeTag,
            void>>;
};

template <typename CellShape>
using cellShapeTag_t = typename CellShapeTag<CellShape>::type;

/**
 * @brief DOF topology by shape tag, topological dimension and order.
 *
 * Only simplex/hypercube are supported in this phase.
 */
template <typename CellShape, int Dim, int OrderSpec>
class DofTopology
{
    static_assert( detail::alwaysFalseV<CellShape>,
                   "DofTopology: unsupported CellShape specialization" );
};

template <int Dim, int OrderSpec>
class DofTopology<SimplexTag, Dim, OrderSpec> : public OrderBase<OrderSpec>
{
    static_assert( Dim >= 0, "DofTopology simplex dimension must be >= 0" );

  public:
    using order_base_type = OrderBase<OrderSpec>;
    using order_base_type::order;
    using order_base_type::runtimeOrder;

    static constexpr int nDim = Dim;
    static constexpr int nOrder = OrderSpec;
    static constexpr bool is_order_dynamic = order_base_type::is_order_dynamic;
    static constexpr bool is_order_static = order_base_type::is_order_static;

    constexpr DofTopology() noexcept = default;
    constexpr explicit DofTopology( RuntimeOrder runtimeOrder ) noexcept
        requires (OrderSpec == Dynamic)
        : order_base_type( runtimeOrder )
    {}
    constexpr explicit DofTopology( uint16_type runtimeOrder ) noexcept
        requires (OrderSpec == Dynamic)
        : order_base_type( runtimeOrder )
    {}
    explicit DofTopology( int runtimeOrder )
        requires (OrderSpec == Dynamic)
        : order_base_type( RuntimeOrder::checked( runtimeOrder ) )
    {}

    [[nodiscard]] constexpr uint16_type localDof() const noexcept
    {
        if constexpr ( OrderSpec == Dynamic )
            return static_cast<uint16_type>( ::Feel::detail::simplexTotal(
                static_cast<uint16_type>( Dim ), this->order() ) );
        else
            return static_cast<uint16_type>( ::Feel::detail::simplexDofStatic<Dim, OrderSpec>() );
    }
    [[nodiscard]] constexpr uint16_type dofPerVertex() const noexcept
    {
        if constexpr ( OrderSpec == Dynamic )
            return ::Feel::detail::simplexPerVertex( static_cast<uint16_type>( Dim ), this->order() );
        else
            return ::Feel::detail::simplexPerVertex( static_cast<uint16_type>( Dim ),
                                                     static_cast<uint16_type>( OrderSpec ) );
    }
    [[nodiscard]] constexpr uint16_type dofPerEdge() const noexcept
    {
        if constexpr ( OrderSpec == Dynamic )
            return ::Feel::detail::simplexPerEdge( static_cast<uint16_type>( Dim ), this->order() );
        else
            return ::Feel::detail::simplexPerEdge( static_cast<uint16_type>( Dim ),
                                                   static_cast<uint16_type>( OrderSpec ) );
    }
    [[nodiscard]] constexpr uint16_type dofPerFace() const noexcept
    {
        if constexpr ( OrderSpec == Dynamic )
            return ::Feel::detail::simplexPerFace( static_cast<uint16_type>( Dim ), this->order() );
        else
            return ::Feel::detail::simplexPerFace( static_cast<uint16_type>( Dim ),
                                                   static_cast<uint16_type>( OrderSpec ) );
    }
    [[nodiscard]] constexpr uint16_type dofPerVolume() const noexcept
    {
        if constexpr ( OrderSpec == Dynamic )
            return ::Feel::detail::simplexPerVolume( static_cast<uint16_type>( Dim ), this->order() );
        else
            return ::Feel::detail::simplexPerVolume( static_cast<uint16_type>( Dim ),
                                                     static_cast<uint16_type>( OrderSpec ) );
    }

    [[nodiscard]] constexpr uint16_type dofPerEntity( uint16_type entityDim,
                                                       uint16_type entityIndex = 0 ) const noexcept
    {
        (void)entityIndex;
        switch ( entityDim )
        {
        case 0:
            return dofPerVertex();
        case 1:
            return dofPerEdge();
        case 2:
            return dofPerFace();
        case 3:
            return dofPerVolume();
        default:
            return 0;
        }
    }

    [[nodiscard]] constexpr uint16_type runtimeLocalDof() const noexcept { return localDof(); }
    [[nodiscard]] constexpr uint16_type runtimeDofPerVertex() const noexcept { return dofPerVertex(); }
    [[nodiscard]] constexpr uint16_type runtimeDofPerEdge() const noexcept { return dofPerEdge(); }
    [[nodiscard]] constexpr uint16_type runtimeDofPerFace() const noexcept { return dofPerFace(); }
    [[nodiscard]] constexpr uint16_type runtimeDofPerVolume() const noexcept { return dofPerVolume(); }
};

template <int Dim, int OrderSpec>
class DofTopology<HypercubeTag, Dim, OrderSpec> : public OrderBase<OrderSpec>
{
    static_assert( Dim >= 0, "DofTopology hypercube dimension must be >= 0" );

  public:
    using order_base_type = OrderBase<OrderSpec>;
    using order_base_type::order;
    using order_base_type::runtimeOrder;

    static constexpr int nDim = Dim;
    static constexpr int nOrder = OrderSpec;
    static constexpr bool is_order_dynamic = order_base_type::is_order_dynamic;
    static constexpr bool is_order_static = order_base_type::is_order_static;

    constexpr DofTopology() noexcept = default;
    constexpr explicit DofTopology( RuntimeOrder runtimeOrder ) noexcept
        requires (OrderSpec == Dynamic)
        : order_base_type( runtimeOrder )
    {}
    constexpr explicit DofTopology( uint16_type runtimeOrder ) noexcept
        requires (OrderSpec == Dynamic)
        : order_base_type( runtimeOrder )
    {}
    explicit DofTopology( int runtimeOrder )
        requires (OrderSpec == Dynamic)
        : order_base_type( RuntimeOrder::checked( runtimeOrder ) )
    {}

    [[nodiscard]] constexpr uint16_type localDof() const noexcept
    {
        if constexpr ( OrderSpec == Dynamic )
            return static_cast<uint16_type>( ::Feel::detail::hypercubeTotal(
                static_cast<uint16_type>( Dim ), this->order() ) );
        else
            return static_cast<uint16_type>( ::Feel::detail::hypercubeDofStatic<Dim, OrderSpec>() );
    }
    [[nodiscard]] constexpr uint16_type dofPerVertex() const noexcept
    {
        if constexpr ( OrderSpec == Dynamic )
            return ::Feel::detail::hypercubePerVertex( static_cast<uint16_type>( Dim ), this->order() );
        else
            return ::Feel::detail::hypercubePerVertex( static_cast<uint16_type>( Dim ),
                                                       static_cast<uint16_type>( OrderSpec ) );
    }
    [[nodiscard]] constexpr uint16_type dofPerEdge() const noexcept
    {
        if constexpr ( OrderSpec == Dynamic )
            return ::Feel::detail::hypercubePerEdge( static_cast<uint16_type>( Dim ), this->order() );
        else
            return ::Feel::detail::hypercubePerEdge( static_cast<uint16_type>( Dim ),
                                                     static_cast<uint16_type>( OrderSpec ) );
    }
    [[nodiscard]] constexpr uint16_type dofPerFace() const noexcept
    {
        if constexpr ( OrderSpec == Dynamic )
            return ::Feel::detail::hypercubePerFace( static_cast<uint16_type>( Dim ), this->order() );
        else
            return ::Feel::detail::hypercubePerFace( static_cast<uint16_type>( Dim ),
                                                     static_cast<uint16_type>( OrderSpec ) );
    }
    [[nodiscard]] constexpr uint16_type dofPerVolume() const noexcept
    {
        if constexpr ( OrderSpec == Dynamic )
            return ::Feel::detail::hypercubePerVolume( static_cast<uint16_type>( Dim ), this->order() );
        else
            return ::Feel::detail::hypercubePerVolume( static_cast<uint16_type>( Dim ),
                                                       static_cast<uint16_type>( OrderSpec ) );
    }

    [[nodiscard]] constexpr uint16_type dofPerEntity( uint16_type entityDim,
                                                       uint16_type entityIndex = 0 ) const noexcept
    {
        (void)entityIndex;
        switch ( entityDim )
        {
        case 0:
            return dofPerVertex();
        case 1:
            return dofPerEdge();
        case 2:
            return dofPerFace();
        case 3:
            return dofPerVolume();
        default:
            return 0;
        }
    }

    [[nodiscard]] constexpr uint16_type runtimeLocalDof() const noexcept { return localDof(); }
    [[nodiscard]] constexpr uint16_type runtimeDofPerVertex() const noexcept { return dofPerVertex(); }
    [[nodiscard]] constexpr uint16_type runtimeDofPerEdge() const noexcept { return dofPerEdge(); }
    [[nodiscard]] constexpr uint16_type runtimeDofPerFace() const noexcept { return dofPerFace(); }
    [[nodiscard]] constexpr uint16_type runtimeDofPerVolume() const noexcept { return dofPerVolume(); }
};

/**
 * @brief Convenience alias when the cell type is known instead of the shape tag.
 */
template <typename CellShape, int Dim, int OrderSpec>
using DofTopologyFor = DofTopology<cellShapeTag_t<CellShape>, Dim, OrderSpec>;

} // namespace Feel

#endif /* FEELPP_FEELPOLY_DOFTOPOLOGY_HPP */
