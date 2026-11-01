/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
class DofTopology : public OrderBase<OrderSpec>
{
    static_assert( Dim >= 0, "DofTopology dimension must be >= 0" );
    static_assert(
        std::is_same_v<CellShape, SimplexTag> ||
            std::is_same_v<CellShape, HypercubeTag>,
        "DofTopology supports only SimplexTag and HypercubeTag" );

  public:
    using order_base_type = OrderBase<OrderSpec>;
    using cardinality_type = std::conditional_t<
        std::is_same_v<CellShape, SimplexTag>,
        detail::SimplexCardinality,
        detail::HypercubeCardinality>;
    static_assert(
        OrderSpec == Dynamic ||
            cardinality_type::total( Dim, OrderSpec ) <=
                detail::maxLocalIndexCardinality,
        "Static nodal DoF cardinality exceeds the uint16_type local-index limit" );
    using order_base_type::order;

    static constexpr int nDim = Dim;
    static constexpr int order_spec = order_base_type::order_spec;
    static constexpr int nOrder = order_spec;
    static constexpr bool is_order_dynamic = order_base_type::is_order_dynamic;
    static constexpr bool is_order_static = order_base_type::is_order_static;

    constexpr DofTopology() noexcept = default;
    constexpr explicit DofTopology( RuntimeOrder runtimeOrder )
        requires (OrderSpec == Dynamic)
        : order_base_type(
              detail::validatedNodalOrder<cardinality_type>(
                  Dim, runtimeOrder ) )
    {}
    constexpr explicit DofTopology( uint16_type runtimeOrder )
        requires (OrderSpec == Dynamic)
        : DofTopology( RuntimeOrder{ runtimeOrder } )
    {}
    explicit DofTopology( int runtimeOrder )
        requires (OrderSpec == Dynamic)
        : DofTopology( RuntimeOrder::checked( runtimeOrder ) )
    {}

    [[nodiscard]] constexpr uint16_type localDof() const
    {
        return detail::checkedLocalIndexValue(
            cardinality_type::total( Dim, this->order() ) );
    }
    [[nodiscard]] constexpr uint16_type dofPerVertex() const noexcept
    {
        return cardinality_type::perVertex( Dim, this->order() );
    }
    [[nodiscard]] constexpr uint16_type dofPerEdge() const noexcept
    {
        return cardinality_type::perEdge( Dim, this->order() );
    }
    [[nodiscard]] constexpr uint16_type dofPerFace() const noexcept
    {
        return cardinality_type::perFace( Dim, this->order() );
    }
    [[nodiscard]] constexpr uint16_type dofPerVolume() const noexcept
    {
        return cardinality_type::perVolume( Dim, this->order() );
    }

    [[nodiscard]] constexpr uint16_type dofPerEntity( uint16_type entityDim,
                                                       uint16_type entityIndex = 0 ) const noexcept
    {
        (void)entityIndex;
        return cardinality_type::onEntity( Dim, entityDim, this->order() );
    }

};

/**
 * @brief Convenience alias when the cell type is known instead of the shape tag.
 */
template <typename CellShape, int Dim, int OrderSpec>
using DofTopologyFor = DofTopology<cellShapeTag_t<CellShape>, Dim, OrderSpec>;

} // namespace Feel

#endif /* FEELPP_FEELPOLY_DOFTOPOLOGY_HPP */
