/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file geometrydescriptor.hpp
 * @brief Uniform static/runtime geometric discretization descriptors.
 */
#if !defined( FEELPP_FEELPOLY_GEOMETRYDESCRIPTOR_HPP )
#define FEELPP_FEELPOLY_GEOMETRYDESCRIPTOR_HPP 1

#include <feel/feelpoly/geometrycoefficients.hpp>
#include <feel/feelpoly/order.hpp>

#include <cstdint>
#include <stdexcept>
#include <type_traits>

namespace Feel::geometry
{
/**
 * @class GeometryDescriptor
 * @brief Mesh-block geometry metadata shared by static and runtime mappings.
 *
 * A descriptor is stored once by a homogeneous mesh block. Static descriptors
 * are empty order policies, while dynamic descriptors store exactly one
 * runtime order. Entity connectivity therefore remains independent of the
 * polynomial-order representation.
 *
 * @tparam Dim topological dimension
 * @tparam OrderSpec compile-time order or `Dynamic`
 * @tparam CellKind reference-cell topology family
 * @tparam Representation coefficient representation
 */
template <int Dim, int OrderSpec, GeometryCellKind CellKind,
          GeometryRepresentationKind Representation =
              GeometryRepresentationKind::NodalLagrange>
class GeometryDescriptor : private OrderBase<OrderSpec>
{
    using order_base_type = OrderBase<OrderSpec>;

    static_assert( Dim >= 0, "Geometry dimension must be non-negative" );
    static_assert( OrderSpec == Dynamic || OrderSpec >= 1,
                   "Static geometric order must be at least one" );
    using cardinality_type = std::conditional_t<
        CellKind == GeometryCellKind::Simplex,
        ::Feel::detail::SimplexCardinality,
        ::Feel::detail::HypercubeCardinality>;
    static_assert(
        OrderSpec == Dynamic ||
            cardinality_type::total( Dim, OrderSpec ) <=
                ::Feel::detail::maxLocalIndexCardinality,
        "Static geometry cardinality exceeds the uint16_type local-index limit" );

    /** Validate the strictly positive geometry-order invariant. */
    [[nodiscard]] static RuntimeOrder
    validatedOrder( RuntimeOrder order )
    {
        if ( order.value < 1 )
            throw std::invalid_argument(
                "Runtime geometric order must be at least one" );
        order = ::Feel::detail::validatedNodalOrder<cardinality_type>(
            Dim, order );
        (void)::Feel::detail::checkedLocalIndexValue(
            static_cast<size_type>( order.value - 1 ) * Dim );
        return order;
    }

  public:
    /** Entity-local coefficient index; guarded before every narrowing. */
    using local_coefficient_index_type = uint16_type;
    /** Largest supported entity-local coefficient cardinality. */
    static constexpr size_type max_local_coefficient_count =
        ::Feel::detail::maxLocalIndexCardinality;
    /** Authoritative compile-time order specification. */
    static constexpr int order_spec = order_base_type::order_spec;
    /** Compile-time order or the dynamic sentinel. */
    static constexpr int nOrder = order_spec;
    /** Topological dimension of the described block. */
    static constexpr int nDim = Dim;
    /** Reference-cell topology family. */
    static constexpr GeometryCellKind cellKind = CellKind;
    /** Default geometry coefficient representation. */
    static constexpr GeometryRepresentationKind default_representation =
        Representation;
    /** Whether the order is supplied at runtime. */
    static constexpr bool is_order_dynamic = order_base_type::is_order_dynamic;
    /** Whether the order is supplied at compile time. */
    static constexpr bool is_order_static = order_base_type::is_order_static;

    /** @brief Construct a compile-time geometry descriptor. */
    explicit constexpr GeometryDescriptor(
        GeometryRepresentationKind representation = Representation )
        requires( is_order_static )
        : M_representation( validatedRepresentation( representation ) )
    {}

    /**
     * @brief Construct a runtime geometry descriptor.
     * @param order strictly positive geometric polynomial order
     */
    explicit GeometryDescriptor(
        RuntimeOrder order,
        GeometryRepresentationKind representation = Representation )
        requires( is_order_dynamic )
        : order_base_type( validatedOrder( order ) ),
          M_representation( validatedRepresentation( representation ) )
    {}

    /** @brief Return the effective geometric polynomial order. */
    [[nodiscard]] constexpr uint16_type order() const noexcept
    {
        return order_base_type::order();
    }

    /** @brief Return the number of geometry coefficients on a cell. */
    [[nodiscard]] constexpr size_type coefficientCount() const noexcept
    {
        if constexpr ( CellKind == GeometryCellKind::Simplex )
            return ::Feel::detail::simplexTotal( Dim, this->order() );
        else
            return ::Feel::detail::hypercubeTotal( Dim, this->order() );
    }

    /**
     * @brief Return a geometry-aware measure quadrature order.
     *
     * The mapping Jacobian has degree `order()-1` in each reference
     * direction. This preserves the historical Feel++ measure requirement and
     * computes it from the runtime descriptor when necessary.
     */
    [[nodiscard]] constexpr uint16_type
    measureQuadratureOrder() const
    {
        return ::Feel::detail::checkedLocalIndexValue(
            static_cast<size_type>( this->order() - 1 ) * Dim );
    }

    /** @brief Return the mathematical meaning of stored coefficients. */
    [[nodiscard]] GeometryRepresentationKind representation() const noexcept
    {
        return M_representation;
    }

    /** @brief Return the mathematical meaning of stored coefficients. */
    [[nodiscard]] constexpr CoefficientSemantics
    coefficientSemantics() const noexcept
    {
        if ( M_representation ==
             GeometryRepresentationKind::NodalLagrange )
            return CoefficientSemantics::InterpolationPoint;
        return CoefficientSemantics::ControlPoint;
    }

  private:
    [[nodiscard]] static constexpr GeometryRepresentationKind
    validatedRepresentation( GeometryRepresentationKind representation )
    {
        if ( !isValidGeometryRepresentation( representation ) )
            throw std::invalid_argument(
                "Unsupported geometry coefficient representation" );
        return representation;
    }

    GeometryRepresentationKind M_representation = Representation;
};

/**
 * @brief Select the geometry cell kind associated with a Feel++ shape.
 * @tparam Shape simplex or hypercube shape type
 */
template <typename Shape>
inline constexpr GeometryCellKind geometryCellKind =
    Shape::is_simplex ? GeometryCellKind::Simplex
                      : GeometryCellKind::Hypercube;

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_GEOMETRYDESCRIPTOR_HPP
