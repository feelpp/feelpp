/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file order.hpp
 * @brief Order tags and DOF calculators for static/dynamic orders.
 * @author Christophe Prud'homme
 * @date 2026-01-02
 */
#ifndef FEELPP_FEELPOLY_ORDER_HPP
#define FEELPP_FEELPOLY_ORDER_HPP 1

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include <feel/feelcore/feeltypes.hpp>

namespace Feel
{

/**
 * @brief Eigen-style dynamic sentinel for polynomial order.
 *
 * Using -1 as the sentinel value matches Eigen's convention for dynamic sizes.
 * This allows type-level distinction between static and dynamic orders:
 * - Simplex<2, 1> : static P1 geometry
 * - Simplex<2, Dynamic> : dynamic order geometry
 */
inline constexpr int Dynamic = -1;

/**
 * @brief Order tag with static value.
 * @tparam N compile-time polynomial order, or `Dynamic` for runtime order
 */
template <int N>
struct order_t
{
    /** Compile-time order carried by the tag. */
    static constexpr int value = N;
    /** Whether this tag denotes runtime order. */
    static constexpr bool is_dynamic = (N == Dynamic);
    /** Whether this tag denotes compile-time order. */
    static constexpr bool is_static = !is_dynamic;
};

/**
 * @brief Type-safe wrapper for runtime order specification.
 *
 * Used to construct types with dynamic order at the type level:
 * @code
 * Simplex<2, Dynamic> simplex(RuntimeOrder(2));  // P2 geometry
 * Lagrange<Dynamic, Scalar>::apply<2,2>::type fe(RuntimeOrder(3));  // P3 element
 * @endcode
 */
struct RuntimeOrder
{
    /** Validated non-negative runtime order value. */
    uint16_type value;

    /**
     * @brief Construct from an unsigned order.
     * @param v non-negative polynomial order
     */
    constexpr explicit RuntimeOrder( uint16_type v ) noexcept : value( v ) {}

    /**
     * @brief Construct from a signed order with range validation.
     * @param v polynomial order in the `uint16_type` range
     * @throws std::invalid_argument if `v` is negative or exceeds the storage range
     */
    constexpr explicit RuntimeOrder( int v )
        : value( v >= 0 && v <= std::numeric_limits<uint16_type>::max()
                     ? static_cast<uint16_type>( v )
                     : throw std::invalid_argument( "RuntimeOrder is outside uint16_type range" ) )
    {}

    /**
     * @brief Construct a validated runtime-order tag from a signed integer.
     * @param v polynomial order in the `uint16_type` range
     * @return validated runtime-order tag
     * @throws std::invalid_argument if `v` is outside the supported range
     */
    [[nodiscard]] static RuntimeOrder checked( int v )
    {
        return RuntimeOrder( v );
    }

    /** @brief Return the stored order value. */
    [[nodiscard]] constexpr operator uint16_type() const noexcept { return value; }
};

/**
 * @brief Common order metadata and runtime storage for static/dynamic order types.
 *
 * Primary template handles static orders (`Order >= 0`).
 * Dynamic order (`Order == Dynamic`) has a dedicated specialization below.
 */
template <int Order>
class OrderBase
{
    static_assert( Order >= 0, "OrderBase static order must be non-negative" );

  public:
    /** Compile-time polynomial order. */
    static constexpr int nOrder = Order;
    /** Static-order policies never carry runtime order. */
    static constexpr bool is_order_dynamic = false;
    /** Static-order policies expose their order at compile time. */
    static constexpr bool is_order_static = true;

    /** @brief Construct the compile-time order policy. */
    constexpr OrderBase() noexcept = default;
    /**
     * @brief Accept a runtime tag without changing a compile-time order.
     *
     * This overload permits generic policy-based construction while retaining
     * the compile-time order as the single source of truth.
     */
    constexpr explicit OrderBase( RuntimeOrder /*unused*/ ) noexcept {}

    /** @brief Return the compile-time order as the common unsigned order type. */
    [[nodiscard]] static constexpr uint16_type staticOrder() noexcept
    {
        return static_cast<uint16_type>( Order );
    }
    /** @brief Return the effective polynomial order. */
    [[nodiscard]] constexpr uint16_type order() const noexcept { return staticOrder(); }
    /** @brief Return the effective order through the uniform runtime API. */
    [[nodiscard]] constexpr uint16_type runtimeOrder() const noexcept { return order(); }
};

/**
 * @brief Dynamic-order specialization of OrderBase.
 */
template <>
class OrderBase<Dynamic>
{
  public:
    /** Dynamic-order sentinel. */
    static constexpr int nOrder = Dynamic;
    /** Dynamic policies carry an order in each object. */
    static constexpr bool is_order_dynamic = true;
    /** Dynamic policies do not expose a concrete compile-time order. */
    static constexpr bool is_order_static = false;

    /** @brief Construct a dynamic policy with order zero. */
    constexpr OrderBase() noexcept = default;
    /** @brief Construct a dynamic policy from a type-safe runtime order. */
    constexpr explicit OrderBase( RuntimeOrder order ) noexcept : M_runtimeOrder( order ) {}
    /** @brief Construct a dynamic policy from an unsigned runtime order. */
    constexpr explicit OrderBase( uint16_type order ) noexcept : M_runtimeOrder( order ) {}
    /**
     * @brief Construct a dynamic policy from a checked signed order.
     * @throws std::invalid_argument if `order` is outside the supported range
     */
    explicit OrderBase( int order ) : M_runtimeOrder( RuntimeOrder::checked( order ) ) {}

    /** @brief Return the effective polynomial order. */
    [[nodiscard]] constexpr uint16_type order() const noexcept { return M_runtimeOrder.value; }
    /** @brief Return the effective order through the uniform runtime API. */
    [[nodiscard]] constexpr uint16_type runtimeOrder() const noexcept { return order(); }

  protected:
    /** @brief Replace the stored runtime order with a validated tag. */
    constexpr void setOrder( RuntimeOrder order ) noexcept { M_runtimeOrder = order; }
    /** @brief Replace the stored runtime order from an unsigned value. */
    constexpr void setOrder( uint16_type order ) noexcept { M_runtimeOrder = RuntimeOrder( order ); }
    /** @brief Replace the stored runtime order from a checked signed value. */
    void setOrderChecked( int order ) { M_runtimeOrder = RuntimeOrder::checked( order ); }

    /** Runtime order carried only by the dynamic specialization. */
    RuntimeOrder M_runtimeOrder{0};
};

/**
 * @brief Detect whether a type exposes an `is_order_dynamic` flag.
 */
template <typename T>
inline constexpr bool hasOrderDynamicFlag = requires {
    std::remove_cvref_t<T>::is_order_dynamic;
};

/**
 * @brief Query whether a type uses runtime polynomial order.
 *
 * Types not exposing `is_order_dynamic` are considered static by default.
 */
template <typename T>
inline constexpr bool orderIsDynamic =
    hasOrderDynamicFlag<T> && std::remove_cvref_t<T>::is_order_dynamic;

/**
 * @brief Query whether a type uses compile-time polynomial order.
 */
template <typename T>
inline constexpr bool orderIsStatic = !orderIsDynamic<T>;

namespace detail
{
constexpr size_type powInt( size_type base, int exp )
{
    size_type result = 1;
    for ( int i = 0; i < exp; ++i )
        result *= base;
    return result;
}

constexpr size_type binomial( int n, int k )
{
    if ( k < 0 || k > n )
        return 0;
    if ( k > n - k )
        k = n - k;
    size_type result = 1;
    for ( int i = 1; i <= k; ++i )
        result = result * static_cast<size_type>( n - k + i ) / static_cast<size_type>( i );
    return result;
}

template <int Dim, int Order>
constexpr size_type simplexDofStatic()
{
    static_assert( Dim >= 0, "Dimension must be non-negative" );
    static_assert( Order >= 0, "Order must be non-negative" );
    return binomial( Dim + Order, Dim );
}

template <int Dim, int Order>
constexpr size_type hypercubeDofStatic()
{
    static_assert( Dim >= 0, "Dimension must be non-negative" );
    static_assert( Order >= 0, "Order must be non-negative" );
    return powInt( static_cast<size_type>( Order + 1 ), Dim );
}

template <int Dim>
inline size_type simplexDofDynamic( int order )
{
    if ( order < 0 )
        return 0;
    return binomial( Dim + order, Dim );
}

template <int Dim>
inline size_type hypercubeDofDynamic( int order )
{
    if ( order < 0 )
        return 0;
    return powInt( static_cast<size_type>( order + 1 ), Dim );
}

//
// Runtime point calculations for dynamic-order simplices
//

/**
 * @brief Compute number of points for simplex of given dimension and order.
 */
constexpr size_type simplexTotal( uint16_type dim, uint16_type order )
{
    if ( order == 0 )
        return 1;
    switch ( dim )
    {
    case 0:
        return 1;
    case 1:
        return order + 1;
    case 2:
        return ( order + 1 ) * ( order + 2 ) / 2;
    case 3:
        return ( order + 1 ) * ( order + 2 ) * ( order + 3 ) / 6;
    default:
        return 0;
    }
}

/**
 * @brief Points per vertex for simplex.
 *
 * For 0D simplex, P0 has one point dof on the unique point/vertex.
 * For higher dimensions, P0 has no vertex dofs.
 */
constexpr uint16_type simplexPerVertex( uint16_type dim, uint16_type order )
{
    if ( dim == 0 )
        return 1;
    return ( order == 0 ) ? 0 : 1;
}

/**
 * @brief Backward-compatible overload (assumes dim > 0 semantics).
 */
constexpr uint16_type simplexPerVertex( uint16_type order )
{
    return simplexPerVertex( 1, order );
}

/**
 * @brief Points per edge (interior edge points).
 */
constexpr uint16_type simplexPerEdge( uint16_type dim, uint16_type order )
{
    if ( dim < 1 )
        return 0;
    if ( order == 0 )
        return ( dim == 1 ) ? 1 : 0;
    if ( order < 2 )
        return 0;
    return order - 1;
}

/**
 * @brief Points per face (interior face points for triangular faces).
 */
constexpr uint16_type simplexPerFace( uint16_type dim, uint16_type order )
{
    if ( dim < 2 )
        return 0;
    if ( order == 0 )
        return ( dim == 2 ) ? 1 : 0;
    if ( order < 3 )
        return 0;
    return ( order - 1 ) * ( order - 2 ) / 2;
}

/**
 * @brief Points per volume (interior volume points).
 */
constexpr uint16_type simplexPerVolume( uint16_type dim, uint16_type order )
{
    if ( dim < 3 )
        return 0;
    if ( order == 0 )
        return ( dim == 3 ) ? 1 : 0;
    if ( order < 4 )
        return 0;
    return ( order - 1 ) * ( order - 2 ) * ( order - 3 ) / 6;
}

//
// Runtime point calculations for dynamic-order hypercubes
//

/**
 * @brief Compute number of points for hypercube of given dimension and order.
 * Formula: (order+1)^dim
 */
constexpr size_type hypercubeTotal( uint16_type dim, uint16_type order )
{
    return powInt( static_cast<size_type>( order + 1 ), dim );
}

/**
 * @brief Points per vertex for hypercube.
 *
 * For 0D hypercube, Q0 has one point dof on the unique point/vertex.
 * For higher dimensions, Q0 has no vertex dofs.
 */
constexpr uint16_type hypercubePerVertex( uint16_type dim, uint16_type order )
{
    if ( dim == 0 )
        return 1;
    return ( order == 0 ) ? 0 : 1;
}

/**
 * @brief Backward-compatible overload (assumes dim > 0 semantics).
 */
constexpr uint16_type hypercubePerVertex( uint16_type order )
{
    return hypercubePerVertex( 1, order );
}

/**
 * @brief Points per edge (interior edge points) for hypercube.
 * Formula: order - 1 (for order >= 2), else special handling for order 0/1
 */
constexpr uint16_type hypercubePerEdge( uint16_type dim, uint16_type order )
{
    if ( dim < 1 )
        return 0;
    if ( order == 0 )
        return ( dim == 1 ) ? 1 : 0;  // Q0 has 1 point on the element
    if ( order < 2 )
        return 0;
    return order - 1;
}

/**
 * @brief Points per face (interior face points) for hypercube.
 * Formula: (order - 1)^2 for order >= 2, else special handling
 */
constexpr uint16_type hypercubePerFace( uint16_type dim, uint16_type order )
{
    if ( dim < 2 )
        return 0;
    if ( order == 0 )
        return ( dim == 2 ) ? 1 : 0;  // Q0 in 2D has 1 interior point
    if ( order < 2 )
        return 0;
    return static_cast<uint16_type>( ( order - 1 ) * ( order - 1 ) );
}

/**
 * @brief Points per volume (interior volume points) for hypercube.
 * Formula: (order - 1)^3 for order >= 2, else special handling
 */
constexpr uint16_type hypercubePerVolume( uint16_type dim, uint16_type order )
{
    if ( dim < 3 )
        return 0;
    if ( order == 0 )
        return ( dim == 3 ) ? 1 : 0;  // Q0 in 3D has 1 interior point
    if ( order < 2 )
        return 0;
    return static_cast<uint16_type>( ( order - 1 ) * ( order - 1 ) * ( order - 1 ) );
}

// Backward-compatible snake_case wrappers.
constexpr size_type pow_int( size_type base, int exp ) { return powInt( base, exp ); }

template <int Dim, int Order>
constexpr size_type simplex_dof_static()
{
    return simplexDofStatic<Dim, Order>();
}

template <int Dim, int Order>
constexpr size_type hypercube_dof_static()
{
    return hypercubeDofStatic<Dim, Order>();
}

template <int Dim>
inline size_type simplex_dof_dynamic( int order )
{
    return simplexDofDynamic<Dim>( order );
}

template <int Dim>
inline size_type hypercube_dof_dynamic( int order )
{
    return hypercubeDofDynamic<Dim>( order );
}

constexpr size_type simplex_num_points( uint16_type dim, uint16_type order )
{
    return simplexTotal( dim, order );
}

constexpr uint16_type simplex_pts_per_vertex( uint16_type dim, uint16_type order )
{
    return simplexPerVertex( dim, order );
}

constexpr uint16_type simplex_pts_per_vertex( uint16_type order )
{
    return simplexPerVertex( order );
}

constexpr uint16_type simplex_pts_per_edge( uint16_type dim, uint16_type order )
{
    return simplexPerEdge( dim, order );
}

constexpr uint16_type simplex_pts_per_face( uint16_type dim, uint16_type order )
{
    return simplexPerFace( dim, order );
}

constexpr uint16_type simplex_pts_per_volume( uint16_type dim, uint16_type order )
{
    return simplexPerVolume( dim, order );
}

constexpr size_type hypercube_num_points( uint16_type dim, uint16_type order )
{
    return hypercubeTotal( dim, order );
}

constexpr uint16_type hypercube_pts_per_vertex( uint16_type dim, uint16_type order )
{
    return hypercubePerVertex( dim, order );
}

constexpr uint16_type hypercube_pts_per_vertex( uint16_type order )
{
    return hypercubePerVertex( order );
}

constexpr uint16_type hypercube_pts_per_edge( uint16_type dim, uint16_type order )
{
    return hypercubePerEdge( dim, order );
}

constexpr uint16_type hypercube_pts_per_face( uint16_type dim, uint16_type order )
{
    return hypercubePerFace( dim, order );
}

constexpr uint16_type hypercube_pts_per_volume( uint16_type dim, uint16_type order )
{
    return hypercubePerVolume( dim, order );
}

} // namespace detail

/**
 * @brief DOF calculator for static or dynamic orders.
 */
template <int Dim, typename Order>
struct DofCalculator;

/**
 * @brief Static order specialization.
 */
template <int Dim, int Order>
struct DofCalculator<Dim, order_t<Order>>
{
    static constexpr size_type simplex()
    {
        return ::Feel::detail::simplexDofStatic<Dim, Order>();
    }

    static constexpr size_type hypercube()
    {
        return ::Feel::detail::hypercubeDofStatic<Dim, Order>();
    }
};

/**
 * @brief Dynamic order specialization.
 */
template <int Dim>
struct DofCalculator<Dim, order_t<Dynamic>>
{
    static size_type simplex( int order )
    {
        return ::Feel::detail::simplexDofDynamic<Dim>( order );
    }

    static size_type hypercube( int order )
    {
        return ::Feel::detail::hypercubeDofDynamic<Dim>( order );
    }
};

} // namespace Feel

#endif /* FEELPP_FEELPOLY_ORDER_HPP */
