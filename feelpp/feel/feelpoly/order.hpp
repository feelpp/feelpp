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
struct OrderTag
{
    /** Unmodified compile-time order specification. */
    static constexpr int order_spec = N;
    /** Compile-time order carried by the tag. */
    static constexpr int value = order_spec;
    /** Whether this tag denotes runtime order. */
    static constexpr bool is_dynamic = (order_spec == Dynamic);
    /** Whether this tag denotes compile-time order. */
    static constexpr bool is_static = !is_dynamic;
};

/** @brief Compatibility spelling of OrderTag for existing clients. */
template<int N>
using order_t = OrderTag<N>;

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
    /** Authoritative, unmodified template order specification. */
    static constexpr int order_spec = Order;
    /** Compatibility alias for the compile-time polynomial order. */
    static constexpr int nOrder = order_spec;
    /** Static-order policies never carry runtime order. */
    static constexpr bool is_order_dynamic = false;
    /** Static-order policies expose their order at compile time. */
    static constexpr bool is_order_static = true;

    /** @brief Construct the compile-time order policy. */
    constexpr OrderBase() noexcept = default;
    /**
     * @brief Validate a runtime tag against the compile-time order.
     *
     * This overload permits generic policy-based construction without silently
     * accepting a configuration value that disagrees with the static type.
     * Static order remains storage-free and is still the single source of
     * truth.
     *
     * @throws std::invalid_argument when the supplied order differs from
     * `Order`.
     */
    constexpr explicit OrderBase( RuntimeOrder runtimeOrder )
    {
        if ( runtimeOrder.value != static_cast<uint16_type>( Order ) )
            throw std::invalid_argument(
                "RuntimeOrder disagrees with compile-time order" );
    }

    /** @brief Return the compile-time order as the common unsigned order type. */
    [[nodiscard]] static constexpr uint16_type staticOrder() noexcept
    {
        return static_cast<uint16_type>( order_spec );
    }
    /** @brief Return the effective polynomial order. */
    [[nodiscard]] constexpr uint16_type order() const noexcept { return staticOrder(); }
};

/**
 * @brief Dynamic-order specialization of OrderBase.
 */
template <>
class OrderBase<Dynamic>
{
  public:
    /** Authoritative, unmodified template order specification. */
    static constexpr int order_spec = Dynamic;
    /** Compatibility alias for the dynamic-order sentinel. */
    static constexpr int nOrder = order_spec;
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
 * @tparam T Order provider, optionally cv/ref-qualified.
 */
template <typename T>
inline constexpr bool hasOrderDynamicFlag = requires {
    std::remove_cvref_t<T>::is_order_dynamic;
};

/**
 * @brief Query whether a type uses runtime polynomial order.
 * @tparam T Order provider, optionally cv/ref-qualified.
 *
 * Types not exposing `is_order_dynamic` are considered static by default.
 */
template <typename T>
inline constexpr bool orderIsDynamic = []
{
    if constexpr ( hasOrderDynamicFlag<T> )
        return bool( std::remove_cvref_t<T>::is_order_dynamic );
    else
        return false;
}();

/**
 * @brief Query whether a type uses compile-time polynomial order.
 * @tparam T Order provider, optionally cv/ref-qualified.
 */
template <typename T>
inline constexpr bool orderIsStatic = !orderIsDynamic<T>;

/**
 * @brief Return the public compile-time order exposed by an order-aware type.
 *
 * Public FE shells use `nOrder_v`, common order policies use `order_spec`, and
 * legacy static providers expose only `nOrder`. Keeping that compatibility
 * priority here prevents consumers from repeating it.
 */
template <typename T>
[[nodiscard]] consteval int
staticOrderValue()
{
    using type = std::remove_cvref_t<T>;
    static_assert( orderIsStatic<type>,
                   "staticOrderValue requires a compile-time-order type" );

    if constexpr ( requires { type::nOrder_v; } )
        return type::nOrder_v;
    else if constexpr ( requires { type::order_spec; } )
        return type::order_spec;
    else
    {
        static_assert( requires { type::nOrder; },
                       "static order provider exposes no order metadata" );
        return type::nOrder;
    }
}

/** @brief Return the signed specification without erasing runtime identity. */
template<class T>
[[nodiscard]] consteval int orderSpecification()
{
    if constexpr ( orderIsDynamic<T> ) return Dynamic;
    else return staticOrderValue<T>();
}

/** @brief Return the default construction order for a static/dynamic type. */
template <typename T>
[[nodiscard]] constexpr RuntimeOrder
defaultOrder() noexcept
{
    if constexpr ( orderIsDynamic<T> )
        return RuntimeOrder{ 0 };
    else
        return RuntimeOrder{
            static_cast<uint16_type>( staticOrderValue<T>() ) };
}

/**
 * @brief Validate an effective order against an order-aware type.
 *
 * Dynamic types retain the supplied value. Static types require equality with
 * their public compile-time order and remain storage-free.
 */
template <typename T>
[[nodiscard]] constexpr RuntimeOrder
validatedOrder( RuntimeOrder order )
{
    if constexpr ( orderIsStatic<T> )
    {
        if ( order.value !=
             static_cast<uint16_type>( staticOrderValue<T>() ) )
            throw std::invalid_argument(
                "RuntimeOrder disagrees with compile-time order" );
    }
    return order;
}

/**
 * @brief Construct an order-aware object without duplicating static/dynamic
 * dispatch.
 *
 * Runtime-order types receive the effective order explicitly. Compile-time
 * types retain their default construction and therefore carry no runtime
 * storage or branch.
 *
 * @tparam T object type exposing the common order policy
 * @param order effective order, validated for compile-time types
 * @return an object initialized consistently with its order policy
 */
template <typename T>
[[nodiscard]] T
makeOrderAware( RuntimeOrder order )
{
    if constexpr ( orderIsDynamic<T> )
        return T{ order };
    else
    {
        (void)validatedOrder<T>( order );
        return T{};
    }
}

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

/** Shape selection for the shared nodal-cardinality kernel. */
enum class NodalCellKind
{
    Simplex,
    Hypercube
};

/**
 * @brief Stateless cardinality formulas shared by static and dynamic layouts.
 *
 * Shape-independent validation, P0 placement, and entity dispatch live once
 * here. The `if constexpr` branches contain only the two genuinely different
 * mathematical formulas.
 */
template<NodalCellKind CellKind>
struct NodalCardinality
{
    [[nodiscard]] static constexpr size_type total( int cellDim, int order ) noexcept
    {
        if ( cellDim < 0 || order < 0 )
            return 0;
        if constexpr ( CellKind == NodalCellKind::Simplex )
            return binomial( cellDim + order, cellDim );
        else
            return powInt( static_cast<size_type>( order + 1 ), cellDim );
    }

    [[nodiscard]] static constexpr uint16_type onEntity(
        int cellDim, int entityDim, int order ) noexcept
    {
        if ( cellDim < 0 || entityDim < 0 || entityDim > cellDim || order < 0 )
            return 0;
        if ( order == 0 )
            return entityDim == cellDim ? 1 : 0;
        if ( entityDim == 0 )
            return 1;

        if constexpr ( CellKind == NodalCellKind::Simplex )
        {
            if ( order <= entityDim )
                return 0;
            return static_cast<uint16_type>( binomial( order - 1, entityDim ) );
        }
        else
        {
            if ( order < 2 )
                return 0;
            return static_cast<uint16_type>(
                powInt( static_cast<size_type>( order - 1 ), entityDim ) );
        }
    }

    [[nodiscard]] static constexpr uint16_type perVertex( int cellDim, int order ) noexcept
    {
        return onEntity( cellDim, 0, order );
    }
    [[nodiscard]] static constexpr uint16_type perEdge( int cellDim, int order ) noexcept
    {
        return onEntity( cellDim, 1, order );
    }
    [[nodiscard]] static constexpr uint16_type perFace( int cellDim, int order ) noexcept
    {
        return onEntity( cellDim, 2, order );
    }
    [[nodiscard]] static constexpr uint16_type perVolume( int cellDim, int order ) noexcept
    {
        return onEntity( cellDim, 3, order );
    }
};

using SimplexCardinality = NodalCardinality<NodalCellKind::Simplex>;
using HypercubeCardinality = NodalCardinality<NodalCellKind::Hypercube>;

/** Largest point/DoF cardinality addressable by the local-index API. */
inline constexpr size_type maxLocalIndexCardinality =
    std::numeric_limits<uint16_type>::max();

/**
 * @brief Narrow a derived count only after checking the local-index range.
 * @throws std::length_error when `value` cannot be represented by
 * `uint16_type`.
 */
[[nodiscard]] constexpr uint16_type
checkedLocalIndexValue( size_type value )
{
    if ( value > maxLocalIndexCardinality )
        throw std::length_error(
            "Derived point/DoF count exceeds the uint16_type local-index limit" );
    return static_cast<uint16_type>( value );
}

/**
 * @brief Validate that a nodal order fits the local point/DoF index model.
 */
template <typename Cardinality>
[[nodiscard]] constexpr RuntimeOrder
validatedNodalOrder( int cellDim, RuntimeOrder order,
                     size_type components = 1 )
{
    (void)checkedLocalIndexValue(
        components * Cardinality::total( cellDim, order.value ) );
    return order;
}

template <int Dim, int Order>
constexpr size_type simplexDofStatic()
{
    static_assert( Dim >= 0, "Dimension must be non-negative" );
    static_assert( Order >= 0, "Order must be non-negative" );
    return SimplexCardinality::total( Dim, Order );
}

template <int Dim, int Order>
constexpr size_type hypercubeDofStatic()
{
    static_assert( Dim >= 0, "Dimension must be non-negative" );
    static_assert( Order >= 0, "Order must be non-negative" );
    return HypercubeCardinality::total( Dim, Order );
}

template <int Dim>
constexpr size_type simplexDofDynamic( int order ) noexcept
{
    return SimplexCardinality::total( Dim, order );
}

template <int Dim>
constexpr size_type hypercubeDofDynamic( int order ) noexcept
{
    return HypercubeCardinality::total( Dim, order );
}

constexpr size_type simplexTotal( uint16_type dim, uint16_type order ) noexcept
{
    return SimplexCardinality::total( dim, order );
}

constexpr uint16_type simplexPerVertex( uint16_type dim, uint16_type order ) noexcept
{
    return SimplexCardinality::perVertex( dim, order );
}

constexpr uint16_type simplexPerVertex( uint16_type order ) noexcept
{
    return SimplexCardinality::perVertex( 1, order );
}

constexpr uint16_type simplexPerEdge( uint16_type dim, uint16_type order ) noexcept
{
    return SimplexCardinality::perEdge( dim, order );
}

constexpr uint16_type simplexPerFace( uint16_type dim, uint16_type order ) noexcept
{
    return SimplexCardinality::perFace( dim, order );
}

constexpr uint16_type simplexPerVolume( uint16_type dim, uint16_type order ) noexcept
{
    return SimplexCardinality::perVolume( dim, order );
}

constexpr size_type hypercubeTotal( uint16_type dim, uint16_type order ) noexcept
{
    return HypercubeCardinality::total( dim, order );
}

constexpr uint16_type hypercubePerVertex( uint16_type dim, uint16_type order ) noexcept
{
    return HypercubeCardinality::perVertex( dim, order );
}

constexpr uint16_type hypercubePerVertex( uint16_type order ) noexcept
{
    return HypercubeCardinality::perVertex( 1, order );
}

constexpr uint16_type hypercubePerEdge( uint16_type dim, uint16_type order ) noexcept
{
    return HypercubeCardinality::perEdge( dim, order );
}

constexpr uint16_type hypercubePerFace( uint16_type dim, uint16_type order ) noexcept
{
    return HypercubeCardinality::perFace( dim, order );
}

constexpr uint16_type hypercubePerVolume( uint16_type dim, uint16_type order ) noexcept
{
    return HypercubeCardinality::perVolume( dim, order );
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
struct DofCalculator
{
    static_assert( requires {
        Order::value;
        Order::is_static;
        Order::is_dynamic;
    }, "DofCalculator requires an order tag" );

    static constexpr size_type simplex()
        requires( Order::is_static )
    {
        return ::Feel::detail::SimplexCardinality::total( Dim, Order::value );
    }

    static constexpr size_type hypercube()
        requires( Order::is_static )
    {
        return ::Feel::detail::HypercubeCardinality::total( Dim, Order::value );
    }

    static constexpr size_type simplex( int order ) noexcept
        requires( Order::is_dynamic )
    {
        return ::Feel::detail::SimplexCardinality::total( Dim, order );
    }

    static constexpr size_type hypercube( int order ) noexcept
        requires( Order::is_dynamic )
    {
        return ::Feel::detail::HypercubeCardinality::total( Dim, order );
    }
};

} // namespace Feel

#endif /* FEELPP_FEELPOLY_ORDER_HPP */
