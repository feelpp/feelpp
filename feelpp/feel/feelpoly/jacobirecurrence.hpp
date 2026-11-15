/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2005-2006 EPFL
 * SPDX-FileCopyrightText: 2006-2012 Joseph Fourier University
 * SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file jacobirecurrence.hpp
 * @brief Shared Jacobi recurrence for scalar and Eigen batch evaluation.
 */
#ifndef FEELPP_FEELPOLY_JACOBIRECURRENCE_HPP
#define FEELPP_FEELPOLY_JACOBIRECURRENCE_HPP 1
#include <feel/feelpoly/pointtabulation.hpp>

namespace Feel
{
namespace detail
{
/** @brief Coefficients in Pn=(constant+linear*x)Pn-1-previous*Pn-2. */
template<typename T>
struct JacobiCoefficients
{
    T constant; ///< Constant multiplier of Pn-1.
    T linear;   ///< Coordinate multiplier of Pn-1.
    T previous; ///< Multiplier of Pn-2.
};

/** @brief First-degree Jacobi coefficients, shared by scalar and batch paths. */
template<typename T>
[[nodiscard]] constexpr std::array<T, 2> jacobiFirstCoefficients( T a, T b )
{ return { ( a - b ) / T( 2 ), ( a + b + T( 2 ) ) / T( 2 ) }; }

/** @brief One normalized three-term recurrence, shared by every order mode. */
template<typename T>
[[nodiscard]] JacobiCoefficients<T> jacobiCoefficients( int n, T a, T b )
{
    T const k = T( n ), sum = a + b, twice = T( 2 ) * k + sum;
    T const denominator = T( 2 ) * k * ( k + sum ) * ( twice - T( 2 ) );
    return { ( twice - T( 1 ) ) * ( a * a - b * b ) / denominator,
             ( twice - T( 2 ) ) * ( twice - T( 1 ) ) * twice / denominator,
             T( 2 ) * ( k + a - T( 1 ) ) * ( k + b - T( 1 ) ) * twice / denominator };
}

/** @brief Rolling-state scalar evaluation; static degree remains a template input. */
template<int OrderSpec, typename T>
[[nodiscard]] T jacobiValue( int runtimeDegree, T a, T b, T x )
{
    int const degree = OrderSpec == Dynamic ? runtimeDegree : OrderSpec;
    if ( degree == 0 ) return T( 1 );
    auto const first = jacobiFirstCoefficients( a, b );
    T previous = T( 1 ), current = first[0] + first[1] * x;
    for ( int n = 2; n <= degree; ++n )
    {
        auto const c = jacobiCoefficients( n, a, b );
        T const next = ( c.constant + c.linear * x ) * current - c.previous * previous;
        previous = current;
        current = next;
    }
    return current;
}

/** @brief Derivative multiplier 2^-m (n+a+b+1)_m, shared by scalar/batch paths. */
template<unsigned DerivativeOrder, typename T>
[[nodiscard]] constexpr T jacobiDerivativeFactor( int n, T a, T b )
{
    T result = T( 1 );
    for ( unsigned j = 0; j < DerivativeOrder; ++j )
        result *= ( T( n ) + a + b + T( 1 ) + T( j ) ) / T( 2 );
    return result;
}

/** @brief Fill degree rows using one coordinate row, including empty batches. */
template<int OrderSpec, typename T, typename Output, typename X>
void jacobiRows( int runtimeDegree, T a, T b, Output& output, X const& x )
{
    int const degree = OrderSpec == Dynamic ? runtimeDegree : OrderSpec;
    output.row( 0 ).setOnes();
    if ( degree == 0 ) return;
    auto const first = jacobiFirstCoefficients( a, b );
    output.row( 1 ).array() = first[0] + first[1] * x;
    for ( int n = 2; n <= degree; ++n )
    {
        auto const c = jacobiCoefficients( n, a, b );
        output.row( n ).array() = ( c.constant + c.linear * x ) * output.row( n - 1 ).array()
                                  - c.previous * output.row( n - 2 ).array();
    }
}

/** @brief Fill a prepared value/derivative table through the same shifted recurrence.
 * @tparam DerivativeOrder Mathematical derivative order, including zero for values.
 * @tparam OrderSpec Compile-time degree or Dynamic.
 * @param runtimeDegree Effective degree, retained dynamically when OrderSpec is Dynamic.
 * @param a Jacobi weight exponent at x=1.
 * @param b Jacobi weight exponent at x=-1.
 * @param target Correctly sized writable table with degree+1 rows and sample columns.
 * @param x Coordinate array row, snapshotted or proven disjoint from target.
 */
template<unsigned DerivativeOrder, int OrderSpec, typename T, typename Output, typename X>
void jacobiTableRows( int runtimeDegree, T a, T b, Output& target, X const& x )
{
    static_assert( DerivativeOrder <= std::numeric_limits<uint16_type>::max() );
    int const degree = OrderSpec == Dynamic ? runtimeDegree : OrderSpec;
    if constexpr ( DerivativeOrder == 0 ) jacobiRows<OrderSpec>( degree, a, b, target, x );
    else
    {
        target.topRows( std::min( degree + 1, static_cast<int>( DerivativeOrder ) ) ).setZero();
        if constexpr ( OrderSpec == Dynamic || OrderSpec >= static_cast<int>( DerivativeOrder ) )
        {
            if ( degree < static_cast<int>( DerivativeOrder ) ) return;
            constexpr int shifted = OrderSpec == Dynamic ? Dynamic : OrderSpec - static_cast<int>( DerivativeOrder );
            auto lower = target.bottomRows( degree - DerivativeOrder + 1 );
            jacobiRows<shifted>( degree - DerivativeOrder, a + T( DerivativeOrder ), b + T( DerivativeOrder ), lower, x );
            for ( int n = DerivativeOrder; n <= degree; ++n )
                target.row( n ) *= jacobiDerivativeFactor<DerivativeOrder>( n, a, b );
        }
    }
}

} // namespace detail

/** @brief Jacobi polynomial Pn^(alpha,beta) with shared static/runtime kernels.
 * @tparam OrderSpec Compile-time degree or `Dynamic`.
 * @tparam T Scalar type.
 *
 * Single-point evaluation uses two rolling scalars. Batched values are degree
 * rows by sample columns. In-place batch evaluation reuses Eigen storage and
 * only snapshots points when an expression or output alias requires ownership.
 */
template<int OrderSpec, typename T = double>
class JacobiPolynomial : public OrderBase<OrderSpec>
{
    using super = OrderBase<OrderSpec>;
    static_assert( OrderSpec == Dynamic || OrderSpec <= std::numeric_limits<uint16_type>::max() );
public:
    /** Scalar type used by values and recurrence coefficients. */
    using value_type = T;
    /** Batch row extent, distinct from the polynomial degree. */
    static constexpr int rowExtent = OrderSpec == Dynamic ? Dynamic : OrderSpec + 1;
    /** Owning batch representation with bounded static extents. */
    template<int Samples = Dynamic, std::size_t InlineBytes = ::Feel::detail::defaultInlineMatrixBytes>
    using TableFor = ::Feel::detail::EigenMatrixBuffer<T, rowExtent, Samples, InlineBytes>;
    /** Table with an independent runtime sample count. */
    using Table = TableFor<>;

    /** @brief Construct a static-degree Jacobi polynomial, defaulting to Legendre. */
    explicit JacobiPolynomial( T alpha = T( 0 ), T beta = T( 0 ) )
        requires( OrderSpec != Dynamic ) : super(), M_alpha( alpha ), M_beta( beta ) {}
    /** @brief Construct or validate an effective degree for either order policy. */
    explicit JacobiPolynomial( RuntimeOrder degree, T alpha = T( 0 ), T beta = T( 0 ) )
        : super( degree ), M_alpha( alpha ), M_beta( beta ) {}
    /** @brief Effective polynomial degree. */
    [[nodiscard]] constexpr int degree() const noexcept
    { if constexpr ( OrderSpec == Dynamic ) return this->order(); else return OrderSpec; }
    /** @brief Exponent of the Jacobi weight factor (1-x). */
    [[nodiscard]] T alpha() const noexcept { return M_alpha; }
    /** @brief Exponent of the Jacobi weight factor (1+x). */
    [[nodiscard]] T beta() const noexcept { return M_beta; }
    /** @brief Evaluate Pn(x) using rolling state and the selected degree policy. */
    [[nodiscard]] T value( T x ) const
    { return ::Feel::detail::jacobiValue<OrderSpec>( degree(), M_alpha, M_beta, x ); }
    /** @brief Mathematical scalar evaluation syntax. */
    [[nodiscard]] T operator()( T x ) const { return value( x ); }
    /** @brief Evaluate an arbitrary derivative through the shifted Jacobi identity. */
    template<unsigned DerivativeOrder = 1>
    [[nodiscard]] T derivative( T x ) const
    {
        static_assert( DerivativeOrder <= std::numeric_limits<uint16_type>::max() );
        if constexpr ( OrderSpec != Dynamic && OrderSpec < static_cast<int>( DerivativeOrder ) ) return T( 0 );
        else
        {
            if ( degree() < static_cast<int>( DerivativeOrder ) ) return T( 0 );
            constexpr int shifted = OrderSpec == Dynamic ? Dynamic : OrderSpec - static_cast<int>( DerivativeOrder );
            return ::Feel::detail::jacobiDerivativeFactor<DerivativeOrder>( degree(), M_alpha, M_beta ) *
                ::Feel::detail::jacobiValue<shifted>( degree() - DerivativeOrder,
                    M_alpha + T( DerivativeOrder ), M_beta + T( DerivativeOrder ), x );
        }
    }
    /** @brief Allocate and evaluate P0 through Pn at row or column vector points. */
    template<typename Derived>
    [[nodiscard]] auto values( Eigen::MatrixBase<Derived> const& points ) const
    {
        TableFor<::Feel::detail::pointSampleExtent<Derived>()> result( degree() + 1, ::Feel::detail::pointCount( points ), degree() + 1, Dynamic );
        valuesInto( result, points );
        return result;
    }
    /** @brief Allocate and evaluate derivatives of P0 through Pn. */
    template<unsigned DerivativeOrder = 1, typename Derived>
    [[nodiscard]] auto derivatives( Eigen::MatrixBase<Derived> const& points ) const
    {
        TableFor<::Feel::detail::pointSampleExtent<Derived>()> result( degree() + 1, ::Feel::detail::pointCount( points ), degree() + 1, Dynamic );
        derivativesInto<DerivativeOrder>( result, points );
        return result;
    }
    /** @brief Evaluate values into reusable storage, supporting strided/aliased points. */
    template<::Feel::detail::EigenBuffer Buffer, typename Derived>
        requires std::same_as<typename Buffer::value_type, T>
    void valuesInto( Buffer& output, Eigen::MatrixBase<Derived> const& points ) const
    { tabulateInto<0>( output, points ); }
    /** @brief Evaluate derivatives into reusable storage with the same alias contract. */
    template<unsigned DerivativeOrder = 1, ::Feel::detail::EigenBuffer Buffer, typename Derived>
        requires std::same_as<typename Buffer::value_type, T>
    void derivativesInto( Buffer& output, Eigen::MatrixBase<Derived> const& points ) const
    { tabulateInto<DerivativeOrder>( output, points ); }

    /** @brief Fill a correctly sized Eigen destination, including a mapped legacy owner. */
    template<typename Output, typename Derived>
        requires( ( Output::Flags & Eigen::DirectAccessBit ) != 0 )
    void valuesInto( Eigen::MatrixBase<Output>& output, Eigen::MatrixBase<Derived> const& points ) const
    {
        ::Feel::detail::EigenTabulationDestination<Output> destination( output );
        tabulateInto<0>( destination, points );
    }

    /** @brief Fill a caller-owned derivative table with the same alias and shape checks. */
    template<unsigned DerivativeOrder = 1, typename Output, typename Derived>
        requires( ( Output::Flags & Eigen::DirectAccessBit ) != 0 )
    void derivativesInto( Eigen::MatrixBase<Output>& output, Eigen::MatrixBase<Derived> const& points ) const
    {
        ::Feel::detail::EigenTabulationDestination<Output> destination( output );
        tabulateInto<DerivativeOrder>( destination, points );
    }

private:
    /** @brief Shared orientation, ownership and shape handling for both batch operations. */
    template<unsigned DerivativeOrder, typename Buffer, typename Derived>
    void tabulateInto( Buffer& output, Eigen::MatrixBase<Derived> const& points ) const
    {
        static_assert( DerivativeOrder <= std::numeric_limits<uint16_type>::max() );
        ::Feel::detail::tabulatePointRows( output, points, degree() + 1,
            [&]( auto& target, auto const& x )
            {
                ::Feel::detail::jacobiTableRows<DerivativeOrder, OrderSpec>( degree(), M_alpha, M_beta, target, x );
            } );
    }
    /** Jacobi weight exponent at x=1. */
    T M_alpha;
    /** Jacobi weight exponent at x=-1. */
    T M_beta;
};
} // namespace Feel
#endif // FEELPP_FEELPOLY_JACOBIRECURRENCE_HPP
