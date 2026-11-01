/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file bernsteinbezier.hpp
 * @brief Bernstein/Bezier geometry representation on supported reference cells.
 */
#if !defined( FEELPP_FEELPOLY_BERNSTEINBEZIER_HPP )
#define FEELPP_FEELPOLY_BERNSTEINBEZIER_HPP 1

#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelpoly/geometrycoefficients.hpp>
#include <feel/feelpoly/order.hpp>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace Feel::geometry
{
namespace detail
{
/** @brief Evaluate an integer binomial coefficient in the scalar type. */
template <typename T>
[[nodiscard]] constexpr T binomialCoefficient( uint16_type n, uint16_type k )
{
    if ( k > n )
        return T( 0 );
    k = std::min<uint16_type>( k, n - k );
    T result = T( 1 );
    for ( uint16_type i = 1; i <= k; ++i )
        result = result * T( n - k + i ) / T( i );
    return result;
}

/** @brief Evaluate a non-negative integer power without a generic `pow` call. */
template <typename T>
[[nodiscard]] constexpr T integerPower( T base, int exponent )
{
    T result = T( 1 );
    for ( int i = 0; i < exponent; ++i )
        result *= base;
    return result;
}

/** @brief Evaluate a univariate Bernstein polynomial on `[0,1]`. */
template <typename T>
[[nodiscard]] T bernstein1d( int index, int degree, T t )
{
    if ( index < 0 || index > degree || degree < 0 )
        return T( 0 );
    return binomialCoefficient<T>( static_cast<uint16_type>( degree ),
                                   static_cast<uint16_type>( index ) ) *
           integerPower( t, index ) * integerPower( T( 1 ) - t, degree - index );
}

/** @brief Enumerate simplex multi-indices of total degree `order`. */
template <int Dim>
[[nodiscard]] std::vector<std::array<uint16_type, Dim + 1>>
simplexMultiIndices( uint16_type order )
{
    using index_type = std::array<uint16_type, Dim + 1>;
    std::vector<index_type> result;
    index_type current{};
    auto enumerate = [&]( auto&& self, int coordinate, uint16_type remaining ) -> void
    {
        if ( coordinate == Dim )
        {
            current[coordinate] = remaining;
            result.push_back( current );
            return;
        }
        for ( int value = remaining; value >= 0; --value )
        {
            current[coordinate] = static_cast<uint16_type>( value );
            self( self, coordinate + 1,
                  static_cast<uint16_type>( remaining - value ) );
        }
    };
    enumerate( enumerate, 0, order );
    return result;
}

/** @brief Enumerate tensor-product multi-indices bounded by `order`. */
template <int Dim>
[[nodiscard]] std::vector<std::array<uint16_type, Dim>>
hypercubeMultiIndices( uint16_type order )
{
    using index_type = std::array<uint16_type, Dim>;
    std::vector<index_type> result;
    index_type current{};
    auto enumerate = [&]( auto&& self, int coordinate ) -> void
    {
        if ( coordinate == Dim )
        {
            result.push_back( current );
            return;
        }
        for ( uint16_type value = 0; value <= order; ++value )
        {
            current[coordinate] = value;
            self( self, coordinate + 1 );
        }
    };
    enumerate( enumerate, 0 );
    return result;
}

/** @brief Evaluate one simplex Bernstein polynomial from barycentric values. */
template <int Dim, typename T>
[[nodiscard]] T simplexBernstein(
    std::array<int, Dim + 1> const& alpha, int degree,
    Eigen::Matrix<T, Dim + 1, 1> const& lambda )
{
    if ( degree < 0 )
        return T( 0 );
    int sum = 0;
    T coefficient = T( 1 );
    int remaining = degree;
    for ( int i = 0; i < Dim + 1; ++i )
    {
        if ( alpha[i] < 0 )
            return T( 0 );
        sum += alpha[i];
        coefficient *= binomialCoefficient<T>(
            static_cast<uint16_type>( remaining ),
            static_cast<uint16_type>( alpha[i] ) );
        remaining -= alpha[i];
    }
    if ( sum != degree )
        return T( 0 );
    T result = coefficient;
    for ( int i = 0; i < Dim + 1; ++i )
        result *= integerPower( lambda( i ), alpha[i] );
    return result;
}
} // namespace detail

//!
//! \class BernsteinBezier
//! \brief Bernstein basis and Bezier control-point representation.
//!
//! \tparam Dim reference dimension
//! \tparam Order static polynomial order or `Dynamic`
//! \tparam Entity reference-entity family
//! \tparam T scalar type
template <int Dim, int Order = 1,
          template <int, int, int> class Entity = Simplex, typename T = double>
class BernsteinBezier : public OrderBase<Order>
{
    static_assert( Order == Dynamic || Order >= 1,
                   "A geometry representation order must be at least one" );
    static_assert( Feel::is_simplex_v<Entity<Dim, 1, Dim>> ||
                       Feel::is_hypercube_v<Entity<Dim, 1, Dim>>,
                   "BernsteinBezier currently supports simplex and hypercube cells only" );
    using order_base_type = OrderBase<Order>;
    using multi_index_type = std::conditional_t<
        Feel::is_simplex_v<Entity<Dim, 1, Dim>>,
        std::array<uint16_type, Dim + 1>, std::array<uint16_type, Dim>>;
    using multi_indices_type = std::vector<multi_index_type>;

  public:
    /** Scalar value type. */
    using value_type = T;
    /** Reference convex used by the representation. */
    using convex_type = Entity<Dim, ( Order == Dynamic ? 1 : Order ), Dim>;
    /** Matrix whose columns are reference points. */
    using points_type = Eigen::Matrix<T, Dim, Eigen::Dynamic>;
    /** Matrix whose rows are basis functions and columns are points. */
    using basis_matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    /** @brief Return the mathematical coefficient semantics. */
    [[nodiscard]] static consteval CoefficientSemantics coefficientSemantics() noexcept
    {
        return CoefficientSemantics::ControlPoint;
    }
    /** @brief Return the reference dimension. */
    [[nodiscard]] static consteval int dimension() noexcept { return Dim; }
    /** @brief Return whether the polynomial order is selected at runtime. */
    [[nodiscard]] static consteval bool isOrderDynamic() noexcept { return Order == Dynamic; }
    /** @brief Return whether the polynomial order is known at compile time. */
    [[nodiscard]] static consteval bool isOrderStatic() noexcept { return !isOrderDynamic(); }
    /** @brief Return the template order, including `Dynamic`. */
    [[nodiscard]] static consteval int staticOrder() noexcept { return Order; }
    /** @brief Return whether exact de Casteljau subdivision is supported. */
    [[nodiscard]] static consteval bool supportsSubdivision() noexcept { return Dim == 1; }
    /** @brief Return whether coefficient convex-hull bounds are available. */
    [[nodiscard]] static consteval bool supportsConvexHull() noexcept { return true; }
    /** @brief Return whether the reference entity is a simplex. */
    [[nodiscard]] static consteval bool isSimplex() noexcept
    {
        return Feel::is_simplex_v<Entity<Dim, 1, Dim>>;
    }
    /** @brief Return whether the reference entity is a hypercube. */
    [[nodiscard]] static consteval bool isHypercube() noexcept
    {
        return Feel::is_hypercube_v<Entity<Dim, 1, Dim>>;
    }

    /** @brief Construct a statically ordered representation. */
    BernsteinBezier()
        requires( Order != Dynamic )
        : M_multiIndices( makeMultiIndices( Order ) )
    {
    }

    /**
     * @brief Construct a representation from a runtime order.
     * @param order requested polynomial order
     * @throws std::invalid_argument for order zero or a static-order mismatch
     */
    explicit BernsteinBezier( RuntimeOrder order )
        : order_base_type( order ), M_multiIndices( makeMultiIndices( order.value ) )
    {
        if ( order.value < 1 )
            throw std::invalid_argument( "A geometry representation order must be at least one" );
        if constexpr ( Order != Dynamic )
            if ( order.value != Order )
                throw std::invalid_argument( "Runtime order disagrees with static geometry representation" );
    }

    /** @brief Return the effective polynomial order. */
    [[nodiscard]] constexpr uint16_type order() const noexcept
    {
        return order_base_type::order();
    }

    /** @brief Return the number of control coefficients. */
    [[nodiscard]] size_type coefficientCount() const noexcept
    {
        if constexpr ( isSimplex() )
            return Feel::detail::simplexTotal( Dim, this->order() );
        else if constexpr ( isHypercube() )
            return Feel::detail::hypercubeTotal( Dim, this->order() );
    }

    /** @brief Return a stable human-readable representation name. */
    [[nodiscard]] static constexpr std::string_view name() noexcept
    {
        return "BernsteinBezier";
    }

    /**
     * @brief Evaluate every Bernstein basis function.
     * @param points reference points stored by column
     * @return basis-by-point evaluation matrix
     */
    [[nodiscard]] basis_matrix_type basisValues( points_type const& points ) const
    {
        if constexpr ( isSimplex() )
            return this->simplexValues( points );
        else
            return this->hypercubeValues( points );
    }

    /**
     * @brief Evaluate reference gradients of every basis function.
     * @param points reference points stored by column
     * @return one basis-by-point matrix per reference derivative
     */
    [[nodiscard]] std::array<basis_matrix_type, Dim>
    basisGradients( points_type const& points ) const
    {
        if constexpr ( isSimplex() )
            return this->simplexGradients( points );
        else
            return this->hypercubeGradients( points );
    }

    /**
     * @brief Evaluate reference Hessians of every basis function.
     * @param points reference points stored by column
     * @return row-major `(i,j)` array of basis-by-point derivative matrices
     */
    [[nodiscard]] std::array<basis_matrix_type, Dim * Dim>
    basisHessians( points_type const& points ) const
    {
        if constexpr ( isSimplex() )
            return this->simplexHessians( points );
        else
            return this->hypercubeHessians( points );
    }

    /**
     * @brief Return control-point indices retained on a simplex facet.
     * @param facet facet identified by its vanishing barycentric coordinate
     */
    [[nodiscard]] std::vector<size_type>
    facetCoefficientIndices( SimplexFacet facet ) const
        requires( isSimplex() )
    {
        if ( facet.omittedBarycentricCoordinate > Dim )
            throw std::out_of_range( "Invalid simplex barycentric coordinate" );
        std::vector<size_type> result;
        for ( size_type i = 0; i < M_multiIndices.size(); ++i )
            if ( M_multiIndices[i][facet.omittedBarycentricCoordinate] == 0 )
                result.push_back( i );
        return result;
    }

    /**
     * @brief Return control-point indices retained on a hypercube facet.
     * @param facet coordinate axis and lower/upper endpoint
     */
    [[nodiscard]] std::vector<size_type>
    facetCoefficientIndices( HypercubeFacet facet ) const
        requires( !isSimplex() )
    {
        if ( facet.axis >= Dim )
            throw std::out_of_range( "Invalid hypercube facet axis" );
        std::vector<size_type> result;
        uint16_type const endpoint = facet.upper ? this->order() : 0;
        for ( size_type i = 0; i < M_multiIndices.size(); ++i )
            if ( M_multiIndices[i][facet.axis] == endpoint )
                result.push_back( i );
        return result;
    }

  private:
    /** Evaluate all simplex Bernstein basis functions. */
    [[nodiscard]] basis_matrix_type simplexValues( points_type const& points ) const
    {
        basis_matrix_type result( M_multiIndices.size(), points.cols() );
        for ( Eigen::Index q = 0; q < points.cols(); ++q )
        {
            Eigen::Matrix<T, Dim + 1, 1> lambda;
            lambda( 0 ) = T( 1 );
            for ( int d = 0; d < Dim; ++d )
            {
                lambda( d + 1 ) = ( points( d, q ) + T( 1 ) ) / T( 2 );
                lambda( 0 ) -= lambda( d + 1 );
            }
            for ( size_type a = 0; a < M_multiIndices.size(); ++a )
            {
                std::array<int, Dim + 1> alpha{};
                for ( int i = 0; i < Dim + 1; ++i )
                    alpha[i] = M_multiIndices[a][i];
                result( a, q ) = detail::simplexBernstein<Dim>( alpha, this->order(), lambda );
            }
        }
        return result;
    }

    /** Evaluate all simplex Bernstein first derivatives. */
    [[nodiscard]] std::array<basis_matrix_type, Dim>
    simplexGradients( points_type const& points ) const
    {
        std::array<basis_matrix_type, Dim> result;
        for ( auto& derivative : result )
            derivative.resize( M_multiIndices.size(), points.cols() );
        T const factor = T( this->order() );
        for ( Eigen::Index q = 0; q < points.cols(); ++q )
        {
            Eigen::Matrix<T, Dim + 1, 1> lambda;
            lambda( 0 ) = T( 1 );
            for ( int d = 0; d < Dim; ++d )
            {
                lambda( d + 1 ) = ( points( d, q ) + T( 1 ) ) / T( 2 );
                lambda( 0 ) -= lambda( d + 1 );
            }
            for ( size_type a = 0; a < M_multiIndices.size(); ++a )
                for ( int d = 0; d < Dim; ++d )
                {
                    std::array<int, Dim + 1> lower{};
                    for ( int i = 0; i < Dim + 1; ++i )
                        lower[i] = M_multiIndices[a][i];
                    --lower[d + 1];
                    T const positive = detail::simplexBernstein<Dim>(
                        lower, this->order() - 1, lambda );
                    ++lower[d + 1];
                    --lower[0];
                    T const negative = detail::simplexBernstein<Dim>(
                        lower, this->order() - 1, lambda );
                    result[d]( a, q ) = factor * ( positive - negative ) / T( 2 );
                }
        }
        return result;
    }

    /** Evaluate all simplex Bernstein second derivatives. */
    [[nodiscard]] std::array<basis_matrix_type, Dim * Dim>
    simplexHessians( points_type const& points ) const
    {
        std::array<basis_matrix_type, Dim * Dim> result;
        for ( auto& derivative : result )
            derivative.resize( M_multiIndices.size(), points.cols() );
        T const factor = T( this->order() ) * T( this->order() - 1 ) / T( 4 );
        for ( Eigen::Index q = 0; q < points.cols(); ++q )
        {
            Eigen::Matrix<T, Dim + 1, 1> lambda;
            lambda( 0 ) = T( 1 );
            for ( int d = 0; d < Dim; ++d )
            {
                lambda( d + 1 ) = ( points( d, q ) + T( 1 ) ) / T( 2 );
                lambda( 0 ) -= lambda( d + 1 );
            }
            for ( size_type a = 0; a < M_multiIndices.size(); ++a )
                for ( int d = 0; d < Dim; ++d )
                    for ( int e = 0; e < Dim; ++e )
                    {
                        std::array<int, Dim + 1> alpha{};
                        for ( int i = 0; i < Dim + 1; ++i )
                            alpha[i] = M_multiIndices[a][i];
                        auto term = [&]( int first, int second )
                        {
                            auto lower = alpha;
                            --lower[first];
                            --lower[second];
                            return detail::simplexBernstein<Dim>(
                                lower, this->order() - 2, lambda );
                        };
                        result[d * Dim + e]( a, q ) = factor * ( term( d + 1, e + 1 ) - term( d + 1, 0 ) -
                                                                 term( 0, e + 1 ) + term( 0, 0 ) );
                    }
        }
        return result;
    }

    /** Evaluate all tensor-product Bernstein basis functions. */
    [[nodiscard]] basis_matrix_type hypercubeValues( points_type const& points ) const
    {
        basis_matrix_type result( M_multiIndices.size(), points.cols() );
        for ( Eigen::Index q = 0; q < points.cols(); ++q )
            for ( size_type a = 0; a < M_multiIndices.size(); ++a )
            {
                T value = T( 1 );
                for ( int d = 0; d < Dim; ++d )
                    value *= detail::bernstein1d<T>(
                        M_multiIndices[a][d], this->order(),
                        ( points( d, q ) + T( 1 ) ) / T( 2 ) );
                result( a, q ) = value;
            }
        return result;
    }

    /** Evaluate all tensor-product Bernstein first derivatives. */
    [[nodiscard]] std::array<basis_matrix_type, Dim>
    hypercubeGradients( points_type const& points ) const
    {
        std::array<basis_matrix_type, Dim> result;
        for ( auto& derivative : result )
            derivative.resize( M_multiIndices.size(), points.cols() );
        for ( Eigen::Index q = 0; q < points.cols(); ++q )
            for ( size_type a = 0; a < M_multiIndices.size(); ++a )
                for ( int d = 0; d < Dim; ++d )
                {
                    T value = T( this->order() ) / T( 2 );
                    for ( int e = 0; e < Dim; ++e )
                    {
                        T const t = ( points( e, q ) + T( 1 ) ) / T( 2 );
                        if ( e == d )
                            value *= detail::bernstein1d<T>( M_multiIndices[a][e] - 1,
                                                             this->order() - 1, t ) -
                                     detail::bernstein1d<T>( M_multiIndices[a][e],
                                                             this->order() - 1, t );
                        else
                            value *= detail::bernstein1d<T>(
                                M_multiIndices[a][e], this->order(), t );
                    }
                    result[d]( a, q ) = value;
                }
        return result;
    }

    /** Evaluate all tensor-product Bernstein second derivatives. */
    [[nodiscard]] std::array<basis_matrix_type, Dim * Dim>
    hypercubeHessians( points_type const& points ) const
    {
        std::array<basis_matrix_type, Dim * Dim> result;
        for ( auto& derivative : result )
            derivative.resize( M_multiIndices.size(), points.cols() );
        for ( Eigen::Index q = 0; q < points.cols(); ++q )
            for ( size_type a = 0; a < M_multiIndices.size(); ++a )
                for ( int d = 0; d < Dim; ++d )
                    for ( int e = 0; e < Dim; ++e )
                    {
                        T value = d == e
                                      ? T( this->order() ) * T( this->order() - 1 ) / T( 4 )
                                      : T( this->order() ) * T( this->order() ) / T( 4 );
                        for ( int f = 0; f < Dim; ++f )
                        {
                            T const t = ( points( f, q ) + T( 1 ) ) / T( 2 );
                            int const index = M_multiIndices[a][f];
                            if ( f == d && f == e )
                                value *= detail::bernstein1d<T>( index - 2, this->order() - 2, t ) -
                                         T( 2 ) * detail::bernstein1d<T>( index - 1, this->order() - 2, t ) +
                                         detail::bernstein1d<T>( index, this->order() - 2, t );
                            else if ( f == d || f == e )
                                value *= detail::bernstein1d<T>( index - 1, this->order() - 1, t ) -
                                         detail::bernstein1d<T>( index, this->order() - 1, t );
                            else
                                value *= detail::bernstein1d<T>( index, this->order(), t );
                        }
                        result[d * Dim + e]( a, q ) = value;
                    }
        return result;
    }

    /** @brief Build the immutable multi-index table for an effective order. */
    [[nodiscard]] static multi_indices_type makeMultiIndices( uint16_type order )
    {
        if constexpr ( Feel::is_simplex_v<Entity<Dim, 1, Dim>> )
            return detail::simplexMultiIndices<Dim>( order );
        else if constexpr ( isHypercube() )
            return detail::hypercubeMultiIndices<Dim>( order );
    }

    /** Basis multi-indices reused by all evaluations. */
    multi_indices_type M_multiIndices;
};

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_BERNSTEINBEZIER_HPP
