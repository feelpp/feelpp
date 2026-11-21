/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2005-2006 EPFL
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file expansionkernels.hpp
 * @brief Eigen coordinate maps and shared degree-aware simplex scaling powers.
 */
#ifndef FEELPP_FEELPOLY_EXPANSIONKERNELS_HPP
#define FEELPP_FEELPOLY_EXPANSIONKERNELS_HPP 1
#include <feel/feelpoly/pointtabulation.hpp>
#include <cmath>
#include <limits>

namespace Feel
{
/** @brief Cartesian/collapsed simplex maps with one implementation for all batch sizes.
 * @tparam Dim Reference dimension, two or three.
 * @tparam T Scalar type.
 * @note Legacy collapsed-vertex representatives are retained for quadrature
 * and coordinate-map clients. The simplex modal kernel uses regular coordinates.
 */
template<int Dim, typename T = double>
class SimplexCoordinates
{
    static_assert( Dim == 2 || Dim == 3 );
public:
    /** Fixed spatial rows and independent sample columns. */
    using Points = ::Feel::detail::EigenMatrixBuffer<T, Dim, Dynamic>;
    /** Tolerance matching the established collapsed-coordinate convention. */
    static constexpr T boundaryTolerance = T( 64 ) * std::numeric_limits<T>::epsilon();

    /** @brief Map Cartesian simplex points to collapsed coordinates. */
    template<typename Derived>
    [[nodiscard]] static auto toCollapsed( Eigen::MatrixBase<Derived> const& points )
    {
        if ( points.rows() != Dim ) throw std::invalid_argument( "simplex point dimension mismatch" );
        ::Feel::detail::EigenMatrixBuffer<T, Dim, Derived::ColsAtCompileTime> result( Dim, points.cols() );
        auto mapped = result.matrix();
        for ( Eigen::Index q = 0; q < points.cols(); ++q )
        {
            if constexpr ( Dim == 2 )
            {
                T const denominator = T( 1 ) - points( 1, q );
                mapped( 0, q ) = std::abs( denominator ) <= boundaryTolerance ? T( -1 ) :
                    T( 2 ) * ( T( 1 ) + points( 0, q ) ) / denominator - T( 1 );
            }
            else
            {
                T const first = points( 1, q ) + points( 2, q );
                T const second = T( 1 ) - points( 2, q );
                mapped( 0, q ) = std::abs( first ) <= boundaryTolerance ? T( 1 ) :
                    -T( 2 ) * ( T( 1 ) + points( 0, q ) ) / first - T( 1 );
                mapped( 1, q ) = std::abs( second ) <= boundaryTolerance ? T( -1 ) :
                    T( 2 ) * ( T( 1 ) + points( 1, q ) ) / second - T( 1 );
            }
            mapped( Dim - 1, q ) = points( Dim - 1, q );
        }
        return result;
    }

    /** @brief Map collapsed coordinates to the Cartesian reference simplex. */
    template<typename Derived>
    [[nodiscard]] static auto toCartesian( Eigen::MatrixBase<Derived> const& points )
    {
        if ( points.rows() != Dim ) throw std::invalid_argument( "simplex point dimension mismatch" );
        ::Feel::detail::EigenMatrixBuffer<T, Dim, Derived::ColsAtCompileTime> result( Dim, points.cols() );
        auto mapped = result.matrix();
        for ( Eigen::Index q = 0; q < points.cols(); ++q )
        {
            mapped( Dim - 1, q ) = points( Dim - 1, q );
            T factor = T( 1 );
            for ( int d = Dim - 2; d >= 0; --d )
            {
                factor *= ( T( 1 ) - points( d + 1, q ) ) / T( 2 );
                mapped( d, q ) = ( T( 1 ) + points( d, q ) ) * factor - T( 1 );
            }
        }
        return result;
    }
};

/** @brief Powers ((1-x)/2)^k for static or runtime maximum degree.
 * @tparam OrderSpec Maximum power or `Dynamic`.
 * @tparam T Scalar type.
 */
template<int OrderSpec, typename T = double>
class ScalingPowers : public OrderBase<OrderSpec>
{
    using super = OrderBase<OrderSpec>;
    static_assert( OrderSpec == Dynamic || OrderSpec <= std::numeric_limits<uint16_type>::max() );
public:
    /** Degree rows by independently sized sample columns. */
    template<int Samples = Dynamic>
    using TableFor = ::Feel::detail::EigenMatrixBuffer<T, OrderSpec == Dynamic ? Dynamic : OrderSpec + 1, Samples>;
    /** Scaling rows with a runtime sample count. */
    using Table = TableFor<>;
    /** @brief Construct or validate the maximum degree. */
    explicit ScalingPowers( RuntimeOrder degree = RuntimeOrder{ OrderSpec == Dynamic ? 0 : OrderSpec } )
        : super( degree ) {}
    /** @brief Evaluate the power table into a new owner. */
    template<typename Derived>
    [[nodiscard]] auto values( Eigen::MatrixBase<Derived> const& points ) const
    {
        TableFor<::Feel::detail::pointSampleExtent<Derived>()> result( this->order() + 1, ::Feel::detail::pointCount( points ), this->order() + 1, Dynamic );
        valuesInto( result, points );
        return result;
    }
    /** @brief Evaluate the power table in reusable storage, preserving alias safety. */
    template<::Feel::detail::EigenBuffer Buffer, typename Derived>
        requires std::same_as<typename Buffer::value_type, T>
    void valuesInto( Buffer& output, Eigen::MatrixBase<Derived> const& points ) const
    { tabulateInto( output, points ); }

    /** @brief Fill a correctly sized Eigen or mapped legacy destination directly. */
    template<typename Output, typename Derived>
        requires( ( Output::Flags & Eigen::DirectAccessBit ) != 0 )
    void valuesInto( Eigen::MatrixBase<Output>& output, Eigen::MatrixBase<Derived> const& points ) const
    {
        ::Feel::detail::EigenTabulationDestination<Output> destination( output );
        tabulateInto( destination, points );
    }

private:
    /** @brief Share the full-write scaling recurrence across owning and borrowed outputs. */
    template<typename Buffer, typename Derived>
    void tabulateInto( Buffer& output, Eigen::MatrixBase<Derived> const& points ) const
    {
        ::Feel::detail::tabulatePointRows( output, points, this->order() + 1,
            [&]( auto& target, auto const& x )
            {
                target.row( 0 ).setOnes();
                if ( this->order() == 0 ) return;
                target.row( 1 ).array() = ( T( 1 ) - x ) / T( 2 );
                for ( int n = 2; n <= this->order(); ++n )
                    target.row( n ).array() = target.row( n - 1 ).array() * target.row( 1 ).array();
            } );
    }
};
} // namespace Feel
#endif // FEELPP_FEELPOLY_EXPANSIONKERNELS_HPP
