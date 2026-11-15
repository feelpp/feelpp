/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file pointtabulation.hpp
 * @brief Shared orientation, aliasing and reusable storage for univariate batches.
 */
#ifndef FEELPP_FEELPOLY_POINTTABULATION_HPP
#define FEELPP_FEELPOLY_POINTTABULATION_HPP 1
#include <feel/feelpoly/eigenextent.hpp>
#include <algorithm>
#include <cstdint>

namespace Feel::detail
{
/** @brief Adapt a caller-owned Eigen destination without resizing its storage.
 * @tparam Derived Mutable Eigen matrix or view, including a mapped legacy owner.
 * @note This adapter borrows storage only for the duration of one tabulation.
 */
template<typename Derived>
class EigenTabulationDestination
{
public:
    /** Scalar type of the borrowed output. */
    using value_type = typename Derived::Scalar;
    /** @brief Borrow a writable destination whose owner outlives the operation. */
    explicit EigenTabulationDestination( Eigen::MatrixBase<Derived>& output ) : M_output( output.derived() ) {}
    /** @brief Require the caller's destination to have the complete output shape. */
    void resizeForOverwrite( Eigen::Index rows, Eigen::Index cols )
    {
        (void)matrixElementCount( rows, cols, sizeof( value_type ) );
        if ( M_output.rows() != rows || M_output.cols() != cols )
            throw std::invalid_argument( "polynomial tabulation destination dimensions disagree" );
    }
    /** @brief Borrow the mutable matrix without allocating another owner. */
    [[nodiscard]] Derived& matrix() noexcept { return M_output; }
private:
    /** Caller-owned output; never resized or retained after the call. */
    Derived& M_output;
};

/** @brief Preserve a statically known vector length independently of polynomial degree. */
template<typename Points>
[[nodiscard]] consteval int pointSampleExtent()
{
    if constexpr ( Points::RowsAtCompileTime == 1 ) return Points::ColsAtCompileTime;
    else if constexpr ( Points::ColsAtCompileTime == 1 ) return Points::RowsAtCompileTime;
    else if constexpr ( Points::RowsAtCompileTime == Dynamic ) return Points::ColsAtCompileTime;
    else if constexpr ( Points::ColsAtCompileTime == Dynamic ) return Points::RowsAtCompileTime;
    else return Dynamic;
}

/** @brief Validate vector orientation before deriving the sample count or allocating output. */
template<typename Derived>
[[nodiscard]] Eigen::Index pointCount( Eigen::MatrixBase<Derived> const& points )
{
    if ( points.rows() == 1 ) return points.cols();
    if ( points.cols() == 1 ) return points.rows();
    throw std::invalid_argument( "sample coordinates must be a row or column vector" );
}

/** @brief Conservatively check overlap with a direct-access, possibly strided output. */
template<typename Points, typename Output>
[[nodiscard]] bool pointStorageOverlaps( Points const& points, Output const& output )
{
    if ( points.size() == 0 || output.size() == 0 ) return false;
    auto const* first = points.data();
    auto const stride = points.rows() == 1 ? points.colStride() : points.rowStride();
    auto const* last = first + ( points.size() - 1 ) * stride;
    auto const begin = std::min( reinterpret_cast<std::uintptr_t>( first ), reinterpret_cast<std::uintptr_t>( last ) );
    auto const end = std::max( reinterpret_cast<std::uintptr_t>( first ), reinterpret_cast<std::uintptr_t>( last ) ) + sizeof( typename Points::Scalar );
    auto const rowOffset = ( output.rows() - 1 ) * output.rowStride();
    auto const colOffset = ( output.cols() - 1 ) * output.colStride();
    auto const* outFirst = output.data() + std::min<Eigen::Index>( 0, rowOffset ) + std::min<Eigen::Index>( 0, colOffset );
    auto const* outLast = output.data() + std::max<Eigen::Index>( 0, rowOffset ) + std::max<Eigen::Index>( 0, colOffset );
    auto const outBegin = reinterpret_cast<std::uintptr_t>( outFirst );
    auto const outEnd = reinterpret_cast<std::uintptr_t>( outLast ) + sizeof( typename Output::Scalar );
    return begin < outEnd && outBegin < end;
}

/** @brief Apply a callable to a point row, preserving static vector orientation. */
template<typename Points, typename Operation>
void withPointRow( Points const& points, Operation&& operation )
{
    if constexpr ( Points::RowsAtCompileTime == 1 ) operation( points );
    else if constexpr ( Points::ColsAtCompileTime == 1 ) operation( points.transpose() );
    else if ( points.rows() == 1 ) operation( points.row( 0 ) );
    else operation( points.col( 0 ).transpose() );
}

/** @brief Validate and prepare one batch, snapshotting only expression/aliased points.
 * @param output Owner or destination whose complete shape is checked before writing.
 * @param points Row or column vector, possibly strided or aliased with output.
 * @param rows Number of mathematical output rows.
 * @param operation Callable receiving a mutable output view and a coordinate array row.
 */
template<typename Buffer, typename Derived, typename Operation>
void tabulatePointRows( Buffer& output, Eigen::MatrixBase<Derived> const& points,
                        Eigen::Index rows, Operation&& operation )
{
    using T = typename Buffer::value_type;
    auto const samples = pointCount( points );
    auto evaluate = [&]( auto const& source )
    {
        output.resizeForOverwrite( rows, source.size() );
        decltype(auto) target = output.matrix();
        withPointRow( source, [&]( auto const& row )
        { operation( target, row.template cast<T>().array() ); } );
    };
    bool snapshot = true;
    if constexpr ( Derived::Flags & Eigen::DirectAccessBit )
        snapshot = pointStorageOverlaps( points.derived(), output.matrix() );
    if ( snapshot )
    {
        Eigen::Matrix<T, 1, Dynamic> saved( 1, samples );
        withPointRow( points.derived(), [&]( auto const& row )
        { saved = row.template cast<T>(); } );
        evaluate( saved );
    }
    else evaluate( points.derived() );
}
} // namespace Feel::detail
#endif // FEELPP_FEELPOLY_POINTTABULATION_HPP
