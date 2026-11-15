/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file referenceoperators.hpp
 * @brief Immutable reference differentiation data and bounded requested-order reuse.
 */
#ifndef FEELPP_FEELPOLY_REFERENCEOPERATORS_HPP
#define FEELPP_FEELPOLY_REFERENCEOPERATORS_HPP 1

#include <array>
#include <feel/feelalg/eigenmap.hpp>
#include <feel/feelpoly/eigenextent.hpp>
#include <feel/feelpoly/polynomialsolve.hpp>
#include <memory>
#include <mutex>
#include <string_view>

namespace Feel::detail
{
/** @brief Borrow a matrix or native buffer without retaining a temporary expression. */
template <typename Matrix>
decltype( auto ) referenceOperatorView( Matrix& matrix )
{
    if constexpr ( EigenMatrix<Matrix> )
        return ( matrix );
    else if constexpr ( requires { matrix.matrix(); } )
        return matrix.matrix();
    else
        return eigenMap( matrix );
}

/** @brief Own every reference direction for one effective polynomial degree.
 * @tparam Matrix Owning compatibility representation, including scalar/layout policy.
 * @tparam Dimension Fixed reference dimension.
 */
template <typename Matrix, int Dimension>
struct ReferenceDerivativeOperators
{
    /** Effective degree, independent of the descriptor's seed. */
    uint16_type degree;
    /** Fixed direction collection, immutable after publication. */
    std::array<Matrix, Dimension> matrices;

    /** @brief Borrow a checked direction while its owner is retained. */
    [[nodiscard]] Matrix const& matrix( uint16_type direction ) const& { return matrices.at( direction ); }
    /** @brief Reject a borrowed matrix from a temporary owner. */
    Matrix const& matrix( uint16_type ) const&& = delete;

    /** @brief Count coefficient bytes for bounded retention, excluding object headers. */
    [[nodiscard]] std::size_t coefficientBytes() const
    {
        std::size_t bytes = 0;
        for ( auto const& matrix : matrices )
        {
            decltype( auto ) values = referenceOperatorView( matrix );
            auto const count = matrixElementCount( values.rows(), values.cols(), sizeof( typename Matrix::value_type ) );
            auto const added = checkedProduct( static_cast<std::size_t>( count ), sizeof( typename Matrix::value_type ) );
            if ( added > std::numeric_limits<std::size_t>::max() - bytes )
                throw std::length_error( "reference operator byte count overflows" );
            bytes += added;
        }
        return bytes;
    }
};

/** @brief Retain a bounded set of immutable requested-order configurations.
 * @tparam Family Type encoding family, dimension, scalar, normalization and storage.
 * @tparam Data Immutable data with a coefficientBytes() query.
 * @tparam Capacity Maximum indexed effective degrees.
 * @tparam RetainedBytes Maximum strongly retained coefficient bytes per family.
 * @note Live callers retain evicted data. Oversized data use weak reuse only;
 * static reference-returning facades retain their own permanent handle.
 */
template <typename Family, typename Data, std::size_t Capacity = 16,
          std::size_t RetainedBytes = 32 * 1024 * 1024>
class ReferenceOperatorCache
{
    static_assert( Capacity > 0 );
    /** One weak index and optional budgeted keep-alive owner. */
    struct Entry
    {
        /** Effective degree in this slot. */
        uint16_type degree = 0;
        /** Discover caller-owned data after keep-alive eviction. */
        std::weak_ptr<Data const> weak;
        /** Retain repeated value-returning compatibility requests within budget. */
        std::shared_ptr<Data const> retained;
        /** Coefficient bytes charged to the retention budget. */
        std::size_t bytes = 0;
    };

  public:
    /** @brief Find or construct every direction, publishing only complete checked data. */
    template <typename Build>
    [[nodiscard]] static std::shared_ptr<Data const> get( RuntimeOrder degree, Build&& build )
    {
        // Serialize first construction; exceptions preserve every existing entry.
        std::lock_guard lock( M_mutex );
        std::size_t slot = M_next;
        for ( std::size_t i = 0; i < Capacity; ++i )
        {
            if ( auto owner = M_entries[i].weak.lock() )
            {
                if ( M_entries[i].degree == degree.value ) return owner;
            }
            else
                slot = i;
        }
        auto owner = std::make_shared<Data const>( std::forward<Build>( build )() );
        if ( owner->degree != degree.value )
            throw std::invalid_argument( "reference operator degree disagrees with its cache key" );
        auto const bytes = owner->coefficientBytes();

        // Drop only cache ownership, never a caller's live matrices or views.
        release( M_entries[slot] );
        if ( bytes <= RetainedBytes )
            for ( std::size_t i = 0; i < Capacity && M_bytes > RetainedBytes - bytes; ++i )
                release( M_entries[( M_next + i ) % Capacity] );
        auto& entry = M_entries[slot];
        entry.degree = degree.value;
        entry.weak = owner;
        if ( bytes <= RetainedBytes )
        {
            entry.retained = owner;
            entry.bytes = bytes;
            M_bytes += bytes;
        }
        M_next = ( slot + 1 ) % Capacity;
        return owner;
    }

  private:
    /** @brief Drop a keep-alive without invalidating its indexed weak handle. */
    static void release( Entry& entry )
    {
        M_bytes -= entry.bytes;
        entry.retained.reset();
        entry.bytes = 0;
    }
    /** Fixed index avoids an unbounded runtime map. */
    inline static std::array<Entry, Capacity> M_entries;
    /** Next replacement/retention scan position. */
    inline static std::size_t M_next = 0;
    /** Strongly retained coefficient bytes. */
    inline static std::size_t M_bytes = 0;
    /** Protect lookup, construction and budget accounting. */
    inline static std::mutex M_mutex;
};

/** @brief Solve every D_i V = partial_i V with one interpolation factorization.
 * @note Eigen performs the numerical work; compatibility outputs are materialized once.
 */
template <int Dimension, typename Matrix, typename Values, typename Gradients>
[[nodiscard]] auto interpolateReferenceOperators( RuntimeOrder degree, Values const& values,
                                                  Gradients const& gradients,
                                                  typename Matrix::value_type cleanupTolerance,
                                                  std::string_view context )
{
    using T = typename Matrix::value_type;
    ReferenceDerivativeOperators<Matrix, Dimension> result{ degree.value, {} };
    decltype( auto ) modal = referenceOperatorView( values );
    PolynomialLuSolver<T> solver( modal, context );
    for ( int direction = 0; direction < Dimension; ++direction )
    {
        decltype( auto ) derivative = referenceOperatorView( gradients[direction] );
        auto solution = solver.solveRight( derivative );
        Matrix matrix( solution.rows(), solution.cols() );
        decltype( auto ) target = referenceOperatorView( matrix );
        auto clean = [cleanupTolerance]( T value )
        {
            using std::abs;
            return abs( value ) < cleanupTolerance ? T( 0 ) : value;
        };
        target = solution.unaryExpr( clean );
        result.matrices[direction] = std::move( matrix );
    }
    return result;
}
} // namespace Feel::detail
#endif // FEELPP_FEELPOLY_REFERENCEOPERATORS_HPP
