/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file polynomialjet.hpp
 * @brief Compact value/gradient/symmetric-Hessian algebra for polynomial recurrences.
 */
#ifndef FEELPP_FEELPOLY_POLYNOMIALJET_HPP
#define FEELPP_FEELPOLY_POLYNOMIALJET_HPP 1
#include <array>
#include <cstddef>

namespace Feel::detail
{
/** @brief Truncated derivative algebra with actual, unscaled partial derivatives.
 * @tparam Dim Fixed reference dimension.
 * @tparam Order Highest derivative order, zero through two.
 * @tparam T Scalar type.
 *
 * Slots are value, coordinate-ordered first partials, then upper-triangular
 * second partials (00,01,...,11,12,...). Mixed derivatives have one slot.
 * Multiplication applies the product rule, including the factor two for a
 * repeated direction. There is no coordinate division or dynamic allocation.
 */
template<int Dim, unsigned Order, typename T>
struct PolynomialJet
{
    static_assert( Dim >= 1 && Dim <= 3 && Order <= 2 );
    /** Number of retained scalar derivatives, including the value. */
    static constexpr int componentCount = 1 + ( Order >= 1 ? Dim : 0 ) +
        ( Order >= 2 ? Dim * ( Dim + 1 ) / 2 : 0 );
    /** Actual derivatives in the documented compact order. */
    std::array<T, componentCount> data{};

    /** @brief Construct a constant with zero derivatives. */
    explicit PolynomialJet( T value = T( 0 ) ) { data[0] = value; }
    /** @brief Seed a reference coordinate and its unit first partial. */
    [[nodiscard]] static PolynomialJet coordinate( T value, int direction )
    {
        PolynomialJet result( value );
        if constexpr ( Order >= 1 ) result.data[1 + direction] = T( 1 );
        return result;
    }
    /** @brief Return the value. */
    [[nodiscard]] T value() const { return data[0]; }
    /** @brief Return a first partial for a validated coordinate direction. */
    [[nodiscard]] T gradient( int direction ) const requires( Order >= 1 )
    { return data[1 + direction]; }
    /** @brief Map a validated symmetric direction pair to its derivative slot. */
    [[nodiscard]] static constexpr int hessianSlot( int i, int j )
    {
        if ( i > j ) { int const first = i; i = j; j = first; }
        return 1 + Dim + i * Dim - i * ( i + 1 ) / 2 + j;
    }
    /** @brief Return an actual second partial, without factorial scaling. */
    [[nodiscard]] T hessian( int i, int j ) const requires( Order >= 2 )
    { return data[hessianSlot( i, j )]; }
    /** @brief Add two jets componentwise. */
    friend PolynomialJet operator+( PolynomialJet const& a, PolynomialJet const& b )
    {
        PolynomialJet result;
        for ( int c = 0; c < componentCount; ++c ) result.data[c] = a.data[c] + b.data[c];
        return result;
    }
    /** @brief Subtract two jets componentwise. */
    friend PolynomialJet operator-( PolynomialJet const& a, PolynomialJet const& b )
    {
        PolynomialJet result;
        for ( int c = 0; c < componentCount; ++c ) result.data[c] = a.data[c] - b.data[c];
        return result;
    }
    /** @brief Scale every retained derivative by a constant. */
    friend PolynomialJet operator*( T scale, PolynomialJet const& a )
    {
        PolynomialJet result;
        for ( int c = 0; c < componentCount; ++c ) result.data[c] = scale * a.data[c];
        return result;
    }
    /** @brief Multiply through the same product rule for every recurrence factor. */
    friend PolynomialJet operator*( PolynomialJet const& a, PolynomialJet const& b )
    {
        PolynomialJet result( a.value() * b.value() );
        if constexpr ( Order >= 1 )
            for ( int i = 0; i < Dim; ++i )
                result.data[1 + i] = a.gradient( i ) * b.value() + a.value() * b.gradient( i );
        if constexpr ( Order >= 2 )
            for ( int i = 0; i < Dim; ++i )
                for ( int j = i; j < Dim; ++j )
                    result.data[hessianSlot( i, j )] = a.hessian( i, j ) * b.value() +
                        a.gradient( i ) * b.gradient( j ) + a.gradient( j ) * b.gradient( i ) +
                        a.value() * b.hessian( i, j );
        return result;
    }
};
} // namespace Feel::detail
#endif // FEELPP_FEELPOLY_POLYNOMIALJET_HPP
