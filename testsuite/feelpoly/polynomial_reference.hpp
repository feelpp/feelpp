/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#ifndef FEELPP_TEST_POLYNOMIAL_REFERENCE_HPP
#define FEELPP_TEST_POLYNOMIAL_REFERENCE_HPP

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <array>
#include <stdexcept>
#include <utility>
#include <vector>

namespace polynomial_test
{
/** @brief Decimal arithmetic retaining 100 significant digits in the test oracle. */
using precise_type = boost::multiprecision::cpp_dec_float_100;

/** @brief Independent value, gradient and Hessian of a reference polynomial.
 * @tparam Dim Reference-space dimension.
 *
 * Uses finite binomial expansion, not the production
 * three-term recurrence or inverse collapsed-coordinate transformation.
 * Derivatives use the product rule in 100-digit arithmetic. Conventions:
 * biunit simplex, degree-ordered Dubiner modes, no factorial-scaled derivatives.
 */
template<int Dim>
struct Reference
{
    /** Polynomial value at the reference point. */
    precise_type value = 0;
    /** First derivatives in reference-coordinate order. */
    std::array<precise_type, Dim> gradient{};
    /** Actual second derivatives, indexed by `i * Dim + j`. */
    std::array<precise_type, Dim * Dim> hessian{};

    /** @brief Construct the zero polynomial and its zero derivatives. */
    Reference() = default;
    /** @brief Construct a constant polynomial.
     * @param v Constant value; all derivatives are zero.
     */
    Reference( precise_type const& v ) : value( v ) {}

    /** @brief Seed an independent reference coordinate for differentiation.
     * @param value Coordinate value at the evaluation point.
     * @param direction Coordinate axis, in `[0, Dim)`.
     * @return Coordinate value with a unit derivative along the selected axis.
     */
    static Reference coordinate( double value, int direction )
    {
        Reference result{ precise_type( value ) };
        result.gradient[direction] = 1;
        return result;
    }

    /** @brief Add polynomial values and their derivatives componentwise.
     * @param a First polynomial at the common reference point.
     * @param b Second polynomial at the same point.
     * @return Value, gradient and Hessian of the sum.
     */
    friend Reference operator+( Reference const& a, Reference const& b )
    {
        Reference result{ a.value + b.value };
        for ( int i = 0; i < Dim; ++i )
        {
            result.gradient[i] = a.gradient[i] + b.gradient[i];
            for ( int j = 0; j < Dim; ++j )
                result.hessian[i * Dim + j] =
                    a.hessian[i * Dim + j] + b.hessian[i * Dim + j];
        }
        return result;
    }

    /** @brief Multiply polynomials using the first- and second-order product rules.
     * @param a First polynomial at the common reference point.
     * @param b Second polynomial at the same point.
     * @return Value, gradient and Hessian of the product, without factorial scaling.
     */
    friend Reference operator*( Reference const& a, Reference const& b )
    {
        Reference result{ a.value * b.value };
        for ( int i = 0; i < Dim; ++i )
        {
            result.gradient[i] = a.gradient[i] * b.value + a.value * b.gradient[i];
            for ( int j = 0; j < Dim; ++j )
                result.hessian[i * Dim + j] =
                    a.hessian[i * Dim + j] * b.value + a.value * b.hessian[i * Dim + j] +
                    a.gradient[i] * b.gradient[j] + a.gradient[j] * b.gradient[i];
        }
        return result;
    }
};

/** @brief Lift a constant into the reference derivative algebra.
 * @tparam Dim Reference-space dimension.
 * @param value Constant polynomial value.
 * @return Value with zero gradient and Hessian.
 */
template<int Dim>
Reference<Dim> constant( precise_type const& value )
{
    return Reference<Dim>{ value };
}

/** @brief Evaluate a generalized binomial coefficient by a finite product.
 * @param n Upper argument, which may be nonintegral.
 * @param k Nonnegative integer lower argument.
 * @return Generalized binomial coefficient in high-precision arithmetic.
 */
inline precise_type binomial( precise_type const& n, unsigned k )
{
    precise_type result = 1;
    for ( unsigned i = 0; i < k; ++i )
        result *= ( n - i ) / ( i + 1 );
    return result;
}

/** @brief Raise a reference polynomial to a nonnegative integer power.
 * @tparam Dim Reference-space dimension.
 * @param x Polynomial value and derivatives.
 * @param exponent Integer power; zero produces the constant one.
 * @return Powered polynomial with product-rule derivatives.
 */
template<int Dim>
Reference<Dim> power( Reference<Dim> const& x, unsigned exponent )
{
    auto result = constant<Dim>( 1 );
    for ( unsigned i = 0; i < exponent; ++i )
        result = result * x;
    return result;
}

/** @brief Evaluate `v^n P_n^(a,b)(u/v)` without division, including `v = 0`.
 * @tparam Dim Reference-space dimension.
 * @param n Jacobi degree.
 * @param a First Jacobi parameter.
 * @param b Second Jacobi parameter.
 * @param u Numerator polynomial and its derivatives.
 * @param v Homogenizing polynomial and its derivatives.
 * @return Homogeneous Jacobi value, gradient and Hessian.
 * @see https://dlmf.nist.gov/18.5.E8 for the finite binomial identity.
 */
template<int Dim>
Reference<Dim> homogeneousJacobi( unsigned n, double a, double b,
                                 Reference<Dim> const& u, Reference<Dim> const& v )
{
    auto const minus = constant<Dim>( precise_type( "0.5" ) ) *
                       ( u + constant<Dim>( -1 ) * v );
    auto const plus = constant<Dim>( precise_type( "0.5" ) ) * ( u + v );
    Reference<Dim> result;
    for ( unsigned k = 0; k <= n; ++k )
        result = result + constant<Dim>( binomial( precise_type( n ) + a, k ) *
                                        binomial( precise_type( n ) + b, n - k ) ) *
                          power( minus, n - k ) * power( plus, k );
    return result;
}

/** @brief Evaluate a univariate Jacobi polynomial independently of its recurrence.
 * @param degree Polynomial degree.
 * @param a First Jacobi parameter.
 * @param b Second Jacobi parameter.
 * @param x Reference coordinate.
 * @return Value and first/second derivatives with respect to `x`.
 */
inline Reference<1> jacobi( unsigned degree, double a, double b, double x )
{
    return homogeneousJacobi( degree, a, b, Reference<1>::coordinate( x, 0 ),
                              constant<1>( 1 ) );
}

/** @brief Tabulate independent Dubiner polynomial references at one point.
 * @tparam Dim Reference dimension, restricted to 1, 2 or 3.
 * @param degree Maximum total degree of the expansion.
 * @param point Point in the biunit reference simplex.
 * @param normalized Whether to use the existing orthonormal scaling.
 * @return Modes ordered by total degree, with tetrahedron `k` outer and `j` inner.
 */
template<int Dim>
std::vector<Reference<Dim>> dubiner( unsigned degree, std::array<double, Dim> const& point,
                                    bool normalized )
{
    std::array<Reference<Dim>, Dim> x;
    for ( int d = 0; d < Dim; ++d )
        x[d] = Reference<Dim>::coordinate( point[d], d );
    auto const one = constant<Dim>( 1 );
    auto const half = constant<Dim>( precise_type( "0.5" ) );
    auto const negativeHalf = constant<Dim>( precise_type( "-0.5" ) );
    std::vector<Reference<Dim>> result;
    for ( unsigned total = 0; total <= degree; ++total )
        for ( unsigned k = 0; k <= ( Dim == 3 ? total : 0 ); ++k )
            for ( unsigned j = 0; j <= ( Dim >= 2 ? total - k : 0 ); ++j )
            {
                unsigned const i = total - j - k;
                Reference<Dim> mode;
                precise_type normalization = precise_type( i ) + precise_type( "0.5" );
                if constexpr ( Dim == 1 )
                    mode = homogeneousJacobi( i, 0, 0, x[0], one );
                else if constexpr ( Dim == 2 )
                {
                    mode = homogeneousJacobi( i, 0, 0, x[0] + half * ( one + x[1] ),
                                               half + negativeHalf * x[1] ) *
                           homogeneousJacobi( j, 2 * i + 1, 0, x[1], one );
                    normalization *= i + j + 1;
                }
                else
                {
                    mode = homogeneousJacobi( i, 0, 0, one + x[0] + half * ( x[1] + x[2] ),
                                               negativeHalf * ( x[1] + x[2] ) ) *
                           homogeneousJacobi( j, 2 * i + 1, 0, half + x[1] + half * x[2],
                                               half + negativeHalf * x[2] ) *
                           homogeneousJacobi( k, 2 * ( i + j + 1 ), 0, x[2], one );
                    normalization *= ( i + j + 1 ) * ( precise_type( total ) + precise_type( "1.5" ) );
                }
                if ( normalized )
                    mode = constant<Dim>( sqrt( normalization ) ) * mode;
                result.push_back( std::move( mode ) );
            }
    return result;
}

/** @brief Dispatch a runtime test dimension to a typed callable once.
 * @tparam Function Callable exposing a nullary `operator()<Dim>()`.
 * @param dimension Reference dimension, restricted to 1, 2 or 3.
 * @param function Callable receiving the dimension as a template argument.
 * @throws std::invalid_argument When the dimension is outside the supported set.
 */
template<typename Function>
void visitDimension( unsigned dimension, Function&& function )
{
    switch ( dimension )
    {
    case 1: function.template operator()<1>(); break;
    case 2: function.template operator()<2>(); break;
    case 3: function.template operator()<3>(); break;
    default: throw std::invalid_argument( "test dimension must be 1, 2, or 3" );
    }
}
} // namespace polynomial_test
#endif
