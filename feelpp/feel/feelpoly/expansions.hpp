/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2005-2006 EPFL
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file expansions.hpp
 * @brief Compatibility boundaries for shared Eigen simplex expansion kernels.
 */
#ifndef FEELPP_FEELPOLY_EXPANSIONS_HPP
#define FEELPP_FEELPOLY_EXPANSIONS_HPP 1
#include <feel/feelalg/glas.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feelpoly/jacobi.hpp>
#include <feel/feelpoly/expansionkernels.hpp>
#include <type_traits>

namespace Feel
{
/** Simplex reference shapes accepted by expansion adapters. */
enum TheShape { LINE = 1, TRIANGLE = 2, TETRAHEDRON = 3 };
/** @brief Reference dimension without an MPL dispatch hierarchy. */
template<TheShape Shape>
struct DimFromShape
{
    /** Fixed spatial dimension. */
    static constexpr uint16_type value = static_cast<uint16_type>( Shape );
};

namespace details
{
/** @brief Scalar/batch legacy owner forwarding coordinate maps to Eigen.
 * @tparam Shape Triangle or tetrahedron.
 * @tparam T Scalar type.
 * @tparam Collapsed Map Cartesian points to collapsed coordinates when true.
 * @tparam Batch Accept a point matrix when true, a single node otherwise.
 */
template<TheShape Shape, typename T, bool Collapsed, bool Batch>
class CoordinateMap
{
    /** Fixed spatial dimension. */
    static constexpr int dim = DimFromShape<Shape>::value;
public:
    /** Coordinate scalar type. */
    using value_type = T;
    /** Established scalar-point representation. */
    using node_type = typename node<T>::type;
    /** Established point-batch representation. */
    using matrix_node_type = typename matrix_node<T>::type;
    /** Legacy result representation. */
    using storage_type = std::conditional_t<Batch, matrix_node_type, node_type>;
    /** @brief Initialize an empty batch or a zero scalar point. */
    CoordinateMap() : M_coordinates( makeStorage() ) {}
    /** @brief Map initial coordinates through the ordinary call path. */
    explicit CoordinateMap( storage_type const& points ) : CoordinateMap() { (*this)( points ); }
    /** @brief Map points and retain the materialized legacy result. */
    storage_type const& operator()( storage_type const& points )
    {
        auto transform = [&]( auto const& coordinates )
        {
            auto mapped = [&]
            {
                if constexpr ( Collapsed ) return SimplexCoordinates<dim, T>::toCollapsed( coordinates );
                else return SimplexCoordinates<dim, T>::toCartesian( coordinates );
            }();
            if constexpr ( Batch )
            {
                M_coordinates.resize( dim, mapped.cols() );
                eigenMap( M_coordinates ) = mapped.matrix();
            }
            else eigenMap( M_coordinates ) = mapped.matrix().col( 0 );
        };
        if constexpr ( Batch )
        {
            // A proven single column allows bounded inline coordinates without
            // changing the legacy result type or adding another mapping kernel.
            if ( points.size2() == 1 ) transform( eigenMap<dim, 1>( points ) );
            else transform( eigenMap<dim, Dynamic>( points ) );
        }
        else transform( eigenMap( points ) );
        return M_coordinates;
    }
    /** @brief Result retained by the most recent evaluation. */
    storage_type const& operator()() const noexcept { return M_coordinates; }
private:
    /** @brief Initialize storage without duplicating coordinate maps. */
    static storage_type makeStorage()
    {
        if constexpr ( Batch ) return storage_type( dim, 0 );
        else return storage_type( dim );
    }
    /** Sole legacy result owner. */
    storage_type M_coordinates;
};
/** Scalar Cartesian-coordinate adapter. */
template<TheShape Shape, typename T = double>
using CartesianCoordinates = CoordinateMap<Shape, T, false, false>;
/** Scalar collapsed-coordinate adapter. */
template<TheShape Shape, typename T = double>
using CollapsedCoordinates = CoordinateMap<Shape, T, true, false>;
/** Batch collapsed-coordinate adapter. */
template<TheShape Shape, typename T = double>
using CollapsedPointCoordinates = CoordinateMap<Shape, T, true, true>;
/** Compatibility spelling for quadrature clients. */
template<TheShape Shape, typename T = double> using xi = CartesianCoordinates<Shape, T>;
/** Compatibility spelling for scalar clients. */
template<TheShape Shape, typename T = double> using eta = CollapsedCoordinates<Shape, T>;
/** Compatibility spelling for modal clients. */
template<TheShape Shape, typename T = double> using etas = CollapsedPointCoordinates<Shape, T>;

/** @brief Scalar scaled Jacobi factors used by legacy simplex expansions. */
template<typename T>
class PsiTilde
{
public:
    /** Scalar type of values. */
    using value_type = T;
    /** Runtime Jacobi facade using the native recurrence. */
    using P = Feel::dyna::Jacobi<T>;
    /** @brief Construct a constant factor. */
    PsiTilde() : p( 0 ), M_b( 0 ) {}
    /** @brief Construct an ordinary one-dimensional factor. */
    explicit PsiTilde( int i ) : p( i ), M_b( i ) {}
    /** @brief Construct a triangular factor. */
    PsiTilde( int i, int j ) : p( j, T( 2 * i + 1 ), T( 0 ) ), M_b( i ) {}
    /** @brief Construct a tetrahedral factor. */
    PsiTilde( int i, int j, int k ) : p( k, T( 2 * ( i + j + 1 ) ), T( 0 ) ), M_b( i + j ) {}
    /** @brief Ordinary Jacobi factor. */
    [[nodiscard]] T a( T const& x ) const { return p( x ); }
    /** @brief Scaled Jacobi factor. */
    [[nodiscard]] T b( T const& x ) const { return math::pow( T( .5 ) * ( T( 1 ) - x ), M_b ) * p( x ); }
    /** @brief Tetrahedral spelling of the same scaled factor. */
    [[nodiscard]] T c( T const& x ) const { return b( x ); }
    /** Jacobi factor retained by the legacy interface. */
    P p;
    /** Scaling exponent. */
    T M_b;
};
/** Compatibility spelling for scalar-factor clients. */
template<typename T> using psitilde = PsiTilde<T>;

/** @brief Static scaling adapter using the shared Eigen recurrence. */
template<uint16_type N, typename T = double>
class ScalingFactors
{
public:
    /** Scalar type of factors. */
    using value_type = T;
    /** Established legacy result representation. */
    using matrix_type = ublas::matrix<T>;
    /** Established point-batch representation. */
    using matrix_node_type = typename matrix_node<T>::type;
    /** @brief Evaluate the static recurrence and materialize the boundary once. */
    explicit ScalingFactors( ublas::vector<T> const& points )
        : M_s( ::Feel::detail::makeLegacyPolynomialTable<T>( N + 1, points.size() ) )
    {
        auto output = eigenMap( M_s );
        ScalingPowers<N, T>{}.valuesInto( output, eigenMap( points ) );
    }
    /** @brief Scaling rows retained by the compatibility owner. */
    matrix_type const& operator()() const noexcept { return M_s; }
    /** Legacy result evaluated by ScalingPowers. */
    matrix_type M_s;
};
/** Compatibility spelling for static modal clients. */
template<uint16_type N, typename T = double> using scalings = ScalingFactors<N, T>;
/** @brief Runtime scaling boundary using the same native recurrence. */
template<typename T>
ublas::matrix<T> scalingsRuntime( uint16_type order, ublas::vector<T> const& points )
{
    auto result = ::Feel::detail::makeLegacyPolynomialTable<T>( order + 1, points.size() );
    auto output = eigenMap( result );
    ScalingPowers<Dynamic, T>{ RuntimeOrder{ order } }.valuesInto( output, eigenMap( points ) );
    return result;
}
} // namespace details
} // namespace Feel
#endif // FEELPP_FEELPOLY_EXPANSIONS_HPP
