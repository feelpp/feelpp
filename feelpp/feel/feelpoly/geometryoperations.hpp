/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file geometryoperations.hpp
 * @brief Representation-independent geometry conversion and restriction operations.
 */
#if !defined( FEELPP_FEELPOLY_GEOMETRYOPERATIONS_HPP )
#define FEELPP_FEELPOLY_GEOMETRYOPERATIONS_HPP 1

#include <feel/feelpoly/bernsteinbezier.hpp>
#include <feel/feelpoly/geometrynodallagrange.hpp>

#include <stdexcept>
#include <utility>

namespace Feel::geometry
{
//!
//! \brief Evaluate a represented geometric polynomial.
//! \tparam Representation geometry representation type
//! \param source coefficient table bound to its representation
//! \param points reference points stored by column
//! \return component-by-point physical values
template <typename Representation>
[[nodiscard]] auto evaluate(
    CoefficientTable<Representation> const& source,
    typename Representation::points_type const& points )
{
    typename CoefficientTable<Representation>::matrix_type result =
        source.values() * source.representation().basisValues( points );
    return result;
}

//!
//! \brief Convert nodal interpolation values to Bezier control points.
//! \tparam NodalRepresentation nodal-Lagrange representation type
//! \tparam BezierRepresentation Bernstein/Bezier representation type
//! \param source nodal coefficient table
//! \param target target Bezier representation
//! \return Bezier coefficient table representing the same polynomial map
template <typename NodalRepresentation, typename BezierRepresentation>
    requires( NodalRepresentation::coefficientSemantics() ==
                  CoefficientSemantics::InterpolationPoint &&
              BezierRepresentation::coefficientSemantics() ==
                  CoefficientSemantics::ControlPoint )
[[nodiscard]] auto convert( CoefficientTable<NodalRepresentation> const& source,
                            BezierRepresentation target )
{
    if ( source.representation().order() != target.order() ||
         source.size() != target.coefficientCount() )
        throw std::invalid_argument( "Incompatible geometry representations" );
    auto values = source.representation().toBezier( source.values() );
    return CoefficientTable<BezierRepresentation>( std::move( target ), std::move( values ) );
}

//!
//! \brief Convert Bezier control points to nodal interpolation values.
//! \tparam BezierRepresentation Bernstein/Bezier representation type
//! \tparam NodalRepresentation nodal-Lagrange representation type
//! \param source Bezier coefficient table
//! \param target target nodal representation
//! \return nodal coefficient table representing the same polynomial map
template <typename BezierRepresentation, typename NodalRepresentation>
    requires( BezierRepresentation::coefficientSemantics() ==
                  CoefficientSemantics::ControlPoint &&
              NodalRepresentation::coefficientSemantics() ==
                  CoefficientSemantics::InterpolationPoint )
[[nodiscard]] auto convert( CoefficientTable<BezierRepresentation> const& source,
                            NodalRepresentation target )
{
    if ( source.representation().order() != target.order() ||
         source.size() != target.coefficientCount() )
        throw std::invalid_argument( "Incompatible geometry representations" );
    auto values = target.fromBezier( source.values() );
    return CoefficientTable<NodalRepresentation>( std::move( target ), std::move( values ) );
}

//!
//! \brief Restrict a Bezier coefficient table exactly to a facet.
//! \tparam Representation Bernstein/Bezier representation type
//! \tparam FacetDescriptor simplex or hypercube facet descriptor
//! \param source source coefficient table
//! \param facet facet descriptor
//! \return component-by-facet-control-point matrix
template <typename Representation, typename FacetDescriptor>
[[nodiscard]] auto restrictToFacet( CoefficientTable<Representation> const& source,
                                    FacetDescriptor facet )
{
    auto const indices = source.representation().facetCoefficientIndices( facet );
    typename CoefficientTable<Representation>::matrix_type result(
        source.values().rows(), indices.size() );
    for ( size_type i = 0; i < indices.size(); ++i )
        result.col( i ) = source.values().col( indices[i] );
    return result;
}

//!
//! \brief Subdivide a Bezier curve by the de Casteljau algorithm.
//! \tparam Order static polynomial order or `Dynamic`
//! \tparam Entity one-dimensional reference-entity family
//! \tparam T scalar type
//! \param source Bezier curve coefficient table
//! \param parameter subdivision parameter in `[0,1]`
//! \return left and right control-point tables
template <int Order, template <int, int, int> class Entity, typename T>
[[nodiscard]] auto subdivide(
    CoefficientTable<BernsteinBezier<1, Order, Entity, T>> const& source, T parameter )
{
    if ( parameter < T( 0 ) || parameter > T( 1 ) )
        throw std::out_of_range( "Bezier subdivision parameter must lie in [0,1]" );
    using representation_type = BernsteinBezier<1, Order, Entity, T>;
    using coefficient_type = CoefficientTable<representation_type>;
    using matrix_type = typename coefficient_type::matrix_type;
    Eigen::Index const count = source.values().cols();
    matrix_type work = source.values();
    matrix_type left( source.values().rows(), count );
    matrix_type right( source.values().rows(), count );
    left.col( 0 ) = work.col( 0 );
    right.col( count - 1 ) = work.col( count - 1 );
    for ( Eigen::Index level = 1; level < count; ++level )
    {
        for ( Eigen::Index i = 0; i < count - level; ++i )
            work.col( i ) = ( T( 1 ) - parameter ) * work.col( i ) +
                            parameter * work.col( i + 1 );
        left.col( level ) = work.col( 0 );
        right.col( count - level - 1 ) = work.col( count - level - 1 );
    }
    return std::pair{ coefficient_type( source.representation(), std::move( left ) ),
                      coefficient_type( source.representation(), std::move( right ) ) };
}

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_GEOMETRYOPERATIONS_HPP
