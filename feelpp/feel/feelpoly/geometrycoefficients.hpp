/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file geometrycoefficients.hpp
 * @brief Mathematical coefficient semantics and typed geometry coefficient tables.
 */
#if !defined( FEELPP_FEELPOLY_GEOMETRYCOEFFICIENTS_HPP )
#define FEELPP_FEELPOLY_GEOMETRYCOEFFICIENTS_HPP 1

#include <feel/feelcore/feeltypes.hpp>

#include <Eigen/Core>

#include <stdexcept>
#include <utility>

namespace Feel::geometry
{
//!
//! \enum CoefficientSemantics
//! \brief Mathematical meaning of geometry coefficients.
enum class CoefficientSemantics
{
    /** Coefficients are physical values at interpolation nodes. */
    InterpolationPoint,
    /** Coefficients are Bernstein/Bezier control points. */
    ControlPoint
};

//! \brief Identify a simplex facet by its vanishing barycentric coordinate.
struct SimplexFacet
{
    /** Barycentric coordinate equal to zero on the facet. */
    uint16_type omittedBarycentricCoordinate;
};

//! \brief Identify a hypercube facet by a coordinate axis and endpoint.
struct HypercubeFacet
{
    /** Reference coordinate axis fixed on the facet. */
    uint16_type axis;
    /** Select `+1` when true and `-1` when false. */
    bool upper;
};

//! \class CoefficientTable
//! \brief Coefficient table bound to its mathematical representation.
//! \tparam Representation geometry representation type
template <typename Representation>
class CoefficientTable
{
  public:
    /** Representation type defining coefficient semantics. */
    using representation_type = Representation;
    /** Scalar value type. */
    using value_type = typename representation_type::value_type;
    /** Dynamic component-by-coefficient matrix. */
    using matrix_type = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;

    /**
     * @brief Construct a checked coefficient table.
     * @param representation representation instance defining the cardinality
     * @param values component-by-coefficient values
     */
    CoefficientTable( representation_type representation, matrix_type values )
        : M_representation( std::move( representation ) ), M_values( std::move( values ) )
    {
        if ( M_values.cols() != static_cast<Eigen::Index>( M_representation.coefficientCount() ) )
            throw std::invalid_argument( "Geometry coefficient table has an invalid column count" );
    }

    /** @brief Return the bound representation. */
    [[nodiscard]] representation_type const& representation() const noexcept
    {
        return M_representation;
    }

    /** @brief Return the component-by-coefficient matrix. */
    [[nodiscard]] matrix_type const& values() const noexcept { return M_values; }
    /** @brief Return mutable component-by-coefficient values. */
    [[nodiscard]] matrix_type& values() noexcept { return M_values; }
    /** @brief Return the number of geometry coefficients. */
    [[nodiscard]] size_type size() const noexcept { return M_values.cols(); }

  private:
    /** Representation defining coefficient semantics and cardinality. */
    representation_type M_representation;
    /** Component-by-coefficient values. */
    matrix_type M_values;
};

} // namespace Feel::geometry

#endif // FEELPP_FEELPOLY_GEOMETRYCOEFFICIENTS_HPP
