/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/** @file hcurlinterpolation.hpp @brief Reusable H(curl) interpolation kernels. */

#ifndef FEELPP_FEELPOLY_HCURLINTERPOLATION_HPP
#define FEELPP_FEELPOLY_HCURLINTERPOLATION_HPP 1

#include <type_traits>

#include <feel/feelcore/feeltypes.hpp>

namespace Feel::detail
{

/**
 * Apply edge-tangential point DoFs to an expression.
 *
 * Geometry/tangent acquisition is a policy supplied by the caller, keeping
 * this local kernel common to element and facet interpolation and suitable for
 * static or runtime cardinalities.
 *
 * @tparam Expr expression evaluated at functional support points
 * @tparam LocalVector mutable local coefficient vector
 * @tparam Tangent tangent vector storage
 * @tparam TangentProvider callable providing the physical edge tangent
 * @param expression expression to interpolate
 * @param localValues destination local coefficients
 * @param edgeCount number of edges in the traversed entity
 * @param dofPerEdge number of ordered point DoFs per edge
 * @param tangent reusable tangent workspace
 * @param tangentProvider geometry policy called for each edge
 */
template<typename Expr, typename LocalVector, typename Tangent, typename TangentProvider>
void
interpolateHCurlEdgePointDofs( Expr const& expression,
                              LocalVector& localValues,
                              uint16_type edgeCount,
                              uint16_type dofPerEdge,
                              Tangent& tangent,
                              TangentProvider&& tangentProvider )
{
    using shape_type = typename std::remove_cvref_t<Expr>::shape;
    for ( uint16_type edge = 0; edge < edgeCount; ++edge )
    {
        tangentProvider( edge, tangent );
        for ( uint16_type ordinal = 0; ordinal < dofPerEdge; ++ordinal )
        {
            const uint16_type point = static_cast<uint16_type>( edge*dofPerEdge + ordinal );
            for ( uint16_type component = 0; component < shape_type::M; ++component )
                localValues( point ) += expression.evalq( component, 0, point ) * tangent( component );
        }
    }
}

} // namespace Feel::detail

#endif // FEELPP_FEELPOLY_HCURLINTERPOLATION_HPP
