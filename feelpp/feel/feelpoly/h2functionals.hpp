/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/** @file h2functionals.hpp @brief Derivative functional blocks for Hermite/C1/H2 families. */

#ifndef FEELPP_FEELPOLY_H2FUNCTIONALS_HPP
#define FEELPP_FEELPOLY_H2FUNCTIONALS_HPP 1

#include <feel/feelpoly/pointfunctionals.hpp>

#include <iterator>
#include <vector>

namespace Feel::functional
{

/**
 * Build the value-and-gradient point jet used by Hermite-like and C1/H2
 * elements. Second/mixed and normal derivative blocks can be appended by the
 * family builder while sharing the same point-derivative primitives.
 *
 * @tparam Space scalar primal polynomial space
 * @param space primal space on which the jet functionals act
 * @param points points supporting value and first-derivative evaluations
 * @return ordered value block followed by one derivative block per direction
 */
template<typename Space>
[[nodiscard]] std::vector<Functional<Space>>
makeValueGradientPointFunctionals( Space const& space,
                                   typename Space::points_type const& points )
{
    std::vector<Functional<Space>> result;
    auto values = PointsEvaluation<Space>( space, points );
    result.insert( result.end(), values.begin(), values.end() );
    for ( uint16_type derivative = 0; derivative < Space::nDim; ++derivative )
    {
        auto derivatives = PointsDerivative<Space>( space, derivative, points );
        result.insert( result.end(), derivatives.begin(), derivatives.end() );
    }
    return result;
}

} // namespace Feel::functional

#endif // FEELPP_FEELPOLY_H2FUNCTIONALS_HPP
