/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/** @file tracefunctionals.hpp @brief Directional trace functional primitives. */

#ifndef FEELPP_FEELPOLY_TRACEFUNCTIONALS_HPP
#define FEELPP_FEELPOLY_TRACEFUNCTIONALS_HPP 1

#include <feel/feelpoly/pointfunctionals.hpp>

#include <algorithm>
#include <iterator>
#include <vector>

namespace Feel::functional
{

/**
 * @brief Append point evaluations contracted with a direction.
 * @tparam Space vector- or tensor-valued primal polynomial space
 * @param space primal space on which the functionals act
 * @param direction reference normal, tangent, or arbitrary direction
 * @param points functional support points
 * @param functionals ordered destination dual block
 */
template<typename Space>
void
appendDirectionalPointFunctionals(
    Space const& space,
    typename node<typename Space::value_type>::type const& direction,
    typename Space::points_type const& points,
    std::vector<Functional<Space>>& functionals )
{
    DirectionalComponentPointsEvaluation<Space> block( space, direction, points );
    std::copy( block.begin(), block.end(), std::back_inserter( functionals ) );
}

} // namespace Feel::functional

#endif // FEELPP_FEELPOLY_TRACEFUNCTIONALS_HPP
