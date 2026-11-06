/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/** @file hcurlmomentpolicy.hpp @brief Shared H(curl) moment-order policies. */

#ifndef FEELPP_FEELPOLY_HCURLMOMENTPOLICY_HPP
#define FEELPP_FEELPOLY_HCURLMOMENTPOLICY_HPP 1

#include <feel/feelcore/feeltypes.hpp>

namespace Feel::detail
{

/**
 * @brief Quadrature order for triangular tangential polynomial moments.
 * @param polynomialDegree scalar triangular test-polynomial degree
 * @return quadrature order shared by dual construction and interpolation
 */
[[nodiscard]] constexpr uint16_type
hcurlTriangularFaceMomentQuadratureOrder( int polynomialDegree ) noexcept
{
    return static_cast<uint16_type>( 2*( polynomialDegree+3 ) );
}

} // namespace Feel::detail

#endif // FEELPP_FEELPOLY_HCURLMOMENTPOLICY_HPP
