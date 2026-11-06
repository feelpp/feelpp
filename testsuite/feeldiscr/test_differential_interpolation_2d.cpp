/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#define BOOST_TEST_MODULE differential interpolation 2d
#include "differential_interpolation_checks.hpp"
using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_CASE( differential_matrices_2d )
{
    differential_interpolation_checks::checkDimension<2>();
}
