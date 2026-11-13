/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#define BOOST_TEST_MODULE Piola second operators 2d
#include "piola_second_operator_checks.hpp"
#include <feel/feelcore/testsuite.hpp>
using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_CASE( piola_2d )
{
    natural_second_checks::checkPiolaDimension<2>();
}
