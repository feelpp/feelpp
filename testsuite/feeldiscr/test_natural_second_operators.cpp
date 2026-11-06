/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#define BOOST_TEST_MODULE natural second operators
#include "natural_second_operator_checks.hpp"
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/pchv.hpp>
using namespace Feel;
using namespace Feel::natural_second_checks;
FEELPP_ENVIRONMENT_NO_OPTIONS

template <int Dim>
void checkDimension()
{
    for ( uint16_type geometryOrder : { 1, 2, 3 } )
    {
        auto mesh = box<Dim>( geometryOrder );
        checkPolynomial( Pchv<2>( mesh ) );
        checkPolynomial( Pchv<Dynamic>( mesh, RuntimeOrder{ 2 } ) );
    }
    for ( uint16_type geometryOrder : { 2, 3 } )
    {
        auto mesh = curved<Dim>( geometryOrder );
        for ( uint16_type order : { 0, 1, 3 } )
            checkCurved( Pchv<Dynamic>( mesh, RuntimeOrder{ order } ) );
    }
    checkPolynomial( Pchv<2>( staticBox<Dim>() ) );
}
BOOST_AUTO_TEST_CASE( nodal_2d )
{
    checkDimension<2>();
}
BOOST_AUTO_TEST_CASE( nodal_3d )
{
    checkDimension<3>();
}
