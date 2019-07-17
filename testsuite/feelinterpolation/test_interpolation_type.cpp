/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 18 août 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#define BOOST_TEST_MODULE interpolation_type testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>

BOOST_AUTO_TEST_SUITE( interpolation_type_suite )

BOOST_AUTO_TEST_CASE( test_interpolation_type_1 )
{
    BOOST_MESSAGE( "================================================================================" );
    BOOST_MESSAGE( "== test_interpolation_type" );
    using namespace Feel;
    auto id_c = makeInterpolation( conforming_t() );
    BOOST_CHECK_EQUAL( id_c.isConforming(), true );
    auto grad_c = makeGradientInterpolation( conforming_t() );
    BOOST_CHECK_EQUAL( grad_c.isConforming(), true );

    auto id_nc = makeInterpolation( nonconforming_t() );
    BOOST_CHECK_EQUAL( id_nc.isConforming(), false );
    auto grad_nc = makeGradientInterpolation( nonconforming_t() );
    BOOST_CHECK_EQUAL( grad_nc.isConforming(), false );

    BOOST_CHECK_EQUAL( id_c.searchWithCommunication(), true );
    BOOST_CHECK_EQUAL( id_c.componentsAreSamePoint(), true );
    BOOST_CHECK_EQUAL( id_c.onlyLocalizeOnBoundary(), false );
    BOOST_CHECK_EQUAL( id_c.nbNearNeighborInKdTree(), 15 );

    BOOST_CHECK_EQUAL( id_c.interpolationOperand(), interpolation_operand_type::ID );
    BOOST_CHECK_EQUAL( grad_c.interpolationOperand(), interpolation_operand_type::GRADIENT );
    BOOST_CHECK_EQUAL( id_nc.interpolationOperand(), interpolation_operand_type::ID );
    BOOST_CHECK_EQUAL( grad_nc.interpolationOperand(), interpolation_operand_type::GRADIENT );

}


BOOST_AUTO_TEST_CASE( test_interpolation_type_2 )
{
    BOOST_MESSAGE( "================================================================================" );
    BOOST_MESSAGE( "== test_interpolation_type 2" );
    using namespace Feel;
    auto id_c = makeInterpolation( conforming_t(), false, false, true, 5 );
    BOOST_CHECK_EQUAL( id_c.isConforming(), true );
    auto grad_c = makeGradientInterpolation( conforming_t(), false, false, true, 5 );
    BOOST_CHECK_EQUAL( grad_c.isConforming(), true );

    auto id_nc = makeInterpolation( nonconforming_t(), false, false, true, 5 );
    BOOST_CHECK_EQUAL( id_nc.isConforming(), false );
    auto grad_nc = makeGradientInterpolation( nonconforming_t(), false, false, true, 5 );
    BOOST_CHECK_EQUAL( grad_nc.isConforming(), false );

    BOOST_CHECK_EQUAL( id_c.searchWithCommunication(), false);
    BOOST_CHECK_EQUAL( id_c.componentsAreSamePoint(), false );
    BOOST_CHECK_EQUAL( id_c.onlyLocalizeOnBoundary(), true );
    BOOST_CHECK_EQUAL( id_c.nbNearNeighborInKdTree(), 5 );

    BOOST_CHECK_EQUAL( id_c.interpolationOperand(), interpolation_operand_type::ID );
    BOOST_CHECK_EQUAL( grad_c.interpolationOperand(), interpolation_operand_type::GRADIENT );
    BOOST_CHECK_EQUAL( id_nc.interpolationOperand(), interpolation_operand_type::ID );
    BOOST_CHECK_EQUAL( grad_nc.interpolationOperand(), interpolation_operand_type::GRADIENT );
}

BOOST_AUTO_TEST_SUITE_END()
