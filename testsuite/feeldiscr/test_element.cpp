/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-03-03

  Copyright (C) 2013 Feel++ Consortium

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
/**
   \file test_element.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-03-03
 */

#define USE_BOOST_TEST 1
// Boost.Test

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE element testsuite
#include <testsuite/testsuite.hpp>
#include <feel/feel.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( element )

BOOST_AUTO_TEST_CASE( test_element_1 )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_element_1" );

    auto mesh= unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = Xh->element();
    u.setOnes();
    auto a = form2( _test=Xh, _trial=Xh );
    a = integrate( elements(mesh), idt(u)*id(u));
    a.close();
    auto v = backend()->newVector( Xh );
    auto w = backend()->newVector( Xh );
    v->setOnes();
    a.matrixPtr()->multVector( u, *w );
    BOOST_CHECK_CLOSE( v->dot( w ), 1, 1.5e-13 );
    a.matrixPtr()->multVector( *v, *w );
    BOOST_CHECK_CLOSE( v->dot( w ), 1, 1.5e-13 );

    BOOST_TEST_MESSAGE( "test_element_1 done" );
}

BOOST_AUTO_TEST_SUITE_END()

