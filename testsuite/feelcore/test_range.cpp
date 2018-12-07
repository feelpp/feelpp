/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 17 août 2015
 
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

#define BOOST_TEST_MODULE range testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/range.hpp>


BOOST_AUTO_TEST_SUITE( range_testsuite )

BOOST_AUTO_TEST_CASE( test_range1 )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_range1" );
    
    std::vector<int> v( 10, 0 );
    std::iota( v.begin(), v.end(), 0 );
    for( auto i : range( 10 ) )
        BOOST_CHECK_EQUAL( i, v[i] );
    
    BOOST_TEST_MESSAGE( "test_range1 done" );
}

BOOST_AUTO_TEST_CASE( test_range2 )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_range2" );
    
    std::vector<int> v( 10, 0 );
    std::iota( v.begin(), v.end(), -2 );
    for( auto i : range( -2, 10-2 ) )
        BOOST_CHECK_EQUAL( i, v[i+2] );
    
    BOOST_TEST_MESSAGE( "test_range2 done" );
}

BOOST_AUTO_TEST_CASE( test_range3 )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_range3" );
    
    std::vector<int> v( 10/2, 0 );
    int n = -2;
    std::generate( v.begin(), v.end(), [&n]() { n+=2; return n-2; } );
    auto j = 0;
    for( auto i : range( -2, 10-2, 2 ) )
        BOOST_CHECK_EQUAL( i, v[j++] );
    
    BOOST_TEST_MESSAGE( "test_range3 done" );
}

BOOST_AUTO_TEST_SUITE_END()

