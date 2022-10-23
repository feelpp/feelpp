/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Guillaume Dollé <gdolle@unistra.fr>
 Date: 26 Mar 2015

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
#define BOOST_TEST_MODULE test_worldcomm
#include <feel/feelcore/testsuite.hpp>

#include <vector>
#include <feel/feelcore/enumerate.hpp>

BOOST_AUTO_TEST_SUITE( enumeratesuite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    using namespace Feel;
    std::vector<int> v = {1,2,3,4,5};
    // enumerate the vector, index starts at 0
    for( auto const& [i, e] : enumerate(v) )
    {
        BOOST_CHECK_EQUAL( i, e-1 );
    }
}

BOOST_AUTO_TEST_CASE( test_start )
{
    using namespace Feel;
    std::vector<int> v = {1,2,3,4,5};
    // enumerate the vector, index starts at 1
    for( auto const& [i, e] : enumerate(v,1) )
    {
        BOOST_CHECK_EQUAL( i, e );
    }
}
BOOST_AUTO_TEST_SUITE_END()

