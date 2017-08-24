//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 24 Aug 2017
//! @copyright 2017 Feel++ Consortium
//!

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE submesh testsuite
#include <testsuite.hpp>

#include <feel/feelcore/feel.hpp>

BOOST_AUTO_TEST_SUITE( cppenums )

enum class E
{
    TOTO=1,
    TUTU=2
};
BOOST_AUTO_TEST_CASE( test_to_underlying )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_to_underlying" );

    int toto = to_underlying( E::TOTO );
    BOOST_CHECK_EQUAL( toto, 1 );
    int tutu = to_underlying( E::TUTU );
    BOOST_CHECK_EQUAL( tutu, 2 );
    
    BOOST_TEST_MESSAGE( "test_to_underlying done" );
}

BOOST_AUTO_TEST_CASE( test_to_enum )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_to_enum" );

    E etoto = to_enum<E>( 1 );
    BOOST_CHECK_EQUAL( (int)etoto, (int)E::TOTO );
    E etutu = to_enum<E>( 2 );
    BOOST_CHECK_EQUAL( (int)etutu, (int)E::TUTU );
    
    BOOST_TEST_MESSAGE( "test_to_enum done" );
}

BOOST_AUTO_TEST_SUITE_END()


