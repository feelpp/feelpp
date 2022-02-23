/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2022-02-23

  Copyright (C) 2022 Université de Strasbourg

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
#define BOOST_TEST_MODULE hana testsuite
#include <fmt/core.h>
#include <fmt/compile.h>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>



FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( hanatests )

BOOST_AUTO_TEST_CASE( test_hana_structured_bindings )
{
    auto t = boost::hana::make_tuple( 1, 3.14, "Hello" );
    auto& [a, b, c] = t;
    BOOST_TEST_MESSAGE( "b=" << b );
    BOOST_CHECK_EQUAL( a, 1 );
    BOOST_CHECK_CLOSE( b, 3.14, 1e-14 );
    BOOST_CHECK_EQUAL( c, "Hello" );
}

BOOST_AUTO_TEST_CASE( test_hana_discr )
{
    using namespace Feel;
    auto dimt = hana::make_tuple( hana::int_c<1>, hana::int_c<2>, hana::int_c<3> );
    auto ordert = hana::unpack( hana::make_range( hana::int_c<0>, hana::int_c<FEELPP_INSTANTIATION_ORDER_MAX+1> ),hana::make_tuple);

    auto discretizationt = hana::transform( hana::cartesian_product( hana::make_tuple( dimt, ordert ) ), []( auto x )
                   { return hana::append( x, fmt::format( FMT_COMPILE( "P{}" ), std::decay_t<decltype( hana::at_c<1>( x ) )>::value ) ); } );

    for(int dimension=1; dimension <= 3; ++dimension)
        for ( int order = 0; order < FEELPP_INSTANTIATION_ORDER_MAX + 1; ++order)
        {
            std::string discretization=fmt::format("P{}",order);
            bool found = false;
            hana::for_each( discretizationt,
                   [&discretization, &dimension, &found]( auto const& d ) {
                                auto [_dim, _order,_discr] = d;
                                BOOST_TEST_MESSAGE( fmt::format( "checking dim:{} order:{}, space:{} // dimension: {} discretization: {}\n", _dim, _order, _discr,dimension,discretization ) );
                                if ( _dim == dimension && _discr == discretization )
                                    found = true;

                     } );
            BOOST_CHECK_MESSAGE( found, fmt::format("Discretization {} in dimension {} was not found!", discretization, dimension ) );
        }
}

BOOST_AUTO_TEST_SUITE_END()
