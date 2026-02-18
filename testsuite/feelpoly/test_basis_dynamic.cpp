/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-01-06

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_basis_dynamic.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-06
   \brief Tests for dynamic evaluation in basis classes
 */

#define BOOST_TEST_MODULE test_basis_dynamic
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/dubiner.hpp>
#include <feel/feelpoly/legendre.hpp>

#include <cstddef>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( basis_dynamic_suite )

BOOST_AUTO_TEST_CASE( test_dubiner_dynamic_eval )
{
    using basis_t = Dubiner<2, 2, 10>;
    using points_t = typename basis_t::points_type;

    points_t pts( basis_t::nDim, 3 );
    pts( 0, 0 ) = 0.0;
    pts( 1, 0 ) = 0.0;
    pts( 0, 1 ) = 0.5;
    pts( 1, 1 ) = 0.0;
    pts( 0, 2 ) = 0.0;
    pts( 1, 2 ) = 0.5;

    auto full = basis_t::evaluate( pts );

    BOOST_CHECK_EQUAL( full.size1(), basis_t::convex_type::polyDims( basis_t::nOrder ) );
    BOOST_CHECK_EQUAL( full.size2(), pts.size2() );

    auto full_d = basis_t::derivate( pts );
    BOOST_CHECK_EQUAL( full_d.size(), basis_t::nDim );

    for ( uint16_type order : { uint16_type( 2 ), uint16_type( 5 ), uint16_type( 10 ) } )
    {
        auto reduced = basis_t::evaluate( pts, order );
        BOOST_CHECK_EQUAL( reduced.size1(), basis_t::convex_type::polyDims( order ) );
        BOOST_CHECK_EQUAL( reduced.size2(), pts.size2() );

        auto reduced_d = basis_t::derivate( pts, order );
        BOOST_CHECK_EQUAL( reduced_d.size(), basis_t::nDim );
        for ( std::size_t i = 0; i < reduced_d.size(); ++i )
        {
            BOOST_CHECK_EQUAL( reduced_d[i].size1(), basis_t::convex_type::polyDims( order ) );
            BOOST_CHECK_EQUAL( reduced_d[i].size2(), pts.size2() );
        }
    }
}

BOOST_AUTO_TEST_CASE( test_legendre_dynamic_eval )
{
    using basis_t = Legendre<1, 1, 10>;
    using points_t = typename basis_t::points_type;

    constexpr uint16_type n_dim = 1;
    constexpr uint16_type n_order = 10;
    points_t pts( n_dim, 4 );
    pts( 0, 0 ) = -0.5;
    pts( 0, 1 ) = 0.0;
    pts( 0, 2 ) = 0.5;
    pts( 0, 3 ) = 0.75;

    auto full = basis_t::evaluate( pts );

    BOOST_CHECK_EQUAL( full.size1(), basis_t::convex_type::polyDims( n_order ) );
    BOOST_CHECK_EQUAL( full.size2(), pts.size2() );

    auto full_d = basis_t::derivate( pts );
    BOOST_CHECK_EQUAL( full_d.size(), n_dim );

    for ( uint16_type order : { uint16_type( 2 ), uint16_type( 5 ), uint16_type( 10 ) } )
    {
        auto reduced = basis_t::evaluate( pts, order );
        BOOST_CHECK_EQUAL( reduced.size1(), basis_t::convex_type::polyDims( order ) );
        BOOST_CHECK_EQUAL( reduced.size2(), pts.size2() );

        auto reduced_d = basis_t::derivate( pts, order );
        BOOST_CHECK_EQUAL( reduced_d.size(), n_dim );
        for ( std::size_t i = 0; i < reduced_d.size(); ++i )
        {
            BOOST_CHECK_EQUAL( reduced_d[i].size1(), basis_t::convex_type::polyDims( order ) );
            BOOST_CHECK_EQUAL( reduced_d[i].size2(), pts.size2() );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
