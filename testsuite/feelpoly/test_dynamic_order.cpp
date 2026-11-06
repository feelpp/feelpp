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
   \file test_dynamic_order.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-06
   \brief Tests for dynamic order handling in feelpoly
 */

#define BOOST_TEST_MODULE test_dynamic_order
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/dubiner.hpp>
#include <feel/feelpoly/polynomial.hpp>
#include <feel/feelpoly/polynomialset.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( dynamic_order_suite )

BOOST_AUTO_TEST_CASE( test_polynomialset_static_order )
{
    using poly_type = Dubiner<2, 2, 3>;
    PolynomialSet<poly_type, Scalar> ps;
    static_assert( PolynomialSet<poly_type, Scalar>::is_order_static, "Expected static order" );
    static_assert( PolynomialSet<poly_type, Scalar>::nOrder == 3, "Expected nOrder to match" );
    BOOST_CHECK_EQUAL( ps.degree(), 3 );
}

BOOST_AUTO_TEST_CASE( test_polynomialset_dynamic_order )
{
    using poly_type = Dubiner<2, 2, 3>;
    DynamicPolynomialSet<poly_type> ps_dyn( 10 );
    static_assert( DynamicPolynomialSet<poly_type>::is_order_dynamic, "Expected dynamic order" );
    BOOST_CHECK_EQUAL( ps_dyn.order(), 10 );
    BOOST_CHECK_EQUAL( ps_dyn.degree(), 10 );
}

BOOST_AUTO_TEST_CASE( test_polynomialset_order_conversion )
{
    using poly_type = Dubiner<2, 2, 3>;
    PolynomialSet<poly_type, Scalar> ps;
    auto dyn = ps.toDynamic();
    BOOST_CHECK_EQUAL( dyn.order(), ps.order() );

    auto as_static = dyn.toStatic<3>();
    BOOST_CHECK( static_cast<bool>( as_static ) );
    if ( as_static )
        BOOST_CHECK_EQUAL( as_static->degree(), 3 );

    auto wrong = dyn.toStatic<4>();
    BOOST_CHECK( !wrong );
}

BOOST_AUTO_TEST_CASE( test_polynomial_dynamic_order )
{
    using poly_type = Dubiner<2, 2, 3>;
    using poly_dyn_type = Polynomial<poly_type, Scalar, poly_type::basis_type::matrix_type, Dynamic>;
    poly_dyn_type p_dyn( 10 );
    BOOST_CHECK_EQUAL( p_dyn.order(), 10 );
    auto as_set = p_dyn.toSet( true );
    BOOST_CHECK_EQUAL( as_set.order(), 10 );
}

BOOST_AUTO_TEST_SUITE_END()
