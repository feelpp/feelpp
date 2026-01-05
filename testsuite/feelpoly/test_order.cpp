/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-01-04

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
   \file test_order.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-04
   \brief Tests for order system and DOF calculations in feelpoly
 */

#define BOOST_TEST_MODULE test_order
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/concepts.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( order_suite )

//
// Test order_t basic properties
//
BOOST_AUTO_TEST_CASE( test_order_t_static )
{
    // Static order values
    static_assert( order_t<0>::value == 0, "order_t<0>::value should be 0" );
    static_assert( order_t<1>::value == 1, "order_t<1>::value should be 1" );
    static_assert( order_t<5>::value == 5, "order_t<5>::value should be 5" );
    static_assert( order_t<10>::value == 10, "order_t<10>::value should be 10" );

    // Static/dynamic flags for static orders
    static_assert( order_t<0>::is_static, "order_t<0> should be static" );
    static_assert( order_t<3>::is_static, "order_t<3> should be static" );
    static_assert( !order_t<0>::is_dynamic, "order_t<0> should not be dynamic" );
    static_assert( !order_t<3>::is_dynamic, "order_t<3> should not be dynamic" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_order_t_dynamic )
{
    // Dynamic order
    static_assert( order_t<Dynamic>::value == Dynamic, "order_t<Dynamic>::value should be Dynamic" );
    static_assert( order_t<Dynamic>::value == -1, "order_t<Dynamic>::value should be -1" );

    // Static/dynamic flags for dynamic order
    static_assert( order_t<Dynamic>::is_dynamic, "order_t<Dynamic> should be dynamic" );
    static_assert( !order_t<Dynamic>::is_static, "order_t<Dynamic> should not be static" );

    BOOST_CHECK( true );
}

//
// Test binomial coefficient calculation
//
BOOST_AUTO_TEST_CASE( test_binomial )
{
    // Known binomial coefficients C(n, k) = n! / (k! * (n-k)!)
    // C(0, 0) = 1
    static_assert( detail::binomial( 0, 0 ) == 1, "C(0,0) = 1" );

    // C(n, 0) = 1 for all n
    static_assert( detail::binomial( 1, 0 ) == 1, "C(1,0) = 1" );
    static_assert( detail::binomial( 5, 0 ) == 1, "C(5,0) = 1" );
    static_assert( detail::binomial( 10, 0 ) == 1, "C(10,0) = 1" );

    // C(n, n) = 1 for all n
    static_assert( detail::binomial( 1, 1 ) == 1, "C(1,1) = 1" );
    static_assert( detail::binomial( 5, 5 ) == 1, "C(5,5) = 1" );
    static_assert( detail::binomial( 10, 10 ) == 1, "C(10,10) = 1" );

    // C(n, 1) = n
    static_assert( detail::binomial( 5, 1 ) == 5, "C(5,1) = 5" );
    static_assert( detail::binomial( 10, 1 ) == 10, "C(10,1) = 10" );

    // Pascal's triangle values
    static_assert( detail::binomial( 4, 2 ) == 6, "C(4,2) = 6" );
    static_assert( detail::binomial( 5, 2 ) == 10, "C(5,2) = 10" );
    static_assert( detail::binomial( 6, 3 ) == 20, "C(6,3) = 20" );
    static_assert( detail::binomial( 10, 5 ) == 252, "C(10,5) = 252" );

    // Invalid inputs
    static_assert( detail::binomial( 3, 5 ) == 0, "C(3,5) = 0 (k > n)" );
    static_assert( detail::binomial( 5, -1 ) == 0, "C(5,-1) = 0 (k < 0)" );

    BOOST_CHECK( true );
}

//
// Test pow_int calculation
//
BOOST_AUTO_TEST_CASE( test_pow_int )
{
    // base^0 = 1
    static_assert( detail::pow_int( 2, 0 ) == 1, "2^0 = 1" );
    static_assert( detail::pow_int( 5, 0 ) == 1, "5^0 = 1" );

    // base^1 = base
    static_assert( detail::pow_int( 2, 1 ) == 2, "2^1 = 2" );
    static_assert( detail::pow_int( 5, 1 ) == 5, "5^1 = 5" );

    // Powers of 2
    static_assert( detail::pow_int( 2, 2 ) == 4, "2^2 = 4" );
    static_assert( detail::pow_int( 2, 3 ) == 8, "2^3 = 8" );
    static_assert( detail::pow_int( 2, 4 ) == 16, "2^4 = 16" );
    static_assert( detail::pow_int( 2, 10 ) == 1024, "2^10 = 1024" );

    // Other bases
    static_assert( detail::pow_int( 3, 3 ) == 27, "3^3 = 27" );
    static_assert( detail::pow_int( 4, 3 ) == 64, "4^3 = 64" );
    static_assert( detail::pow_int( 5, 4 ) == 625, "5^4 = 625" );

    BOOST_CHECK( true );
}

//
// Test simplex DOF calculations
// DOF for simplex of dimension d with order p: C(d+p, d) = (d+p)! / (d! * p!)
//
BOOST_AUTO_TEST_CASE( test_simplex_dof_static )
{
    // 1D simplex (segment)
    // P0: 1 DOF, P1: 2 DOFs, P2: 3 DOFs, Pk: k+1 DOFs
    static_assert( detail::simplex_dof_static<1, 0>() == 1, "1D P0: 1 DOF" );
    static_assert( detail::simplex_dof_static<1, 1>() == 2, "1D P1: 2 DOFs" );
    static_assert( detail::simplex_dof_static<1, 2>() == 3, "1D P2: 3 DOFs" );
    static_assert( detail::simplex_dof_static<1, 3>() == 4, "1D P3: 4 DOFs" );
    static_assert( detail::simplex_dof_static<1, 5>() == 6, "1D P5: 6 DOFs" );

    // 2D simplex (triangle)
    // P0: 1, P1: 3, P2: 6, P3: 10, Pk: (k+1)(k+2)/2
    static_assert( detail::simplex_dof_static<2, 0>() == 1, "2D P0: 1 DOF" );
    static_assert( detail::simplex_dof_static<2, 1>() == 3, "2D P1: 3 DOFs" );
    static_assert( detail::simplex_dof_static<2, 2>() == 6, "2D P2: 6 DOFs" );
    static_assert( detail::simplex_dof_static<2, 3>() == 10, "2D P3: 10 DOFs" );
    static_assert( detail::simplex_dof_static<2, 4>() == 15, "2D P4: 15 DOFs" );
    static_assert( detail::simplex_dof_static<2, 5>() == 21, "2D P5: 21 DOFs" );

    // 3D simplex (tetrahedron)
    // P0: 1, P1: 4, P2: 10, P3: 20, Pk: (k+1)(k+2)(k+3)/6
    static_assert( detail::simplex_dof_static<3, 0>() == 1, "3D P0: 1 DOF" );
    static_assert( detail::simplex_dof_static<3, 1>() == 4, "3D P1: 4 DOFs" );
    static_assert( detail::simplex_dof_static<3, 2>() == 10, "3D P2: 10 DOFs" );
    static_assert( detail::simplex_dof_static<3, 3>() == 20, "3D P3: 20 DOFs" );
    static_assert( detail::simplex_dof_static<3, 4>() == 35, "3D P4: 35 DOFs" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_simplex_dof_dynamic )
{
    // Compare dynamic calculation with static values
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<1>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<1>( 1 ), 2 );
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<1>( 2 ), 3 );
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<1>( 5 ), 6 );

    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<2>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<2>( 1 ), 3 );
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<2>( 2 ), 6 );
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<2>( 3 ), 10 );

    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<3>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<3>( 1 ), 4 );
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<3>( 2 ), 10 );
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<3>( 3 ), 20 );

    // Invalid order
    BOOST_CHECK_EQUAL( detail::simplex_dof_dynamic<2>( -1 ), 0 );
}

//
// Test hypercube DOF calculations
// DOF for hypercube of dimension d with order p: (p+1)^d
//
BOOST_AUTO_TEST_CASE( test_hypercube_dof_static )
{
    // 1D hypercube (segment) - same as simplex
    static_assert( detail::hypercube_dof_static<1, 0>() == 1, "1D Q0: 1 DOF" );
    static_assert( detail::hypercube_dof_static<1, 1>() == 2, "1D Q1: 2 DOFs" );
    static_assert( detail::hypercube_dof_static<1, 2>() == 3, "1D Q2: 3 DOFs" );
    static_assert( detail::hypercube_dof_static<1, 3>() == 4, "1D Q3: 4 DOFs" );

    // 2D hypercube (quadrilateral)
    // Q0: 1, Q1: 4, Q2: 9, Q3: 16, Qk: (k+1)^2
    static_assert( detail::hypercube_dof_static<2, 0>() == 1, "2D Q0: 1 DOF" );
    static_assert( detail::hypercube_dof_static<2, 1>() == 4, "2D Q1: 4 DOFs" );
    static_assert( detail::hypercube_dof_static<2, 2>() == 9, "2D Q2: 9 DOFs" );
    static_assert( detail::hypercube_dof_static<2, 3>() == 16, "2D Q3: 16 DOFs" );
    static_assert( detail::hypercube_dof_static<2, 4>() == 25, "2D Q4: 25 DOFs" );

    // 3D hypercube (hexahedron)
    // Q0: 1, Q1: 8, Q2: 27, Q3: 64, Qk: (k+1)^3
    static_assert( detail::hypercube_dof_static<3, 0>() == 1, "3D Q0: 1 DOF" );
    static_assert( detail::hypercube_dof_static<3, 1>() == 8, "3D Q1: 8 DOFs" );
    static_assert( detail::hypercube_dof_static<3, 2>() == 27, "3D Q2: 27 DOFs" );
    static_assert( detail::hypercube_dof_static<3, 3>() == 64, "3D Q3: 64 DOFs" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_hypercube_dof_dynamic )
{
    // Compare dynamic calculation with static values
    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<1>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<1>( 1 ), 2 );
    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<1>( 2 ), 3 );

    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<2>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<2>( 1 ), 4 );
    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<2>( 2 ), 9 );
    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<2>( 3 ), 16 );

    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<3>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<3>( 1 ), 8 );
    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<3>( 2 ), 27 );
    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<3>( 3 ), 64 );

    // Invalid order
    BOOST_CHECK_EQUAL( detail::hypercube_dof_dynamic<2>( -1 ), 0 );
}

//
// Test DofCalculator with static order
//
BOOST_AUTO_TEST_CASE( test_dof_calculator_static )
{
    // 2D simplex, order 3
    using calc_2d_3 = DofCalculator<2, order_t<3>>;
    static_assert( calc_2d_3::simplex() == 10, "2D P3 simplex: 10 DOFs" );
    static_assert( calc_2d_3::hypercube() == 16, "2D Q3 hypercube: 16 DOFs" );

    // 3D simplex, order 2
    using calc_3d_2 = DofCalculator<3, order_t<2>>;
    static_assert( calc_3d_2::simplex() == 10, "3D P2 simplex: 10 DOFs" );
    static_assert( calc_3d_2::hypercube() == 27, "3D Q2 hypercube: 27 DOFs" );

    // 1D, order 5
    using calc_1d_5 = DofCalculator<1, order_t<5>>;
    static_assert( calc_1d_5::simplex() == 6, "1D P5: 6 DOFs" );
    static_assert( calc_1d_5::hypercube() == 6, "1D Q5: 6 DOFs" );

    BOOST_CHECK( true );
}

//
// Test DofCalculator with dynamic order
//
BOOST_AUTO_TEST_CASE( test_dof_calculator_dynamic )
{
    using calc_2d_dyn = DofCalculator<2, order_t<Dynamic>>;
    using calc_3d_dyn = DofCalculator<3, order_t<Dynamic>>;

    // 2D tests
    BOOST_CHECK_EQUAL( calc_2d_dyn::simplex( 0 ), 1 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::simplex( 1 ), 3 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::simplex( 2 ), 6 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::simplex( 3 ), 10 );

    BOOST_CHECK_EQUAL( calc_2d_dyn::hypercube( 0 ), 1 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::hypercube( 1 ), 4 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::hypercube( 2 ), 9 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::hypercube( 3 ), 16 );

    // 3D tests
    BOOST_CHECK_EQUAL( calc_3d_dyn::simplex( 0 ), 1 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::simplex( 1 ), 4 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::simplex( 2 ), 10 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::simplex( 3 ), 20 );

    BOOST_CHECK_EQUAL( calc_3d_dyn::hypercube( 0 ), 1 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::hypercube( 1 ), 8 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::hypercube( 2 ), 27 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::hypercube( 3 ), 64 );
}

//
// Test order_t satisfies concepts
//
BOOST_AUTO_TEST_CASE( test_order_concepts_satisfaction )
{
    // Static orders should satisfy StaticOrder concept
    static_assert( StaticOrder<order_t<0>>, "order_t<0> should satisfy StaticOrder" );
    static_assert( StaticOrder<order_t<1>>, "order_t<1> should satisfy StaticOrder" );
    static_assert( StaticOrder<order_t<5>>, "order_t<5> should satisfy StaticOrder" );

    // Dynamic order should satisfy DynamicOrder concept
    static_assert( DynamicOrder<order_t<Dynamic>>, "order_t<Dynamic> should satisfy DynamicOrder" );

    // Cross-checks
    static_assert( !DynamicOrder<order_t<0>>, "order_t<0> should NOT satisfy DynamicOrder" );
    static_assert( !StaticOrder<order_t<Dynamic>>, "order_t<Dynamic> should NOT satisfy StaticOrder" );

    // All should satisfy HasOrder
    static_assert( HasOrder<order_t<0>>, "order_t<0> should satisfy HasOrder" );
    static_assert( HasOrder<order_t<5>>, "order_t<5> should satisfy HasOrder" );
    static_assert( HasOrder<order_t<Dynamic>>, "order_t<Dynamic> should satisfy HasOrder" );

    BOOST_CHECK( true );
}

//
// Test comparison of simplex vs hypercube DOFs
//
BOOST_AUTO_TEST_CASE( test_simplex_vs_hypercube_dofs )
{
    // For same dimension and order, hypercube always has >= DOFs than simplex
    // (except in 1D where they are equal)

    // 1D: equal
    static_assert( detail::simplex_dof_static<1, 3>() == detail::hypercube_dof_static<1, 3>(),
                   "1D: simplex and hypercube have same DOFs" );

    // 2D: hypercube has more
    static_assert( detail::hypercube_dof_static<2, 2>() > detail::simplex_dof_static<2, 2>(),
                   "2D: hypercube has more DOFs (9 > 6)" );
    static_assert( detail::hypercube_dof_static<2, 3>() > detail::simplex_dof_static<2, 3>(),
                   "2D: hypercube has more DOFs (16 > 10)" );

    // 3D: hypercube has more
    static_assert( detail::hypercube_dof_static<3, 2>() > detail::simplex_dof_static<3, 2>(),
                   "3D: hypercube has more DOFs (27 > 10)" );
    static_assert( detail::hypercube_dof_static<3, 3>() > detail::simplex_dof_static<3, 3>(),
                   "3D: hypercube has more DOFs (64 > 20)" );

    BOOST_CHECK( true );
}

//
// Test high-order DOF calculations
//
BOOST_AUTO_TEST_CASE( test_high_order_dofs )
{
    // Test that high-order calculations don't overflow for reasonable orders
    // 2D simplex P10: C(12, 2) = 66
    static_assert( detail::simplex_dof_static<2, 10>() == 66, "2D P10: 66 DOFs" );

    // 3D simplex P5: C(8, 3) = 56
    static_assert( detail::simplex_dof_static<3, 5>() == 56, "3D P5: 56 DOFs" );

    // 2D hypercube Q10: 11^2 = 121
    static_assert( detail::hypercube_dof_static<2, 10>() == 121, "2D Q10: 121 DOFs" );

    // 3D hypercube Q5: 6^3 = 216
    static_assert( detail::hypercube_dof_static<3, 5>() == 216, "3D Q5: 216 DOFs" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
