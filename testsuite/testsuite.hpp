/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-07-16

  Copyright (C) 2006 EPFL

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
   \file testsuite.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-07-16
 */
#if !defined (__TESTSUITE_HPP__)
#define __TESTSUITE_HPP__ 1

#if defined( USE_BOOST_TEST )
// Boost.Test
//#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/results_reporter.hpp>
#include <boost/test/output_test_stream.hpp>
#include <boost/test/unit_test_log.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/framework.hpp>
#include <boost/test/floating_point_comparison.hpp>
namespace Life
{
namespace ut = boost::unit_test;
static ut::test_suite* testsuite_master = BOOST_TEST_SUITE( "Life Master Test Suite" );
static ut::test_suite* testsuite_core = BOOST_TEST_SUITE( "Life Core Test Suite" );
static ut::test_suite* testsuite_poly = BOOST_TEST_SUITE( "Life Polynomials Test Suite" );
static ut::test_suite* testsuite_mesh = BOOST_TEST_SUITE( "Life Mesh Test Suite" );
static ut::test_suite* testsuite_fem = BOOST_TEST_SUITE( "Life FEM Test Suite" );

}
#endif

#endif //__TESTSUITE_HPP__
