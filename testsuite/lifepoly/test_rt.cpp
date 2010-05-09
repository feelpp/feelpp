/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-05-08

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file test_rt.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-05-08
 */
#define BOOST_TEST_MODULE Raviar-Thomas polynomials test
// Boost.Test
#define USE_TEST 1
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <life/lifepoly/raviartthomas.hpp>

BOOST_AUTO_TEST_SUITE( rt_testsuite )

BOOST_AUTO_TEST_CASE( rt0 )
{
    using namespace Life;
    typedef RaviartThomas<0>::apply<2>::type rt0_type;

    rt0_type rt0;

    rt0_type::points_type pts(2,3);
    pts = rt0.referenceConvex().vertices();

    std::cout << "pts= " << pts << "\n";
    auto eval_at_pts = rt0.evaluate( pts );
    std::cout << "eval at pts= " << eval_at_pts << "\n";

}
BOOST_AUTO_TEST_SUITE_END()


