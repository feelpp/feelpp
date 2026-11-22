/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-06-19

  Copyright (C) 2008-2012 University Joseph Fourier (Grenoble I)

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
   \file test_imsimplex.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-06-19
 */

#define BOOST_TEST_MODULE test_imsimplex
#include <boost/test/data/test_case.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/im.hpp>
#include <algorithm>
#include <cmath>

namespace bdata = boost::unit_test::data;
using namespace Feel;

static const std::vector<uint16_type> kOrders = {0, 1, 2, 3, 5, 10};

template<int Dim>
void check_imsimplex_order( uint16_type order )
{
  IMSimplex<Dim,double> im;
  im.create( order );

  BOOST_TEST_CONTEXT( "dim=" << Dim << " order=" << order )
  {
    BOOST_CHECK( im.nPoints() > 0 );
    BOOST_CHECK_EQUAL( im.weights().size(), im.nPoints() );

    if ( im.weights().size() > 0 )
      BOOST_CHECK_GT( im.weightsSum(), 0.0 );

    if constexpr ( Dim > 1 )
    {
      BOOST_CHECK( im.nFaces() > 0 );
      BOOST_CHECK( im.nPointsOnFace( 0 ) > 0 );
    }
  }
}

template<typename T>
void check_high_order_tetrahedron()
{
  for ( uint16_type order = 10; order <= 21; ++order )
  {
    IMGeneral<3,T,Simplex> im( order );
    BOOST_TEST_CONTEXT( "tetrahedron order=" << order << " scalar bytes=" << sizeof(T) )
    {
      T minimumWeight = im.weight(0);
      long double volume = 0;
      for ( int q = 0; q < im.nPoints(); ++q )
      {
        minimumWeight = std::min(minimumWeight,im.weight(q));
        volume += im.weight(q);
      }
      // Positive weights avoid cancellation in mass/energy assembly and norms.
      BOOST_REQUIRE_GT(minimumWeight,T(0));
      double const tolerance = sizeof(T) == sizeof(float) ? 5e-5 : 5e-12;
      BOOST_CHECK_SMALL(static_cast<double>(volume-4.L/3.L),tolerance);

      // Barycentric monomials span P_order on the reference tetrahedron.
      // Their exact integrals are 8*a!*b!*c!/(a+b+c+3)!.
      for ( int a = 0; a <= order; ++a )
        for ( int b = 0; a+b <= order; ++b )
          for ( int c = 0; a+b+c <= order; ++c )
          {
            long double value = 0;
            for ( int q = 0; q < im.nPoints(); ++q )
            {
              auto const x = (static_cast<long double>(im.points()(0,q))+1.L)/2.L;
              auto const y = (static_cast<long double>(im.points()(1,q))+1.L)/2.L;
              auto const z = (static_cast<long double>(im.points()(2,q))+1.L)/2.L;
              value += static_cast<long double>(im.weight(q))*
                       std::pow(x,a)*std::pow(y,b)*std::pow(z,c);
            }
            auto const exact = 8.L*std::tgamma(a+1.L)*std::tgamma(b+1.L)*
                               std::tgamma(c+1.L)/std::tgamma(a+b+c+4.L);
            BOOST_TEST_CONTEXT( "monomial=" << a << "," << b << "," << c )
              BOOST_CHECK_SMALL(static_cast<double>((value-exact)/exact),tolerance);
          }
    }
  }
}

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( imsimplex_suite )

BOOST_DATA_TEST_CASE( imsimplex_2d_orders, bdata::make( kOrders ), order )
{
  check_imsimplex_order<2>( order );
}

BOOST_DATA_TEST_CASE( imsimplex_3d_orders, bdata::make( kOrders ), order )
{
  check_imsimplex_order<3>( order );
}

BOOST_AUTO_TEST_CASE( high_order_tetrahedron_positive_exactness_double )
{
  check_high_order_tetrahedron<double>();
}

BOOST_AUTO_TEST_CASE( high_order_tetrahedron_positive_exactness_float )
{
  check_high_order_tetrahedron<float>();
}

BOOST_AUTO_TEST_SUITE_END()
