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

    rt0_type::points_type pts(2,5);
    pts( 0, 0 ) = -1./3.; pts( 1, 0 ) = -1./3.;
    pts( 0, 1 ) = -1.; pts( 1, 1 ) = -1.;
    pts( 0, 2 ) =  1.; pts( 1, 2 ) = -1.;
    pts( 0, 3 ) = -1.; pts( 1, 3 ) =  1.;
    pts( 0, 4 ) = -1.; pts( 1, 4 ) = 0.;
    //pts = rt0.referenceConvex().barycenterFaces();

    std::cout << "pts= " << pts << "\n";
    auto eval_at_pts = rt0.evaluate( pts );
    std::cout << "eval at pts= " << eval_at_pts << "\n";
}

BOOST_AUTO_TEST_CASE( rt1 )
{
    using namespace Life;
    typedef RaviartThomas<1>::apply<2>::type rt1_type;
    rt1_type rt1;

    rt1_type::points_type pts(2,1);
    pts( 0, 0 ) = -1./3.; pts( 1, 0 ) = -1./3.;

    //pts = rt1.referenceConvex().barycenterFaces();

    std::cout << "pts= " << pts << "\n";
    auto eval_at_pts = rt1.evaluate( pts );
    std::cout << "eval at pts= " << eval_at_pts << "\n";
    std::cout << "\n" ;
    std::cout << "##################################################################################################################################################"  << "\n" ;

    rt1_type::points_type Pts( 2 , 3 ) ;

    Pts( 0, 0 ) = -1 ; Pts( 1, 0 ) = -1 ;
    Pts( 0, 1 ) =  1.; Pts( 1, 1 ) = -1 ;
    Pts( 0, 2 ) = -1.; Pts( 1, 2 ) =  1 ;

    std::cout << "pts= " << Pts << "\n";
    auto eval_at_Pts = rt1.evaluate( Pts );
    std::cout << "eval at pts= " << eval_at_Pts << "\n";
    std::cout <<"\n" ;

   std::cout << "#####################################################################################################################################################"  << "\n" ;

    rt1_type::points_type opts( 2 , 5 ) ;

    opts( 0, 0 ) = -1.; opts( 1, 0 ) = -1. ;
    opts( 0, 1 ) =  1.; opts( 1, 1 ) = -1. ;
    opts( 0, 2 ) = -1.; opts( 1, 2 ) =  1. ;
    opts( 0, 3 ) = -1.; opts( 1, 3 ) =  0. ;
    opts( 0, 4 ) =  0.; opts( 1, 4 ) =  0. ;

   std::cout << "pts= " << opts << "\n";
    auto eval_at_opts = rt1.evaluate( opts );
    std::cout << "eval at pts= " << eval_at_opts << "\n";




}


BOOST_AUTO_TEST_SUITE_END()


