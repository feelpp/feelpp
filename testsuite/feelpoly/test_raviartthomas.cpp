/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-05-08
 */
#define BOOST_TEST_MODULE Raviar-Thomas polynomials test
// Boost.Test
#define USE_TEST 1
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <feel/feelpoly/raviartthomas.hpp>
#include <polynomial_set.hpp>

BOOST_AUTO_TEST_SUITE( rt_testsuite )

BOOST_AUTO_TEST_CASE( ravth_2d_0 )
{
    using namespace Feel;
    typedef RaviartThomas<0>::apply<2>::type rt0_type;

    rt0_type rt0;

    typedef  boost::numeric::ublas::matrix<double> matrix_type ;
    typedef  boost::numeric::ublas::vector<double> Points_type ;
    typedef  boost::numeric::ublas::vector<Points_type> vectors_type ;

    matrix_type m( 6 ,3 )   ;

    m( 0 ,0 )=m( 0 ,1 )=m( 1 ,0 )=m( 1 ,2 )=m( 2 ,1 )=m( 3 ,0 )=m( 3 ,2 )=m( 4 ,0 )=m( 4 ,1 )=m( 5 ,2 ) = 1./4.;
    m( 0 ,2 )=m( 1 ,1 )=m( 2 ,2 )=m( 3 ,1 )=m( 4 ,2 )=m( 5 ,1 ) = 0. ;
    m( 2 ,0 )=m( 5 ,0 )= -1./4. ;

    Points_type p1( 2 );
    p1( 0 )=-1./3.;
    p1( 1 ) =-1./3.;
    Points_type p2( 2 );
    p2( 0 )=-1.;
    p2( 1 ) =-1.;
    Points_type p3( 2 );
    p3( 0 )= 1.;
    p3( 1 ) =-1.;
    Points_type p4( 2 );
    p4( 0 )=-1.;
    p4( 1 ) = 1.;
    Points_type p5( 2 );
    p5( 0 )=-1.;
    p5( 1 ) = 0.;
    Points_type p6( 2 );
    p6( 0 )=-1./3.;
    p6( 1 ) =-1./3.;
    Points_type p7( 2 );
    p7( 0 )= 0.;
    p7( 1 ) = 0.;
    Points_type p8( 2 );
    p8( 0 )= 0.;
    p8( 1 ) =-1.;
    Points_type p9( 2 );
    p9( 0 )=-1./2.;
    p9( 1 ) = 1./2.;
    Points_type p10( 2 );
    p10( 0 )= 1./2.;
    p10( 1 ) =-1./2.;
    vectors_type P( 10 ) ;

    P( 0 )=p1;
    P( 1 )=p2;
    P( 2 )=p3;
    P( 3 )=p4;
    P( 4 )=p5;
    P( 5 )=p6;
    P( 6 )=p7;
    P( 7 )=p8;
    P( 8 )=p9;
    P( 9 )=p10 ;

    Polynomialset<10,6,3,10> Poly( m ) ;

    BOOST_TEST_MESSAGE(  "**************************************************Analytic method*********************************************" << "\n" );
    BOOST_TEST_MESSAGE( "coefficients= " << m << "\n" );
    BOOST_TEST_MESSAGE( "Points= " << P << "\n" );
    BOOST_TEST_MESSAGE( "evaluate_at_Points= " <<Poly.evaluates_Points ( P ) << "\n" );


    BOOST_TEST_MESSAGE(  "*****************************************Feel method***********************************************************" << "\n" );

    rt0_type::points_type pts( 2,10 );
    pts( 0 ,0 ) = -1./3. ;
    pts( 1 ,0 ) = -1./3. ;
    pts( 0 ,1 ) = -1. ;
    pts( 1 ,1 ) = -1. ;
    pts( 0 ,2 ) =  1. ;
    pts( 1 ,2 ) = -1. ;
    pts( 0 ,3 ) = -1. ;
    pts( 1 ,3 ) =  1. ;
    pts( 0 ,4 ) = -1. ;
    pts( 1 ,4 ) =  0. ;
    pts( 0 ,5 )=-1./3.;
    pts( 1 ,5 ) =-1./3.;
    pts( 0 ,6 )= 0.;
    pts( 1 ,6 ) = 0.;
    pts( 0 ,7 )= 0.;
    pts( 1 ,7 ) =-1.;
    pts( 0 ,8 )=-1./2.;
    pts( 1 ,8 ) = 1./2.;
    pts( 0 ,9 )= 1./2.;
    pts( 1 ,9 ) =-1./2.;

    BOOST_TEST_MESSAGE(  "pts= " << pts << "\n" );
    auto eval_at_pts = rt0.evaluate( pts );
    BOOST_TEST_MESSAGE(  "eval_at_pts= "<< eval_at_pts << "\n" );
    matrix_type diff( Poly.evaluates_Points( P )-eval_at_pts );
    double error = ublas::norm_frobenius( diff );
    BOOST_CHECK_SMALL( error, 1e-10 );


}



BOOST_AUTO_TEST_CASE( ravth_3d_0 )
{
    using namespace Feel;
    typedef RaviartThomas<0>::apply<3>::type rt0_type;
    typedef GeoND<3,Simplex<3, 1, 3> >::point_type point_type;
    typedef  boost::numeric::ublas::matrix<double> matrix_type ;
    typedef  boost::numeric::ublas::vector<double> Points_type ;
    typedef  boost::numeric::ublas::vector<Points_type> vectors_type ;

    GeoND<3,Simplex<3, 1, 3> > tetra;
    point_type V1;
    V1( 0 )=-1;
    V1( 1 )=-1;
    V1( 2 )=-1;
    point_type V2;
    V2( 0 )=1;
    V2( 1 )=-1;
    V2( 2 )=-1;
    point_type V3;
    V3( 0 )=-1;
    V3( 1 )=1;
    V3( 2 )=-1;
    point_type V4;
    V4( 0 )=-1;
    V4( 1 )=-1;
    V4( 2 )=1;
    tetra.setPoint( 0, V1 );
    tetra.setPoint( 1, V2 );
    tetra.setPoint( 2, V3 );
    tetra.setPoint( 3, V4 );
    tetra.update();

    BOOST_TEST_MESSAGE(  "[tetra] barycenter = " << tetra.barycenter() << "\n" );
    BOOST_TEST_MESSAGE(  "[tetra] measure = " << tetra.measure()<< "\n" );
    BOOST_TEST_MESSAGE(  tetra.faceMeasure( 0 ) << "\n" );
    BOOST_TEST_MESSAGE(  tetra.faceMeasure( 1 ) << "\n" );
    BOOST_TEST_MESSAGE(  tetra.faceMeasure( 2 ) << "\n" );
    BOOST_TEST_MESSAGE(  tetra.faceMeasure( 3 ) << "\n" );

    rt0_type rt0;

    matrix_type m( 12 ,4 )   ;

    m( 0 ,0 )=m( 0 ,1 )=m( 1 ,0 )=m( 1 ,2 )=m( 2 ,0 )=m( 2 ,3 )=m( 3 ,1 )=m( 4 ,0 )
                                            =m( 4 ,2 )=m( 5 ,0 )=m( 5 ,3 )=m( 6,0 ) = m( 6,1 ) =m( 7 ,2 ) =m( 8 , 0 )
                                                    =m( 8, 3 )=m( 9,0 ) = m( 9 ,1 ) =m( 10,0 ) =m( 10 ,2 )
                                                            =m( 11,3 ) =  1./( 3*tetra.measure() ) ;
    m( 0 ,2 )=  m( 0 ,3 )= m( 1 ,1 )= m( 1 ,3 )= m( 2 ,1 )=m( 2 ,2 )
                                      =m( 3 ,2 )= m( 3 ,3 )= m( 4 ,1 )=  m( 4 ,3 )=m( 5 ,1 )= m( 5 ,2 )
                                              = m( 6 ,2 )= m( 6 ,3 ) = m( 7 ,1 )= m( 7 ,3 )= m( 8 ,1 )= m( 8 ,2 )=
                                                      m( 9 ,2 )= m( 9 ,3 )= m( 10 ,1 )= m( 10 ,3 )= m( 11 ,1 )= m( 11 ,2 )= 0. ;
    m( 3 ,0 )=m( 7 ,0 )=m( 11 , 0 ) = -1./( 3*tetra.measure() ) ;

    Points_type p1( 3 );
    p1( 0 )=-1.;
    p1( 1 ) =-1.;
    p1( 2 ) = -1. ;
    Points_type p2( 3 );
    p2( 0 )= 1.;
    p2( 1 ) =-1.;
    p2( 2 ) = -1. ;
    Points_type p3( 3 );
    p3( 0 )= -1.;
    p3( 1 ) =1.;
    p3( 2 ) = -1. ;
    Points_type p4( 3 );
    p4( 0 )=-1.;
    p4( 1 ) =-1.;
    p4( 2 ) =  1. ;
    Points_type p5( 3 );
    p5( 0 )=-1.;
    p5( 1 ) = 0.;
    p5( 2 ) = -1. ;
    Points_type p6( 3 );
    p6( 0 )=-1.;
    p6( 1 ) =-1.;
    p6( 2 ) = 0.;
    Points_type p7( 3 );
    p7( 0 )= -1./3.;
    p7( 1 ) =-1./3.;
    p7( 2 ) =-1./3.;
    Points_type p8( 3 );
    p8( 0 )= -1.;
    p8( 1 ) =-1./3.;
    p8( 2 ) =-1./3.;
    Points_type p9( 3 );
    p9( 0 )=-1./3.;
    p9( 1 ) = -1.;
    p9( 2 ) = -1./3.;
    Points_type p10( 3 );
    p10( 0 )= -1./3.;
    p10( 1 ) =-1./3.;
    p10( 2 ) =-1.;

    vectors_type P( 10 ) ;

    P( 0 )=p1;
    P( 1 )=p2;
    P( 2 )=p3;
    P( 3 )=p4;
    P( 4 )=p5;
    P( 5 )=p6;
    P( 6 )=p7;
    P( 7 )=p8;
    P( 8 )=p9;
    P( 9 )=p10 ;
    Polynomialset<20,12,4,10> Poly( m ) ;

    BOOST_TEST_MESSAGE(  "**************************************************Analytic method*********************************************" << "\n" );
    BOOST_TEST_MESSAGE( "coefficients= " << m << "\n" );
    BOOST_TEST_MESSAGE( "\n" );
    BOOST_TEST_MESSAGE( "Points= " << P << "\n" );
    BOOST_TEST_MESSAGE( "\n" );
    BOOST_TEST_MESSAGE( "evaluate_at_Points= " <<Poly.evaluates_Points ( P ) << "\n" );
    BOOST_TEST_MESSAGE( "\n" );

    BOOST_TEST_MESSAGE(  "*****************************************Feel method***********************************************************" << "\n" );

    rt0_type::points_type pts( 3,10 );

    pts( 0 ,0 ) = -1. ;
    pts( 1 ,0 ) = -1. ;
    pts( 2 ,0 ) = -1. ;
    pts( 0 ,1 ) =  1. ;
    pts( 1 ,1 ) = -1. ;
    pts( 2 ,1 ) = -1. ;
    pts( 0 ,2 ) = -1. ;
    pts( 1 ,2 ) =  1. ;
    pts( 2 ,2 ) =-1. ;
    pts( 0 ,3 ) = -1. ;
    pts( 1 ,3 ) = -1. ;
    pts( 2 ,3 ) = 1. ;
    pts( 0 ,4 ) = -1.;
    pts( 1 ,4 ) = 0.  ;
    pts( 2 ,4 ) = -1. ;
    pts( 0 ,5 ) = -1.;
    pts( 1 ,5 ) = -1.;
    pts( 2 ,5 ) =  0. ;
    pts( 0 ,6 )=-1./3.;
    pts( 1 ,6 ) = -1./3.;
    pts( 2 ,6 ) = -1./3. ;
    pts( 0 ,7 )= -1.;
    pts( 1 ,7 ) = -1./3.;
    pts( 2 ,7 ) = -1./3. ;
    pts( 0 ,8 )= -1./3.;
    pts( 1 ,8 ) = -1.;
    pts( 2 ,8 ) = -1./3. ;
    pts( 0 ,9 )= -1./3.;
    pts( 1 ,9 ) =-1./3.;
    pts( 2 ,9 ) = -1. ;

    BOOST_TEST_MESSAGE(  "pts= " << pts << "\n" );
    BOOST_TEST_MESSAGE( "\n" );
    auto eval_at_pts = rt0.evaluate( pts );
    BOOST_TEST_MESSAGE(  "eval_at_pts= "<< eval_at_pts << "\n" );
    double error = ublas::norm_frobenius( Poly.evaluates_Points( P )-eval_at_pts );
    BOOST_CHECK_SMALL ( error, 1e-14 )  ;


}
BOOST_AUTO_TEST_CASE( ravth_2d_1 )
{
    Feel::Assert::setLog( "test_ravth.log" );
    using namespace Feel;
    BOOST_TEST_MESSAGE( "start" );
    typedef RaviartThomas<1>::apply<2>::type rt0_type;
    typedef  boost::numeric::ublas::matrix<double> matrix_type ;
    typedef  boost::numeric::ublas::vector<double> Points_type ;
    typedef  boost::numeric::ublas::vector<Points_type> vectors_type ;

    rt0_type rt0;
    BOOST_TEST_MESSAGE( "RT_1 instantiated" );

    matrix_type m( 16 ,6 )   ;

    /* for(uint k = 0 ; k < 12 ;++ k )
       {
       m(k , 3) = m(k , 4) = m(k ,5) =0 ;
       }
    */
    m( 0 ,0 )=m( 0 ,1 )=m( 1 ,0 )=m( 1 ,2 )=m( 2 ,1 )=m( 3 ,0 )=m( 3 ,2 )
                                            =m( 4 ,0 )=m( 4 ,1 )=m( 5 ,2 )  = 1./4. ;
    m( 2 ,0 )=m( 5 ,0 )= -1./4. ;
    m( 6 ,0 )=m( 6 ,1 )= m( 13 ,0 )=m( 14 ,0 ) = -3./4. ;
    m( 7 ,0 )=m( 7 ,2 )= 3./4. ;

    m( 8 ,0 )=m( 8 ,1 )= m( 10 ,0 )=m( 10 ,1 )= m( 12 ,0 )=m( 15 ,0 ) = -6./4. ;

    m( 9 ,0 )=m( 9 ,2 )= m( 11 ,0 )=m( 11 ,2 )= 6./4. ;

    m( 11 ,1 )= 12./4. ;

    m( 12 ,1 )=m( 15 ,2 )= -8./4. ;
    m( 12 ,2 )=m( 12 ,3 )= m( 15 ,3 )= -4./4. ;
    m( 12 ,4 )=m( 13 ,5 )= m( 14 ,4 )=m( 15 ,5 )= -2./4. ;
    m( 13 ,2 ) = -5./4. ;

    for ( uint k = 0 ; k < 12 ; ++ k )
    {
        m( k , 3 ) = m( k , 4 ) = m( k ,5 ) =0 ;
    }

    BOOST_TEST_MESSAGE( "coeff matrix created " );
    // BOOST_TEST_MESSAGE( " coefficients= " << m << "\n" );

    /*
      m(0 ,0)=m(0 ,1)=m(1 ,0)=m(1 ,2)=m(2 ,1)=m(3 ,0)=m(3 ,2)=m(4 ,0)
      =m(4 ,1)=m(5 ,2) = 1./4.;
      m(0 ,2)=m(1 ,1)=m(2 ,2)=m(3 ,1)=m(4 ,2)=m(5 ,1) = 0. ;
      m(2 ,0)=m(5 ,0)= -1./4. ;
    */

    Points_type p1( 2 );
    p1( 0 )=-1./3.;
    p1( 1 ) =-1./3.;
    Points_type p2( 2 );
    p2( 0 )=-1.;
    p2( 1 ) =-1.;
    Points_type p3( 2 );
    p3( 0 )= 1.;
    p3( 1 ) =-1.;
    Points_type p4( 2 );
    p4( 0 )=-1.;
    p4( 1 ) = 1.;
    Points_type p5( 2 );
    p5( 0 )=-1.;
    p5( 1 ) = 0.;
    Points_type p6( 2 );
    p6( 0 )=-1./3.;
    p6( 1 ) =-1./3.;
    Points_type p7( 2 );
    p7( 0 )= 0.;
    p7( 1 ) = 0.;
    Points_type p8( 2 );
    p8( 0 )= 0.;
    p8( 1 ) =-1.;
    Points_type p9( 2 );
    p9( 0 )=-1./2.;
    p9( 1 ) = 1./2.;
    Points_type p10( 2 );
    p10( 0 )= 1./2.;
    p10( 1 ) =-1./2.;
    vectors_type P( 10 ) ;

    P( 0 )=p1;
    P( 1 )=p2;
    P( 2 )=p3;
    P( 3 )=p4;
    P( 4 )=p5;
    P( 5 )=p6;
    P( 6 )=p7;
    P( 7 )=p8;
    P( 8 )=p9;
    P( 9 )=p10 ;
    BOOST_TEST_MESSAGE( "points created " );

    Polynomialset<10,16,6,10> Poly( m ) ;
    //  Points_type p(2); p(0)= 0.; p(1) = 0.;
    //  vectors_type P(1 ) ;
    //  P(0)=p ;
    BOOST_TEST_MESSAGE( "polyset created " );

    BOOST_TEST_MESSAGE(  "**************************************************Analytic method*********************************************" << "\n" );
    BOOST_TEST_MESSAGE( "coefficients= " << m << "\n" );
    BOOST_TEST_MESSAGE( "Points= " << P << "\n" );
    BOOST_TEST_MESSAGE( "evaluate_at_Points= " <<Poly.evaluates_Points ( P ) << "\n" );


    BOOST_TEST_MESSAGE(  "*****************************************Feel method***********************************************************" << "\n" );

    rt0_type::points_type pts( 2,10 );

    // pts(0 ,0) = 0. ; pts(1 ,0) = 0. ;

    pts( 0 ,0 ) = -1./3. ;
    pts( 1 ,0 ) = -1./3. ;
    pts( 0 ,1 ) = -1. ;
    pts( 1 ,1 ) = -1. ;
    pts( 0 ,2 ) =  1. ;
    pts( 1 ,2 ) = -1. ;
    pts( 0 ,3 ) = -1. ;
    pts( 1 ,3 ) =  1. ;
    pts( 0 ,4 ) = -1. ;
    pts( 1 ,4 ) =  0. ;
    pts( 0 ,5 )=-1./3.;
    pts( 1 ,5 ) =-1./3.;
    pts( 0 ,6 )= 0.;
    pts( 1 ,6 ) = 0.;
    pts( 0 ,7 )= 0.;
    pts( 1 ,7 ) =-1.;
    pts( 0 ,8 )=-1./2.;
    pts( 1 ,8 ) = 1./2.;
    pts( 0 ,9 )= 1./2.;
    pts( 1 ,9 ) =-1./2.;

    BOOST_TEST_MESSAGE(  "pts= " << pts << "\n" );
    auto eval_at_pts = rt0.evaluate( pts );
    BOOST_TEST_MESSAGE(  "eval_at_pts= "<< eval_at_pts << "\n" );
    double error = ublas::norm_frobenius( Poly.evaluates_Points( P )-eval_at_pts );
    BOOST_CHECK_SMALL( error, 1e-10 );


}

BOOST_AUTO_TEST_SUITE_END()


