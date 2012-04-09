/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
       Date: 2005-12-22

  Copyright (C) 2005,2006 EPFL

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
   \file test_qm.cpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2005-12-22
 */

/** The purpose of this file is to test the quadrature methods
    on the simplicies.
**/

// Boost.Test
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
using boost::unit_test::test_suite;

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <fstream>
#include <feel/feelpoly/im.hpp>

namespace mpl = boost::mpl;

/** Streams **/
//@{
std::ofstream PK_x_log( "PK_x_log.dat" );
std::ofstream PK_y_log( "PK_y_log.dat" );

std::ofstream ost_x( "Monom_x.txt" );
std::ofstream ost_y( "Monom_y.txt" );
//@}


double pow( double x, int i )
{
    double res = 1.0;

    for ( int l=1; l <= i; ++l )
        res *= x;

    return res;

}

template<Feel::uint16_type N, typename T>
void PK_Monom_N_opt()
{
    using namespace Feel;

    typedef T value_type;

    GaussLobatto<Simplex<2,1>, N, value_type> im;

    ublas::vector<value_type> x_i( ublas::row( im.points(),0 ) );
    ublas::vector<value_type> y_i( ublas::row( im.points(),1 ) );

    const value_type tol = value_type( 2.5 ) * type_traits<value_type>::epsilon();

    uint16_type Q = ( uint16_type )sqrt( im.nPoints() );

    value_type error_x = 0.0,error_y = 0.0;
    value_type res_x, res_y;
    uint16_type i_x=1,i_y=1;
    value_type sum_x=0.0,sum_y=0.0;

    ost_x << "Number of Points on the triangle : " << Q << "^2 = "<< im.nPoints() <<  std::endl;
    ost_y << "Number of Points on the triangle : " << Q << "^2 = "<< im.nPoints() <<  std::endl;

    do
    {
        if ( i_x%2 == 0 )
            res_x = value_type( 2.0 )/value_type( i_x+1.0 );

        else
            res_x = -value_type( 2.0 )/value_type( i_x+2.0 );

        sum_x = 0.0;

        for ( uint16_type l=0; l< x_i.size(); ++l )
            sum_x += pow( x_i( l ),i_x )*im.weight( l );

        error_x = math::abs( sum_x - res_x );

        ost_x << "i = " << i_x <<" Error = "<< error_x << "\n";
        ++i_x;
    }
    while ( error_x <= tol );

    if ( i_x-2 < 2*Q-3  )
        PK_x_log << "Q = " << Q << " ; i = " << i_x << " ; Error = " << error_x << std::endl;

    do
    {
        if ( i_y%2 == 0 )
            res_y = value_type( 2.0 )/value_type( i_y+1.0 );

        else
            res_y = value_type( 1.0 )/value_type( i_y+1.0 )*( ( 2.0/( value_type( i_y+2.0 ) ) ) - value_type( 2.0 ) );

        sum_y = 0.0;

        for ( uint16_type l=0; l< x_i.size(); ++l )
            sum_y += pow( y_i( l ),i_y )*im.weight( l );

        error_y = math::abs( sum_y - res_y );

        ost_y << "i = " << i_y <<" Error = "<< error_y << "\n";

        ++i_y;
    }
    while ( error_y <= tol );


    if ( i_y-2 < 2*Q-2  )
        PK_y_log << "Q = " << Q << " ; i = " << i_y << " ; Error = " << error_y << std::endl;

    BOOST_CHECK( ( i_x-2 >= 2*Q-3 ) && ( i_y-2 >= 2*Q-2 ) );
}

typedef boost::mpl::list<boost::mpl::int_<0>,
        boost::mpl::int_<1>,
        boost::mpl::int_<2>,
        boost::mpl::int_<10>,
        boost::mpl::int_<40>,
        boost::mpl::int_<80> > test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( PK_Monom_N_opt_double, T, test_types )
{
    PK_Monom_N_opt<T::value,double>();
}


#if 0
template<Feel::int16_type P>
void add_tests( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<P,double> ) ) );
#if defined( FEELPP_HAS_QD_REAL)
    //test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<P,dd_real> ) ) );
    //test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<P,qd_real> ) ) );
#endif

    add_tests<P+2>( test );
}

template<>
void add_tests<80>( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<80,double> ) ) );

#if defined( FEELPP_HAS_QD_REAL)
    //test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<80,dd_real> ) ) );
    //test->add( BOOST_TEST_CASE( ( PK_Monom_N_opt<80,qd_real> ) ) );
#endif
}

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{

    test_suite* test = BOOST_TEST_SUITE( "Integration methods on simplicies test suite" );

    PK_x_log << "This file attempt to save the errors encountered while numerically integrating $x^i$ on the triangle !" << "\n\n";
    PK_y_log << "This file attempt to save the errors encountered while numerically integrating $y^i$ on the triangle !" << "\n\n";

#if defined( FEELPP_HAS_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start( &old_cw );
#endif
    add_tests<2>( test );

    return test;
}
#endif
