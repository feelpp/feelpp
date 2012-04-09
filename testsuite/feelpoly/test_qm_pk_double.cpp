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
using boost::unit_test::test_suite;

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <fstream>
#include <feel/feelpoly/im.hpp>

namespace mpl = boost::mpl;

/** Streams **/
//@{
std::ofstream PK_log( "PK_log.dat" );
std::ofstream ost( "PK_double_results.txt" );
//@}

double pow( double x, int i )
{
    double res = 1.0;

    for ( int l=1; l <= i; ++l )
        res *= x;

    return res;

}

/** Recursive test on 2D - Simplex **/
//@{
template<Feel::uint16_type N, typename T>
void PK_find_N_opt()
//@}
{
    using namespace Feel;

    typedef T value_type;

    Gauss<Simplex<2,1> , N, value_type> im;

    ublas::vector<value_type> x_i( ublas::row( im.points(),0 ) );
    ublas::vector<value_type> y_i( ublas::row( im.points(),1 ) );

    const value_type tol = value_type( 7.0 )*type_traits<value_type>::epsilon();

    ost << "Tolerance = " << tol << "\n";

    uint16_type Q = ( uint16_type )sqrt( im.nPoints() );

    value_type error;
    value_type res;
    uint16_type i=1;
    value_type sum;

    ost << "Nbre of Points on the triangle : " << Q << "^2 = "<< im.nPoints() <<  std::endl;

    do
    {
        if ( i%2 == 0 )
            res = value_type( 2.0 )/value_type( double( i )+1.0 )/value_type( double( i )+1.0 );

        else
            res = value_type( 0.0 );

        sum = value_type( 0.0 );

        for ( uint16_type l=0; l< x_i.size(); ++l )
        {
            value_type p = pow( x_i( l ),i )*pow( y_i( l ),i )*im.weight( l );
            //        std::cout << "Contrib = " << p << "\n";
            sum += p;
        }

#if 0
        std::cout << "i = " << i << "\n";
        std::cout << "res = " << res << "\n";
        std::cout << "sum = " << sum << "\n";
#endif

        error = math::abs( sum - res );
        ost << "i = " << i <<" Error = "<< error << "\n";

        ++i;
    }
    while ( error <= tol );

    if ( i-2 < Q-1  )
        PK_log << "Q = " << Q << " ; i = " << i << " ; Error = " << error << std::endl;

    BOOST_CHECK( i-2 >= Q-1  );

}

template<Feel::int16_type P>
void add_tests( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    test->add( BOOST_TEST_CASE( ( PK_find_N_opt<P,double> ) ) );

#if defined( FEELPP_HAS_QD_REAL)
    test->add( BOOST_TEST_CASE( ( PK_find_N_opt<P,dd_real> ) ) );
    test->add( BOOST_TEST_CASE( ( PK_find_N_opt<P,qd_real> ) ) );
#endif

    add_tests<P+2>( test );
}

template<>
void add_tests<40>( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    test->add( BOOST_TEST_CASE( ( PK_find_N_opt<40,double> ) ) );

#if defined( FEELPP_HAS_QD_REAL)
    test->add( BOOST_TEST_CASE( ( PK_find_N_opt<40,dd_real> ) ) );
    test->add( BOOST_TEST_CASE( ( PK_find_N_opt<40,qd_real> ) ) );
#endif
}

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{

    test_suite* test = BOOST_TEST_SUITE( "Integration methods on simplicies test suite" );

    PK_log << "This file attempt to save the errors encountered while numerically integrating $x^iy^i$ on the triangle !" << "\n\n";
#if defined( FEELPP_HAS_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start( &old_cw );
#endif
    add_tests<20>( test );

    return test;
}
