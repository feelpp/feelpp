/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
       Date: 2005-12-22

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2009-2010 Université de Grenoble 1 (Joseph Fourier)

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
   \file test_quad_order.cpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2005-12-22
 */

/**
    The purpose of this file is to test the quadrature methods
    on tensorised geometries.
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
std::ofstream QK_log( "QK_error.log" );
std::ofstream ost( "QK_results.txt" );

std::ofstream QK_log_3D( "QK_error_3D.log" );
std::ofstream ost_3D( "QK_results_3D.txt" );

std::ofstream PK_log_g( "PK_log_g.dat" );
std::ofstream ost_g( "PK_results_g.dat" );

std::ofstream PK_x_log( "PK_x_log.dat" );
std::ofstream PK_y_log( "PK_y_log.dat" );

std::ofstream ost_x( "Monom_x.txt" );
std::ofstream ost_y( "Monom_y.txt" );

std::ofstream PK_x_log_3D( "PK_x_log_3D.log" );
std::ofstream ost_x_3D( "Monom_x_3D.txt" );

//@}


double pow( double x, int i )
{
    double res = 1.0;

    for ( int l=1; l <= i; ++l )
        res *= x;

    return res;

}

template<Feel::uint16_type N, typename T>
void QK_find_N_opt()
{
    using namespace Feel;
    typedef T value_type;

    Gauss<Hypercube<2,1>, N ,value_type > im;

    ublas::vector<value_type> x_i( ublas::row( im.points(),0 ) );
    ublas::vector<value_type> y_i( ublas::row( im.points(),1 ) );
    const value_type tol = value_type( 7.0 )*type_traits<value_type>::epsilon();

    /* 1D number of Nodes */
    //@{
    uint16_type Q = ( uint16_type )sqrt( im.nPoints() );
    //@}

    value_type error = 0.0;
    value_type res;
    uint16_type i=1;
    value_type sum=0.0;

    ost << "Nbre of Points on the Quadrangle : " << Q << "^2 = "<< im.nPoints() <<  std::endl;

    do
    {
        if ( i%2 == 0 )
            res = value_type( 4.0 )/value_type( i+1.0 )/value_type( i+1.0 );

        else
            res = 0.0;

        sum = 0.0;

        for ( uint16_type l=0; l< x_i.size(); ++l )
        {
            sum += pow( x_i( l ),i )*pow( y_i( l ),i )*im.weight( l );
        }

        error = math::abs( sum - res );
        ost << "i = " << i <<" Error = "<< error << "\n";
        ++i;
    }
    while ( error <= tol );

    if ( i-2 < 2*Q-1  )
        QK_log << "Q = " << Q << " ; i = " << i << " ; Error" << error << std::endl;

    BOOST_CHECK( i-2 >= 2*Q-1  );
}

template<Feel::uint16_type N, typename T>
void QK_find_N_opt_3D()
{
    using namespace Feel;
    typedef T value_type;

    Gauss<Hypercube<3,1>, N ,value_type > im;

    ublas::vector<value_type> x_i( ublas::row( im.points(),0 ) );
    ublas::vector<value_type> y_i( ublas::row( im.points(),1 ) );
    ublas::vector<value_type> z_i( ublas::row( im.points(),2 ) );

    const value_type tol = value_type( 2.5 )*type_traits<value_type>::epsilon();

    /* 1D number of Nodes */
    //@{
    uint16_type Q = ( uint16_type )pow( im.nPoints(),double( 1 )/3 );
    //@}

    value_type error = 0.0;
    value_type res;
    uint16_type i=1;
    value_type sum=0.0;

    ost << "Nbre of Points on the Hexaedra : " << Q << "^3 = "<< im.nPoints() <<  std::endl;

    do
    {
        if ( i%2 == 0 )
            res = value_type( 8.0 )/value_type( i+1.0 )/value_type( i+1.0 )/value_type( i+1.0 );

        else
            res = 0.0;

        sum = 0.0;

        for ( uint16_type l=0; l< x_i.size(); ++l )
        {
            sum += pow( x_i( l ), i )*pow( y_i( l ),i )*pow( z_i( l ),i )*im.weight( l );
        }

        error = math::abs( sum - res );
        ost << "i = " << i <<" Error = "<< error << "\n";
        ++i;
    }
    while ( error <= tol );

    if ( i-2 < 2*Q-1  )
        QK_log_3D << "Q = " << Q << " ; i = " << i << " ; Error" << error << std::endl;

    BOOST_CHECK( i-2 >= 2*Q-1  );
}

/**
 * PK 2D
 *
 */
template<Feel::uint16_type N, typename T>
void PK_Monom_N_opt_gauss_lobatto()
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
template<Feel::uint16_type N, typename T>
void PK_find_N_opt_gauss()
//@}
{
    using namespace Feel;

    typedef T value_type;

    Gauss<Simplex<2,1> , N, value_type> im;

    ublas::vector<value_type> x_i( ublas::row( im.points(),0 ) );
    ublas::vector<value_type> y_i( ublas::row( im.points(),1 ) );

    const value_type tol = value_type( 7.0 )*type_traits<value_type>::epsilon();

    ost_g << "Tolerance = " << tol << "\n";

    uint16_type Q = ( uint16_type )sqrt( im.nPoints() );

    value_type error;
    value_type res;
    uint16_type i=1;
    value_type sum;

    ost_g << "Nbre of Points on the triangle : " << Q << "^2 = "<< im.nPoints() <<  std::endl;

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
        ost_g << "i = " << i <<" Error = "<< error << "\n";

        ++i;
    }
    while ( error <= tol );

    if ( i-2 < Q-1  )
        PK_log_g << "Q = " << Q << " ; i = " << i << " ; Error = " << error << std::endl;

    BOOST_CHECK( i-2 >= Q-1  );

}

/**
 * PK 3D
 */
template<Feel::uint16_type N, typename T>
void PK_Monom_N_opt_3D()
{
    int Q = ( N+3 )/2+1;
    BOOST_TEST_MESSAGE( "Check that the quadrature order in Tetrahedron for N=" << N << " is at least equal to " << 2*Q-3 );
    using namespace Feel;

    typedef T value_type;

    GaussLobatto<Simplex<3,1>, N, value_type> im;

    ublas::vector<value_type> x_i( ublas::row( im.points(),0 ) );

    const value_type tol = value_type( 2.5 ) * type_traits<value_type>::epsilon();



    value_type error_x = 0.0;
    value_type res_x;
    int i_x=1;
    value_type sum_x=0.0;

    ost_x_3D << "Number of Points on the tetraedra : " << Q << "^3 = "<< im.nPoints() <<  std::endl;

    do
    {
        if ( i_x%2 == 0 )
            res_x = value_type( 1.0 )/value_type( i_x+3.0 ) + value_type( 1.0 )/value_type( i_x+1.0 );

        else
            res_x = -value_type( 2.0 )/value_type( i_x+2.0 );

        sum_x = 0.0;

        for ( uint16_type l=0; l< x_i.size(); ++l )
            sum_x += pow( x_i( l ),i_x )*im.weight( l );

        error_x = math::abs( sum_x - res_x );

        //  std::cout << "i = " << i_x <<"\nError = "<< error_x << "\nRes = "<< res_x << " - sum_x = "<< sum_x <<"\n";
        ost_x_3D << "i = " << i_x <<" Error = "<< error_x << "\n";
        ++i_x;
    }
    while ( error_x <= tol );

    if ( i_x-2 < 2*Q-3  )
        PK_x_log_3D << "Q = " << Q << " ; i = " << i_x << " ; Error = " << error_x << std::endl;

    BOOST_TEST_MESSAGE( "check for order " << N
                        << " that " << i_x-2
                        << " is greater than " << 2*Q-3
                        << " : error = " << error_x
                        << " with tolerance = " << tol );
    BOOST_CHECK( ( i_x-2 >= 2*Q-3 ) );
}

/*
 * Testsuite
 */
typedef boost::mpl::list<boost::mpl::int_<2>,
        boost::mpl::int_<4>,
        boost::mpl::int_<8>,
        boost::mpl::int_<16>,
        boost::mpl::int_<32>,
        boost::mpl::int_<40> > test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_QK_find_N_opt_double, T, test_types )
{
    QK_find_N_opt<T::value,double>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( test_QK_find_N_opt_3D_double, T, test_types )
{
    QK_find_N_opt_3D<T::value,double>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( PK_Monom_N_opt_gauss_lobatto_double, T, test_types )
{
    PK_Monom_N_opt_gauss_lobatto<T::value,double>();
}
BOOST_AUTO_TEST_CASE_TEMPLATE( PK_Monom_N_opt_gauss_double, T, test_types )
{
    PK_find_N_opt_gauss<T::value,double>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( PK_Monom_N_opt_3D_double, T, test_types )
{
    PK_Monom_N_opt_3D<T::value,double>();
}
