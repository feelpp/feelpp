/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
       Date: 2006-04-26

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
   \file spectral_3D.cpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2006-04-26
 */

/** This file test the spectral accuracy of the boundary adapted basis in 3D on a single tetrahedra **/

/** Headers **/


// Boost Test

#include <boost/test/unit_test.hpp>

// Boost timer

#include <boost/timer.hpp>

// Boost numeric

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

// Gmm

#include<gmm.h>

// Feel
#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/quadpoint.hpp>
#include <feel/feelpoly/polynomialset.hpp>

using boost::unit_test::test_suite;
using boost::timer;

using namespace std;
using namespace Feel;

//Typedefs

typedef uint16_type int_type;

/**
We solve the Poisson problem :
$$
- \Lambda u = f
$$
**/

/** Streams **/

ofstream B_error( "Data/3D/Boundary_error_H1_3D.dat" );
ofstream B_error_L2( "Data/3D/Boundary_error_L2_3D.dat" );
ofstream B_cond( "Data/3D/Boundary_cond_3D.dat" );
ofstream B_timing( "Data/3D/Boundary_timing_3D.dat" );
ofstream B_tol( "Data/3D/Boundary_tol_3D.dat" );
ofstream B_Order( "Data/3D/Boundary_Order_3D.dat" );

template<int_type N, typename T> void Simp();

template<int_type P>
void add_tests( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Simp<P,value_type> ) ) );

    add_tests<P+2>( test );
}

#if !defined(NDEBUG) // In Debug mode, reducing the number of test cases

template<>
void add_tests<5>( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Simp<5,value_type> ) ) );
}

#else

template<>
void add_tests<9>( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Simp<9,value_type> ) ) );
}
#endif

double pow( int_type N, int_type M )
{
    double res( 1.0 );

    for ( int_type i=1; i <= M ; ++i )
        res*=N;

    return res;
}

/* Main test function */

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv[]*/ )
{
    test_suite* test = BOOST_TEST_SUITE( "3D Spectral Accuracy tests" );

#if defined( FEELPP_HAS_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start( &old_cw );
#endif

    add_tests<1>( test );
    return test;
}

// Problem given parameter

#define cste (-1.0)

template<typename T>
T f( typename node<T>::type const& t )
{
    return ( -exp( cste*( t[0]+t[1]+t[2] ) )*
             (
                 T( 2.0 )*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )
                 + T( 2.0 )*( T( 1.0 )+t[0] )*( T( 1.0 )+t[2] )
                 + T( 2.0 )*( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )

                 + T( 2.0 )*cste*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )*( t[0]+t[1]+t[2]+T( 1.0 ) )
                 + T( 2.0 )*cste*( T( 1.0 )+t[0] )*( T( 1.0 )+t[2] )*( t[0]+t[1]+t[2]+T( 1.0 ) )
                 + T( 2.0 )*cste*( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( t[0]+t[1]+t[2]+T( 1.0 ) )

                 + T( 6.0 )*cste*( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )
                 + T( 3.0 )*cste*cste*( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )*( t[0]+t[1]+t[2]+T( 1.0 ) )
             ) );
}

template<typename T>
T exact_sol( typename node<T>::type const& t )
{
    return ( ( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )*( t[0]+t[1]+t[2]+T( 1.0 ) )*exp( cste*( t[0]+t[1]+t[2] ) ) );
}

template<typename T>
T dexact_sol_x( typename node<T>::type const& t )
{
    return ( exp( cste*( t[0]+t[1]+t[2] ) )*
             (
                 ( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )*( t[0]+t[1]+t[2]+T( 1.0 ) )

                 +( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )

                 + cste*( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )*( t[0]+t[1]+t[2]+T( 1.0 ) )
             ) );
}

template<typename T>
T dexact_sol_y( typename node<T>::type const& t )
{
    return ( exp( cste*( t[0]+t[1]+t[2] ) )*
             (
                 ( T( 1.0 )+t[0] )*( T( 1.0 )+t[2] )*( t[0]+t[1]+t[2]+T( 1.0 ) )

                 +( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )

                 + cste*( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )*( t[0]+t[1]+t[2]+T( 1.0 ) )
             ) );
}

template<typename T>
T dexact_sol_z( typename node<T>::type const& t )
{
    return ( exp( cste*( t[0]+t[1]+t[2] ) )*
             (
                 ( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( t[0]+t[1]+t[2]+T( 1.0 ) )

                 + ( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )

                 + cste*( T( 1.0 )+t[0] )*( T( 1.0 )+t[1] )*( T( 1.0 )+t[2] )*( t[0]+t[1]+t[2]+T( 1.0 ) )
             ) );
}

/** Compute the spectral condition number using a QR algorithm **/

template<typename T>
T cond2( ublas::matrix<T> const& M )
{
    gmm::dense_matrix<T> Q( M.size1(), M.size2() );

    for ( int_type i=0; i < M.size1(); ++i )
        for ( int_type j=0; j < M.size2(); ++j )
        {
            Q( i,j ) = M( i,j );
        }

    return gmm::condition_number( Q );
}

template<typename T>
void PrintMatlab( const ublas::matrix<T> & M )
{
    ofstream str( "Data/Poisson3D.m" );

    str << "function M = Poisson3D()\n";
    str << "% Create the matrix M in order to manipulate it in Matlab\n";

    str << "M = [ ... \n";

    for ( int_type i=0; i < M.size1() ; ++i )
    {
        for ( int_type j=0; j < M.size2() ; ++j )
            str << M( i,j ) << " ";

        str << "\n";
    }

    str <<"];\n";

    str <<"\n return";

}

void timeout( double time )
{
    if ( time < 60 )
        cout << time << " seconds.\n";

    else if ( time < 3600 )
    {
        int_type min = ( int_type )floor( time )/60;
        cout << min << " minutes " << ( time-60*min ) << " seconds.\n";
    }

    else
    {
        int_type hour = ( int_type )floor( time )/3600;
        int_type min = ( int_type )floor( time-hour*3600 )/60;
        cout << hour << " hours " << min << " minutes " << ( time-3600*hour-60*min ) << " seconds.\n";
    }
}


/** All functions return the error **/

/** Simplex resolution **/


template<int_type N, typename T>
void Simp()
{
    cout << "N = " << N << endl;
    B_Order << N << endl;

    typedef T value_type;

    cout << "Construction of the Integration Method ... ";

    timer t_tot;
    timer t;

    const Gauss<Simplex<3,1>, N+N, value_type> Quad;
    cout << t.elapsed() << " seconds.\n";

    cout << "Basis Construction ... ";
    t.restart();
    const BoundaryAdaptedPolynomialSet<3, N, Scalar, value_type, Simplex> BA_Basis;
    timeout( t.elapsed() );

    const int_type Dim_Basis = BA_Basis.basis().coeff().size1();

    /** Matrix construction **/

    /** Stiffness **/
    cout << "\n/********************************/" << endl;
    cout << "Construction of Stiffness Matrix..."  <<endl;

    cout <<"Boundary Adapted Basis Derivation... ";
    t.restart();
    ublas::vector<ublas::matrix<value_type> > Diff( BA_Basis.derivate( Quad.points() ) );
    timeout( t.elapsed() );


    cout << "Wgts construction ... ";
    t.restart();
    ublas::diagonal_matrix<value_type> Wgts( Quad.weights().size() );

    for ( int_type i=0; i < Wgts.size1() ; ++i )
    {
        Wgts( i,i ) = Quad.weights()( i );
    }

    timeout( t.elapsed() );

    cout << "Stiffness numeric construction ... ";

    t.restart();
    ublas::matrix<value_type> A ( ublas::prod( Diff( 0 ), Wgts ) );
    A = ublas::prod( A, ublas::trans( Diff( 0 ) ) );
    ublas::matrix<value_type> B ( ublas::prod( Diff( 1 ), Wgts ) );
    B = ublas::prod( B, ublas::trans( Diff( 1 ) ) );
    ublas::matrix<value_type> C ( ublas::prod( Diff( 2 ), Wgts ) );
    C = ublas::prod( C, ublas::trans( Diff( 2 ) ) );

    ublas::matrix<value_type> Global( A + B + C );

    glas::clean<ublas::matrix<value_type>  >( Global,1e-15 );

    if ( N==10 )
        PrintMatlab<value_type>( Global );

    timeout( t.elapsed() );

    cout <<"Boundary Adapted Basis Evaluation... ";
    t.restart();
    ublas::matrix<value_type> Psi = BA_Basis.evaluate( Quad.points() );
    timeout( t.elapsed() );


    /** Right hand side **/

    cout << "Construction of the RHS ... ";
    t.restart();
    ublas::vector<value_type> F( Quad.points().size2() );

    for ( int_type i=0; i < F.size(); ++i )
    {
        F( i ) = f<value_type>( Quad.point( i ) );
    }

    F= ublas::prod( Wgts,F );
    F= ublas::prod( Psi,F );

    timeout( t.elapsed() );


    /** Imposing Dirichlet boundary conditions **/

    cout << "Imposing Dirichlet boundary conditions ... ";
    t.restart();
    const int_type Dir_Dim = ( N-1 )*( N-2 )*( N-3 )/6;
    const int_type step = 4+6*( N-1 )+2*( N-1 )*( N-2 );

    ublas::matrix<value_type> Global2( Dir_Dim, Dir_Dim );
    ublas::vector<value_type> F2( Dir_Dim );

    for ( int_type i= 0 ; i < Dir_Dim ; ++i )
    {
        for ( int_type j= 0 ; j < Dir_Dim ; ++j )
        {
            Global2( i,j ) = Global( i+step,j+step );
        }

        F2( i ) = F( i+step );
    }

    timeout( t.elapsed() );

    /** Condition number estimation **/

    cout << "Condition number estimation ... ";
    t.restart();
    value_type condition = cond2<value_type>( Global2 );
    timeout( t.elapsed() );

    cout << "Condition number = " << condition << endl;
    B_cond << condition << endl;

    /** System Resolution **/

    cout << "System Resolution ... ";
    t.restart();
    ublas::vector<value_type> Sol( Dim_Basis );
    LU<ublas::matrix<value_type> > lu( Global2 );

    ublas::vector<value_type> Sol2( Dir_Dim );

    Sol2 = lu.solve( F2 );

    /** Expanding Solution on the Dirichlet boundary **/

    for ( int_type i= 0 ; i < step ; ++i )
    {
        Sol( i ) = value_type( 0.0 );
    }

    for ( int_type i= step  ; i < Dir_Dim+step  ; ++i )
    {
        Sol( i ) = Sol2( i-step );
    }

    timeout( t.elapsed() );

    /** Error Estimation **/

    cout << "Error estimation ... ";
    t.restart();
    ublas::vector<value_type> Error( Quad.points().size2() );

    for ( int_type i=0; i< Error.size(); ++i )
        Error( i ) = exact_sol<value_type>( ublas::column( Quad.points(), i ) );

    Error = Error - ublas::prod( ublas::trans( Psi ),Sol );
    Error = ublas::element_prod( Error,Error );

    value_type err =  sqrt( ublas::inner_prod( Error, Quad.weights() ) );
    timeout( t.elapsed() );

    cout << "\nErreur L_2 = " << err << endl;
    B_error_L2 << err << endl;

    /** H¹ Error Estimation **/

    cout << "H^1 Error estimation ... ";
    t.restart();
    ublas::vector<value_type> Error_H1( Quad.points().size2() );

    ublas::vector<value_type> Approx ( Error_H1.size() );
    ublas::vector<value_type> Approx_x( Error_H1.size() );
    ublas::vector<value_type> Approx_y( Error_H1.size() );
    ublas::vector<value_type> Approx_z( Error_H1.size() );

    for ( uint16_type i=0; i< Error_H1.size(); ++i )
    {
        Approx( i )   = exact_sol<value_type>( ublas::column( Quad.points(), i ) );
        Approx_x( i ) = dexact_sol_x<value_type>( ublas::column( Quad.points(), i ) );
        Approx_y( i ) = dexact_sol_y<value_type>( ublas::column( Quad.points(), i ) );
        Approx_z( i ) = dexact_sol_z<value_type>( ublas::column( Quad.points(), i ) );
    }

    ublas::vector<value_type> Delta( Approx - ublas::prod( ublas::trans( Psi ),Sol ) );
    ublas::vector<value_type> Delta_x( Approx_x - ublas::prod( ublas::trans( Diff( 0 ) ),Sol ) );
    ublas::vector<value_type> Delta_y( Approx_y - ublas::prod( ublas::trans( Diff( 1 ) ),Sol ) );
    ublas::vector<value_type> Delta_z( Approx_z - ublas::prod( ublas::trans( Diff( 2 ) ),Sol ) );


    Error_H1 =
        ublas::element_prod( Delta,Delta )
        + ublas::element_prod( Delta_x,Delta_x )
        + ublas::element_prod( Delta_y,Delta_y )
        + ublas::element_prod( Delta_z,Delta_z );

    value_type err_H1 =  sqrt( ublas::inner_prod( Error_H1, Quad.weights() ) );
    timeout( t.elapsed() );

    cout << "\nErreur H^1= " << err_H1 << endl;
    B_error << err_H1 << endl;

    double tol =  10*pow( 1.0/pow( N,N ), 0.6 );
    cout << "\n10*N^{-N*0.6} = " << tol << endl;

    B_tol << tol << endl;
    cout << "\n-------------------------------\n";
    cout << "Total runtime = ";
    double timing( t_tot.elapsed() );
    timeout( timing );
    cout << "-------------------------------\n";
    cout << "\n\n";
    B_timing << timing << endl;

    BOOST_CHECK( ( err_H1 <= tol ) );
}
