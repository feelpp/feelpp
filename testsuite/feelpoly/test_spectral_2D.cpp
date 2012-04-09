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
   \file test_spectral_2D.cpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2006-04-26
 */

/** This file test the spectral accuracy of the Dubiner and boundary adapted basis in 2D on a single triangle **/


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
- \Lambda u + u = f
$$
**/

/** Streams **/

/** To save the results of the execution you must create the directory Data/2D/ **/

ofstream B_error( "Data/2D/Boundary_error_2D.dat" );
ofstream D_error( "Data/2D/Dubiner_error_2D.dat" );
ofstream B_error_L2( "Data/2D/Boundary_error_L2_2D.dat" );
ofstream D_error_L2( "Data/2D/Dubiner_error_L2_2D.dat" );
ofstream B_cond( "Data/2D/Boundary_cond_2D.dat" );
ofstream D_cond( "Data/2D/Dubiner_cond_2D.dat" );
ofstream B_timing( "Data/2D/Boundary_timing_2D.dat" );
ofstream D_timing( "Data/2D/Dubiner_timing_2D.dat" );
ofstream B_tol( "Data/2D/Boundary_tol_2D.dat" );
ofstream D_tol( "Data/2D/Dubiner_tol_2D.dat" );
ofstream B_Order( "Data/2D/Boundary_Order_2D.dat" );
ofstream D_Order( "Data/2D/Dubiner_Order_2D.dat" );


template<int_type N, typename T> void Boundary();
template<int_type N, typename T> void Ortho();

double pow( int_type N, int_type M )
{
    double res( 1.0 );

    for ( int_type i=1; i <= M ; ++i )
        res*=N;

    return res;
}

template<int_type P>
void add_tests( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Boundary<P,value_type> ) ) );
    test->add( BOOST_TEST_CASE( ( Ortho<P,value_type> ) ) );

    add_tests<P+2>( test );
}

#if !defined(NDEBUG) // In Debug mode, reducing the number of test cases

template<>
void add_tests<10>( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Boundary<10,value_type> ) ) );
    test->add( BOOST_TEST_CASE( ( Ortho<10,value_type> ) ) );
}

#else

template<>
void add_tests<16>( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Boundary<16,value_type> ) ) );
    test->add( BOOST_TEST_CASE( ( Ortho<16,value_type> ) ) );
}

#endif

/* Main test function */

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv[]*/ )
{
    test_suite* test = BOOST_TEST_SUITE( "2D spectral accuracy tests" );

#if defined( FEELPP_HAS_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start( &old_cw );
#endif
    add_tests<2>( test );

    return test;
}

// Problem given parameter

#define pi M_PI

template<typename T>
T f( typename node<T>::type const& t )
{
    return ( sin( pi*t[0] )*cos( pi*t[1] )*( 2*pi*pi+1 ) ) ;
}

template<typename T>
T g0( typename node<T>::type const& t )
{
    return ( pi*( cos( pi*t[0] )*cos( pi*t[1] ) - sin( pi*t[0] )*sin( pi*t[1] ) )/( sqrt( 2.0 ) ) );
}

template<typename T>
T g1( typename node<T>::type const& t )
{
    return ( -pi*( cos( pi*t[0] )*cos( pi*t[1] ) ) );
}

template<typename T>
T g2( typename node<T>::type const& t )
{
    return ( pi*( sin( pi*t[0] )*sin( pi*t[1] ) ) );
}

template<typename T>
T exact_sol( typename node<T>::type const& t )
{
    return ( sin( pi*t[0] )*cos( pi*t[1] ) );
}

template<typename T>
T dexact_sol_x( typename node<T>::type const& t )
{
    return ( pi*cos( pi*t[0] )*cos( pi*t[1] ) );
}

template<typename T>
T dexact_sol_y( typename node<T>::type const& t )
{
    return ( -pi*sin( pi*t[0] )*sin( pi*t[1] ) );
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
    ofstream str( "Data/Mat.m" );

    str << "function M = Mat()\n";
    str << "% Create the matrix Mat in order to manipulate it in Matlab\n";

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

/** Boundary adapted resolution **/


template<int_type N, typename T>
void Boundary()
{
    cout << "Basis = Boundary Adapted on Simplex !\n";
    cout << "N = " << N << endl;
    B_Order << N << endl;

    typedef T value_type;

    cout << "Construction of the Integration Method ... ";

    timer t_tot;
    timer t;

    const Gauss<Simplex<2,1>, N+N, value_type> Quad;
    cout << t.elapsed() << " seconds.\n";

    cout << "Basis Construction ... ";
    t.restart();
    const BoundaryAdaptedPolynomialSet<2, N, Scalar, value_type, Simplex> BA_Basis;
    timeout( t.elapsed() );

    const int_type Dim_Basis = BA_Basis.basis().coeff().size1();


    /** Matrix construction **/

    /** Stiffness **/
    cout << "\n/********************************/" << endl;
    cout << "Construction of Stiffness Matrix..."  <<endl;

    cout <<"Boundary Adapted Basis Derivation... ";
    t.restart();
    const ublas::vector<ublas::matrix<value_type> > Diff( BA_Basis.derivate( Quad.points() ) );
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

    ublas::matrix<value_type> Stiff( A + B );
    timeout( t.elapsed() );

    /** Mass **/

    cout << "\n/********************************/" << endl;
    cout << "Construction of the Mass Matrix..." << endl;

    cout <<"Boundary Adapted Basis Evaluation... ";
    t.restart();
    ublas::matrix<value_type> Psi = BA_Basis.evaluate( Quad.points() );
    timeout( t.elapsed() );

    cout << "Mass numeric construction ... ";
    t.restart();
    ublas::matrix<value_type> Mass ( ublas::prod( Psi, Wgts ) );
    Mass = ublas::prod( Mass, ublas::trans( Psi ) );
    timeout( t.elapsed() );

    if ( N == 14 )
    {
#if 0
        PrintMatlab<value_type>( Mass );
#else
        PrintMatlab<value_type>( Stiff );
#endif
    }

    /** Global **/

    cout << "\n/********************************/" << endl;
    cout << "Assembly of the Global matrix ..."<<endl;
    ublas::matrix<value_type> Global( Stiff + Mass );

    /** Right hand side **/

    cout << "Construction of the RHS ... ";
    t.restart();
    ublas::vector<value_type> F( Quad.points().size2() );
    ublas::vector<value_type> G0( Quad.fpoints( 0 ).size2() );
    ublas::vector<value_type> G2( Quad.fpoints( 2 ).size2() );
    ublas::diagonal_matrix<value_type> Wgts0( Quad.weights( 0 ).size() );
    ublas::diagonal_matrix<value_type> Wgts2( Quad.weights( 2 ).size() );

    for ( int_type i=0; i < F.size(); ++i )
    {
        F( i ) = f<value_type>( Quad.point( i ) );
    }

    F= ublas::prod( Wgts,F );
    F= ublas::prod( Psi,F );

    for ( int_type i=0; i < Wgts0.size1() ; ++i )
    {
        Wgts0( i,i ) = Quad.weights( 0 )( i );
        G0( i ) = g0<value_type>( ublas::column( Quad.fpoints( 0 ), i ) );
    }

    G0= ublas::prod( Wgts0,G0 );
    G0= ublas::prod( BA_Basis.evaluate( Quad.fpoints( 0 ) ),G0 );

    for ( int_type i=0; i < Wgts2.size1() ; ++i )
    {
        Wgts2( i,i ) = Quad.weights( 2 )( i );
        G2( i ) = g2<value_type>( ublas::column( Quad.fpoints( 2 ), i ) );
    }

    G2= ublas::prod( Wgts2,G2 );
    G2= ublas::prod( BA_Basis.evaluate( Quad.fpoints( 2 ) ),G2 );

    /** Global right-hand side assembly **/

    F+=G0+G2;
    timeout( t.elapsed() );


    /** Imposing Dirichlet boundary conditions **/

    cout << "Imposing Dirichlet boundary conditions ... ";
    t.restart();
    ublas::matrix<value_type> Global2( Global.size1()-( N-1 ), Global.size2()-( N-1 ) );
    ublas::vector<value_type> F2( F.size()- ( N-1 ) );

    for ( int_type i= 0 ; i < 3+( N-1 ) ; ++i )
    {
        for ( int_type j= 0 ; j < 3+( N-1 ) ; ++j )
        {
            Global2( i,j ) = Global( i,j );
        }

        for ( int_type j=3+2*( N-1 ); j < Global.size1() ; ++j )
        {
            Global2( i,j-( N-1 ) ) = Global( i,j );
        }

        F2( i ) = F( i );
    }

    for ( int_type i=3+2*( N-1 ); i < Global.size1() ; ++i )
    {
        for ( int_type j= 0 ; j < 3+( N-1 ) ; ++j )
        {
            Global2( i-N+1,j ) = Global( i,j );
        }

        for ( int_type j=3+2*( N-1 ); j < Global.size1() ; ++j )
        {
            Global2( i-N+1,j-( N-1 ) ) = Global( i,j );
        }

        F2( i-N+1 ) = F( i );
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

    ublas::vector<value_type> Sol2( Dim_Basis-N+1 );

    Sol2 = lu.solve( F2 );

    /** Expanding Solution on the Dirichlet boundary **/

    for ( int_type i= 0 ; i < 3+( N-1 ) ; ++i )
    {
        Sol( i ) = Sol2( i );
    }

    for ( int_type i= 3+( N-1 ) ; i < 3+2*( N-1 )  ; ++i )
    {
        Sol( i ) = value_type( 0.0 );
    }

    for ( int_type i=3+2*( N-1 ); i < Global.size1() ; ++i )
    {
        Sol( i ) = Sol2( i-N+1 );
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

    for ( uint16_type i=0; i< Error_H1.size(); ++i )
    {
        Approx( i )   = exact_sol<value_type>( ublas::column( Quad.points(), i ) );
        Approx_x( i ) = dexact_sol_x<value_type>( ublas::column( Quad.points(), i ) );
        Approx_y( i ) = dexact_sol_y<value_type>( ublas::column( Quad.points(), i ) );
    }

    ublas::vector<value_type> Delta( Approx - ublas::prod( ublas::trans( Psi ),Sol ) );
    ublas::vector<value_type> Delta_x( Approx_x - ublas::prod( ublas::trans( Diff( 0 ) ),Sol ) );
    ublas::vector<value_type> Delta_y( Approx_y - ublas::prod( ublas::trans( Diff( 1 ) ),Sol ) );


    Error_H1 =
        ublas::element_prod( Delta,Delta )
        + ublas::element_prod( Delta_x,Delta_x )
        + ublas::element_prod( Delta_y,Delta_y );


    value_type err_H1 =  sqrt( ublas::inner_prod( Error_H1, Quad.weights() ) );

    timeout( t.elapsed() );

    cout << "\nErreur H^1= " << err_H1 << endl;
    B_error << err_H1 << endl;

    double tol =  10*pow( 1.0/pow( N,N ), 0.2 );
    cout << "\n10*N^{-N/5} = " << tol << endl;

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


/** Dubiner resolution **/


template<int_type N, typename T>
void Ortho()
{
    cout << "Basis = Dubiner Orthonormal on Simplex !\n";
    cout << "N = " << N << endl;
    D_Order << N << endl;

    typedef T value_type;

    timer t_tot;
    timer t;

    timeout( t.elapsed() );

    cout << "Basis Construction ... ";
    t.restart();
    const OrthonormalPolynomialSet<2, N, Scalar, value_type, Simplex> Dubiner_Basis;
    //const Dubiner<2, N, Normalized<false>, value_type> Dubiner_Basis;
    timeout( t.elapsed() );

    const int_type Dim_Basis = Dubiner_Basis.basis().coeff().size1();

    cout << "Construction of the Integration Method ... \n";
    const Gauss<Simplex<2,1>, N+N, value_type> Quad;

    /** Matrix construction **/

    /** Stiffness **/
    cout << "\n/********************************/" << endl;
    cout << "Construction of Stiffness Matrix..."  <<endl;

    cout <<"Dubiner Basis Derivation... ";
    t.restart();
    const ublas::vector<ublas::matrix<value_type> > Diff( Dubiner_Basis.derivate( Quad.points() ) );
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

#if 1
    cout << "prod ... ";
    ublas::matrix<value_type> A ( ublas::prod( Diff( 0 ), Wgts ) );
    A = ublas::prod( A, ublas::trans( Diff( 0 ) ) );
    ublas::matrix<value_type> B ( ublas::prod( Diff( 1 ), Wgts ) );
    B = ublas::prod( B, ublas::trans( Diff( 1 ) ) );
#else
    cout << "axpy_prod ... ";
    ublas::matrix<value_type> A2( Diff( 0 ).size1(), Wgts.size2() );
    ublas::matrix<value_type> A( A2.size1(), Diff( 0 ).size1() );
    ublas::matrix<value_type> B2( Diff( 1 ).size1(), Wgts.size2() );
    ublas::matrix<value_type> B( B2.size1(), Diff( 1 ).size1() );

    ublas::matrix<value_type> tDiff_0( ublas::trans( Diff( 0 ) ) );
    ublas::matrix<value_type> tDiff_1( ublas::trans( Diff( 1 ) ) );

    ublas::axpy_prod( Diff( 0 ), Wgts, A2 );
    ublas::axpy_prod( A2, tDiff_0, A );
    ublas::axpy_prod( Diff( 1 ), Wgts, B2 );
    ublas::axpy_prod( B2, tDiff_1, B );
#endif

    ublas::matrix<value_type> Stiff( A + B );
    timeout( t.elapsed() );

    /** Mass **/

    cout << "\n/********************************/" << endl;
    cout << "Construction of the Mass Matrix..." << endl;

    cout <<"Dubiner Basis Evaluation... ";
    t.restart();
    ublas::matrix<value_type> Psi = Dubiner_Basis.evaluate( Quad.points() );
    timeout( t.elapsed() );

    cout << "Mass numeric construction ... ";
    t.restart();

#if 1
    cout << "prod ... ";
    ublas::matrix<value_type> Mass ( ublas::prod( Psi, Wgts ) );
    Mass = ublas::prod( Mass, ublas::trans( Psi ) );

#else
    cout << "axpy_prod ... ";
    ublas::matrix<value_type> M0( Psi.size1(), Wgts.size2() );
    ublas::matrix<value_type> Mass( Psi.size1(), Psi.size1() );
    ublas::matrix<value_type> tPsi( ublas::trans( Psi ) );

    ublas::axpy_prod( Psi, Wgts, M0 );
    ublas::axpy_prod( M0, tPsi, Mass );
#endif

    timeout( t.elapsed() );

    /** Global **/

    cout << "\n/********************************/" << endl;
    cout << "Assembly of the Global matrix ..."<<endl;
    ublas::matrix<value_type> Global( Stiff + Mass );

    /** Condition number estimation **/

    cout << "Condition number estimation ... ";
    t.restart();
    value_type condition = cond2<value_type>( Global );
    timeout( t.elapsed() );

    cout << "Condition number = " << condition << endl;
    D_cond << condition << endl;

    /** Right hand side **/

    cout << "Construction of the RHS ... ";
    t.restart();
    ublas::vector<value_type> F( Quad.points().size2() );
    ublas::vector<value_type> G0( Quad.fpoints( 0 ).size2() );
    ublas::vector<value_type> G1( Quad.fpoints( 1 ).size2() );
    ublas::vector<value_type> G2( Quad.fpoints( 2 ).size2() );

    ublas::diagonal_matrix<value_type> Wgts0( Quad.weights( 0 ).size() );
    ublas::diagonal_matrix<value_type> Wgts1( Quad.weights( 1 ).size() );
    ublas::diagonal_matrix<value_type> Wgts2( Quad.weights( 2 ).size() );

    for ( int_type i=0; i < F.size(); ++i )
    {
        F( i ) = f<value_type>( Quad.point( i ) );
    }

    F= ublas::prod( Wgts,F );
    F= ublas::prod( Psi,F );

    for ( int_type i=0; i < Wgts0.size1() ; ++i )
    {
        Wgts0( i,i ) = Quad.weights( 0 )( i );
        G0( i ) = g0<value_type>( ublas::column( Quad.fpoints( 0 ), i ) );
    }

    G0= ublas::prod( Wgts0,G0 );
    G0= ublas::prod( Dubiner_Basis.evaluate( Quad.fpoints( 0 ) ),G0 );

    for ( int_type i=0; i < Wgts1.size1() ; ++i )
    {
        Wgts1( i,i ) = Quad.weights( 1 )( i );
        G1( i ) = g1<value_type>( ublas::column( Quad.fpoints( 1 ), i ) );
    }

    G1= ublas::prod( Wgts1,G1 );
    G1= ublas::prod( Dubiner_Basis.evaluate( Quad.fpoints( 1 ) ),G1 );

    for ( int_type i=0; i < Wgts2.size1() ; ++i )
    {
        Wgts2( i,i ) = Quad.weights( 2 )( i );
        G2( i ) = g2<value_type>( ublas::column( Quad.fpoints( 2 ), i ) );
    }

    G2= ublas::prod( Wgts2,G2 );
    G2= ublas::prod( Dubiner_Basis.evaluate( Quad.fpoints( 2 ) ),G2 );

    /** Global right-hand side assembly **/

    F+=G0+G1+G2;
    timeout( t.elapsed() );

    /** System Resolution **/

    cout << "System Resolution ... ";
    t.restart();
    ublas::vector<value_type> Sol( Dim_Basis );
    LU<ublas::matrix<value_type> > lu( Global );

    Sol = lu.solve( F );

    timeout( t.elapsed() );


    /** L_2 Error Estimation **/

    cout << "L_2 Error estimation ... ";
    t.restart();
    ublas::vector<value_type> Error( Quad.points().size2() );

    for ( int_type i=0; i< Error.size(); ++i )
        Error( i ) = exact_sol<value_type>( ublas::column( Quad.points(), i ) );

    Error = Error - ublas::prod( ublas::trans( Psi ),Sol );
    Error = ublas::element_prod( Error,Error );

    value_type err =  sqrt( ublas::inner_prod( Error, Quad.weights() ) );
    timeout( t.elapsed() );

    cout << "\nErreur L_2= " << err << endl;

    D_error_L2 << err << endl;


    /** H¹ Error Estimation **/

    cout << "H^1 Error estimation ... ";
    t.restart();

    ublas::vector<value_type> Error_H1( Quad.points().size2() );

    ublas::vector<value_type> Approx ( Error_H1.size() );
    ublas::vector<value_type> Approx_x( Error_H1.size() );
    ublas::vector<value_type> Approx_y( Error_H1.size() );

    for ( uint16_type i=0; i< Error_H1.size(); ++i )
    {
        Approx( i )   = exact_sol<value_type>( ublas::column( Quad.points(), i ) );
        Approx_x( i ) = dexact_sol_x<value_type>( ublas::column( Quad.points(), i ) );
        Approx_y( i ) = dexact_sol_y<value_type>( ublas::column( Quad.points(), i ) );
    }

    ublas::vector<value_type> Delta( Approx - ublas::prod( ublas::trans( Psi ),Sol ) );
    ublas::vector<value_type> Delta_x( Approx_x - ublas::prod( ublas::trans( Diff( 0 ) ),Sol ) );
    ublas::vector<value_type> Delta_y( Approx_y - ublas::prod( ublas::trans( Diff( 1 ) ),Sol ) );


    Error_H1 =
        ublas::element_prod( Delta,Delta )
        + ublas::element_prod( Delta_x,Delta_x )
        + ublas::element_prod( Delta_y,Delta_y );


    value_type err_H1 =  sqrt( ublas::inner_prod( Error_H1, Quad.weights() ) );

    timeout( t.elapsed() );

    cout << "\nErreur H^1= " << err_H1 << endl;
    D_error << err_H1 << endl;

    double tol =  10*pow( 1.0/pow( N,N ), 0.2 );
    cout << "\n10*N^{-N/5} = " << tol << endl;

    D_tol << tol << endl;

    cout << "\n-------------------------------\n";
    cout << "Total runtime = ";
    double timing( t_tot.elapsed() );

    timeout( timing );
    cout << "-------------------------------\n";
    cout << "\n\n";
    D_timing << timing << endl;

    BOOST_CHECK( ( err_H1 <= tol ) );
}
