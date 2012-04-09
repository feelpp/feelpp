/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
       Date: 2006-04-27

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
   \file test_spectral_convergence_2D.cpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2006-04-27
 */

/**
 * This file test the spectral accuracy of the Boundary adapted basis in 2D.
 * The 2D-geometry is decomposed into many triangles (h-discretisation).
 * We observe the spectral accuracy (p-reffinement) on different meshes.
 * You can pass the meshSize at run time, for example (h=0.5) :
 *
 * $./test_spectral_convergence_2D 0.5
 *
 * The default value is \f$ h=0.5 \f$.
**/

#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/rational.hpp>

// Boost timer

#include <boost/timer.hpp>

// Boost numeric

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

// Gmm

#include<gmm.h>

#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelalg/matrixublas.hpp>
#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/matrixgmm.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <gmm_iter_solvers.h>

using boost::unit_test::test_suite;
using boost::timer;

using namespace std;
using namespace Feel;

/** Streams **/

/** To save the results of the execution you must create the directory Data/2D/ **/

ofstream S_error( "Data/2D/error_H1_2D.dat" );
ofstream S_error_L2( "Data/2D/error_L2_2D.dat" );
ofstream S_timing( "Data/2D/timing_2D.dat" );
ofstream S_dofs( "Data/2D/dofs_2D.dat" );
ofstream S_iter( "Data/2D/solver_iter_2D.dat" );
ofstream S_order( "Data/2D/polynomial_order_2D.dat" );
ofstream S_tol( "Data/2D/tolerance_2D.dat" );

typedef int16_type int_type;

/**
 * We solve the Poisson problem :
 * \f$  - \Lambda u + u = f  \f$
 * with Neumann boundary conditions.
 **/

template<int_type N, typename T> void Poisson();

double hsize( 0.5 );

void timeout( double time );

template<int_type P>
void add_tests( test_suite* test )
{
    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Poisson<P,value_type> ) ) );

    add_tests<P+2>( test );
}

#if !defined(NDEBUG) // In Debug mode, reducing the number of test cases

#define ORDER 6
#else

#define ORDER 10
#endif

template<>
void add_tests<ORDER>( test_suite* test )
{
    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Poisson<ORDER,value_type> ) ) );
}


/* Main test function */

test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
    test_suite* test = BOOST_TEST_SUITE( "Spectral Accuracy (BA) test using Feel Language" );

    if ( argc > 1 )
    {
        hsize = static_cast<double> ( atof( argv[1] ) );
    }

#if defined( FEELPP_HAS_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start( &old_cw );
#endif

    add_tests<2>( test );
    return test;
}

double pow( int_type N, int_type M )
{
    double res( 1.0 );

    for ( int_type i=1; i <= M ; ++i )
        res*=N;

    return res;
}

template<int_type N, typename T>
void Poisson()
{
    cout << "\n\nN = " << N << endl;
    S_order << N << endl;

    typedef T value_type;

    timer t_tot;
    timer t;

    double meshSize = hsize;
    cout << "hsize = " << meshSize << endl;

    Gmsh __gmsh;
    string fname;
    ostringstream ostr;
    ostringstream nameStr;

    ostr << "h=" << meshSize << ";\n"
         << "Point(1) = {-1, -1,0.0,h};\n"
         << "Point(2) = { 1, -1,0.0,h};\n"
         << "Point(3) = {-1,  1,0.0,h};\n"
         << "Line(1) = {2,3};\n"
         << "Line(2) = {3,1};\n"
         << "Line(3) = {1,2};\n"
         << "Line Loop(4) = {1,2,3};\n"
         << "Plane Surface(5) = {4};\n"
         << "Physical Surface(30) = {5};\n"
         << "Physical Line(31) = {1};\n"
         << "Physical Line(32) = {2};\n"
         << "Physical Line(33) = {3};\n";

    nameStr << "triangle." << meshSize;
    t.restart();
    t_tot.restart();
    cout <<"Mesh generation ... ";
    fname = __gmsh.generate( nameStr.str(), ostr.str() );

    /* Mesh */

    typedef Mesh<GeoEntity<Simplex<2, 1> > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    mesh_ptrtype mesh( new mesh_type );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    timeout( t.elapsed() );

    /** Approximation Space **/

    using namespace Feel::vf;

    cout << "Space and forms generation ...\n";
    t.restart();
    typedef BoundaryAdaptedPolynomialSet<2, N, Scalar, value_type, Simplex> basis_type;
    typedef FunctionSpace<mesh_type, basis_type, value_type > space_type;
    boost::shared_ptr<space_type> Xh( new space_type( mesh ) );

    typedef MatrixGmm<value_type, gmm::row_major> gmm_matrix_type;


#if 0  /** Dof debugging **/

    if ( N == 4 )
    {
        Xh.get()->dof().get()->showMe( cout, true );
    }

#endif

    S_dofs << Xh.get()->nDof() << endl;

    IM_PK<2, N+N, value_type, Gauss> im;

    typename space_type::element_type uM( Xh.get() );
    typename space_type::element_type vM( Xh.get() );

    gmm_matrix_type Mass;
    BilinearForm<space_type, space_type, gmm_matrix_type> m( Xh, Xh, Mass );
    m =  integrate( elements( *mesh ), im, idt( uM )*id( vM ) );
    Mass.close();

    typename space_type::element_type uG( Xh.get() );
    typename space_type::element_type vG( Xh.get() );

    gmm_matrix_type Global;
    BilinearForm<space_type, space_type, gmm_matrix_type> g( Xh, Xh, Global );
    g =  integrate( elements( *mesh ), im, idt( uG )*id( vG ) + dxt( uG )*dx( vG ) + dyt( uG )*dy( vG ) );
    Global.close();

    // -- PROBLEM PARAMETERS --

#if 0 //sin(pi x)cos(pi y)
    value_type pi = 4.0 * math::atan( value_type( 1.0 ) );

    __typeof__( sin( pi*Px() )*cos( pi*Py() ) )
    exact_sol = sin( pi*Px() )*cos( pi*Py() );

    __typeof__( ( 2.0*pi*pi+1.0 )*sin( pi*Px() )*cos( pi*Py() ) )
    f = ( 2.0*pi*pi+1.0 )*sin( pi*Px() )*cos( pi*Py() );

    __typeof__( pi*cos( pi*Px() )*cos( pi*Py() ) )
    grad_x = pi*cos( pi*Px() )*cos( pi*Py() );

    __typeof__( -pi*sin( pi*Px() )*sin( pi*Py() ) )
    grad_y = -pi*sin( pi*Px() )*sin( pi*Py() );

    double tol =  10*pow( 1.0/pow( N,N ), ( 0.4/meshSize ) );
    cout << "\n10*N^{-2N/5h} = " << tol << endl;

#else // exp(-(x+y))
    __typeof__( exp( -Px()-Py() ) )
    exact_sol = exp( -Px()-Py() );

    __typeof__( -exp( -Px()-Py() ) )
    f = -exp( -Px()-Py() );

    __typeof__( -exp( -Px()-Py() ) )
    grad_x = -exp( -Px()-Py() );

    __typeof__( -exp( -Px()-Py() ) )
    grad_y = -exp( -Px()-Py() );

    double tol =  pow( 1.0/pow( N,N ), ( 0.25/meshSize ) );
    cout << "\nN^{-N/4h} = " << tol << endl;
#endif

    VectorUblas<value_type> F( Xh->nDof() );
    LinearForm<space_type, VectorUblas<value_type> > l( Xh, F );

    l = integrate( elements( *mesh ), im, f*id( vG ) )
        + integrate( markedfaces( *mesh,31 ), im, ( Nx()*grad_x + Ny()*grad_y )*id( vG ) )
        + integrate( markedfaces( *mesh,32 ), im, ( Nx()*grad_x + Ny()*grad_y )*id( vG ) )
        + integrate( markedfaces( *mesh,33 ), im, ( Nx()*grad_x + Ny()*grad_y )*id( vG ) );

    timeout( t.elapsed() );

    typename space_type::element_type vP( Xh.get() );
    VectorUblas<value_type> ProjF( Xh->nDof() );
    LinearForm<space_type, VectorUblas<value_type> > pro( Xh, ProjF );


    pro = integrate( elements( *mesh ), im, exact_sol*id( vP ) );

    gmm::iteration projiter( 2.0E-14 );
    projiter.set_noisy( 0 );
    projiter.set_maxiter( 20000 );

    typename space_type::element_type sol_proj( Xh.get() );

#if 0
    // incomplete LU with k fill-in and threshold preconditioner.
    gmm::ildltt_precond<typename gmm_matrix_type::matrix_type> PR2( Mass.mat(), 2, 1e-3 );
#else
    gmm::identity_matrix PR2;
#endif
    size_t restart= 100;
    gmm::gmres( Mass.mat(), sol_proj.container(), ProjF, PR2, restart,  projiter );

    if ( !projiter.converged() )
    {
        cout << "Projection Solver didn't converge" << endl;
    }

    else
    {
        cout << "Projection Solver converged in " << projiter.get_iteration() << " iterations and ";
        timeout( t.elapsed() );
    }

    /** System Resolution **/

    cout << "System iterative resolution ...";
    t.restart();
    typename space_type::element_type sol_ap( Xh.get() );

    gmm::iteration iter( 2.0E-14 );
    iter.set_noisy( 0 );
    iter.set_maxiter( 10000 );

#if 0
    // incomplete LU with k fill-in and threshold preconditioner.
    gmm::ildltt_precond<typename gmm_matrix_type::matrix_type> PR( Global.mat(), 2, 1e-3 );
#else
    gmm::identity_matrix PR;
#endif

    gmm::gmres( Global.mat(), sol_ap.container(), F, PR, restart,  iter );

    if ( !iter.converged() )
    {
        cout << "Problem Solver didn't converge" << endl;
        S_iter << "nan" << endl;
    }

    else
    {
        cout << "Problem Solver converged in " << iter.get_iteration() << " iterations and ";
        S_iter << iter.get_iteration() << endl;
        timeout( t.elapsed() );
    }


    /** Project the exact_solution on the space **/

    cout << "Error Analysis ...";

    typename space_type::element_type error( Xh.get() );
    error = sol_proj - sol_ap;

    value_type ErrH1 = sqrt( g( error,error ) );
    value_type ErrL2 = sqrt( m( error,error ) );

    cout << "\nError in H1-norm = " << ErrH1 << endl;
    cout << "Error in L2-norm = " << ErrL2 << endl;

    cout << "\nTolerance = " << tol << endl;

    S_tol << tol << endl;

    S_error << ErrH1 << endl;

    double time_tot( t_tot.elapsed() );

    cout << "Total resolution time = ";
    timeout( time_tot );

    S_timing << time_tot << endl;

    BOOST_CHECK( ( ErrH1 <= tol ) );
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

