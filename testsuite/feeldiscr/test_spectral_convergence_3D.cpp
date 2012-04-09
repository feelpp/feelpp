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
   \file test_spectral_convergence_3D.cpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2006-04-27
 */

/**
 * This file test the 3D basis.
 * It is yet not complete, but is a good point
 * to start with 3D-spectral in Feel.
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

/** To save the results of the execution you must create the directory Data/3D/ **/

/** Streams **/

ofstream S_error( "Data/3D/error_3D.dat" );
ofstream S_timing( "Data/3D/timing_3D.dat" );
ofstream S_dofs( "Data/3D/dofs_3D.dat" );
ofstream S_iter( "Data/3D/solver_iter_3D.dat" );
ofstream S_tol( "Data/3D/tol_3D.dat" );

typedef int16_type int_type;

void timeout( double time );

template<int_type N, typename T> void Poisson();

using boost::unit_test::test_suite;

template<int_type P>
void add_tests( test_suite* test )
{
    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Poisson<P,value_type> ) ) );

    add_tests<P+1>( test );
}

template<>
void add_tests<2>( test_suite* test )
{
    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Poisson<2,value_type> ) ) );
}

double hsize( 0.5 );

/* Main test function */

test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
    test_suite* test = BOOST_TEST_SUITE( "Spectral Accuracy (BA) test using Feel Language" );

#if defined(FEELPP_HAS_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start( &old_cw );
#endif

    if ( argc > 1 )
    {
        hsize = static_cast<double> ( atof( argv[1] ) );
    }

    add_tests<1>( test );
    return test;
}



template<int_type N, typename T>
void Poisson()
{
    using namespace std;
    cout << "N = " << N << std::endl;
    typedef T value_type;

    timer t_tot;
    timer t;

    Gmsh __gmsh;
    string fname;
    ostringstream ostr;
    ostringstream nameStr;

    double meshSize=hsize;

    std::cout << "3D case ...\n\n";

#if 1 // tetraedra
    ostr << "h=" << meshSize << ";\n"
         << "Point(1) = {-1.0,-1.0,-1.0,h};\n"
         << "Point(2) = {-1.0, 1.0,-1.0,h};\n"
         << "Point(3) = { 1.0,-1.0,-1.0,h};\n"
         << "Point(4) = {-1.0,-1.0, 1.0,h};\n"
         << "Line (1) = {1, 4};\n"
         << "Line (2) = {4, 3};\n"
         << "Line (3) = {3, 1};\n"
         << "Line (4) = {1, 2};\n"
         << "Line (5) = {4, 2};\n"
         << "Line (6) = {3, 2};\n"
         << "Line Loop (108) = {5, -6, -2};\n"
         << "Plane Surface (8) = {108};\n"
         << "Line Loop (109) = {6, -4, -3};\n"
         << "Plane Surface (10) = {109};\n"
         << "Line Loop (110) = {4, -5, -1};\n"
         << "Plane Surface (12) = {110};\n"
         << "Line Loop (111) = {3, 1, 2};\n"
         << "Plane Surface (14) = {111};\n"
         << "Surface Loop (15) = {8, 12, 10, 14};\n"
         << "Volume (26) = {15};\n"
#if 0
         << "Physical Line(12) = {1};\n"
         << "Physical Line(13) = {2};\n"
         << "Physical Line(14) = {3};\n"
         << "Physical Line(15) = {4};\n"
         << "Physical Line(16) = {5};\n"
         << "Physical Line(17) = {6};\n"
#endif
         << "Physical Surface(21) = {8};\n"
         << "Physical Surface(22) = {10};\n"
         << "Physical Surface(23) = {12};\n"
         << "Physical Surface(24) = {14};\n"
         << "Physical Volume (30) = {26};\n";

    nameStr << "tetra." << meshSize;
    t.restart();
    t_tot.restart();
    cout <<"Mesh generation ... ";
    fname = __gmsh.generate( nameStr.str(), ostr.str() );

#else
    ostr << "h=" << meshSize << ";\n"
         << "Point(1) = {-1,-1,-1,h};\n"
         << "Point(2) = {-1, 1,-1,h};\n"
         << "Point(3) = { 1, 1,-1,h};\n"
         << "Point(4) = { 1,-1,-1,h};\n"
         << "Line(1) = {1,4};\n"
         << "Line(2) = {4,3};\n"
         << "Line(3) = {3,2};\n"
         << "Line(4) = {2,1};\n"
         << "Line Loop(5) = {3,4,1,2};\n"
         << "Plane Surface(6) = {5};\n"
         << "Extrude Surface {6, {0,0,2}};\n"
         << "Physical Surface(10) = {15,23,6,28};\n"
         << "Physical Surface(20) = {19,27};\n"
         << "Surface Loop(31) = {28,15,-6,19,23,27};\n"
         << "Volume(1) = {31};\n"
         << "Physical Volume(2) = {1};\n";

    t.restart();
    t_tot.restart();
    nameStr << "cube." << meshSize;
    cout <<"Mesh generation ... ";
    fname = __gmsh.generate( nameStr.str(), ostr.str() );

#endif

    /* Mesh */

    typedef Mesh<GeoEntity<Simplex<3, 1> > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    mesh_ptrtype mesh( new mesh_type );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    timeout( t.elapsed() );

    /** Approximation Space **/

    using namespace Feel::vf;

    std::cout << "Space definition...\n";
    t.restart();

    /* Choose your basis ! */

#if 0 // Continuous Basis type
#if 0
    typedef fem::Lagrange<3, N, Scalar, Continuous, value_type, Simplex> basis_type;
#else
    typedef BoundaryAdaptedPolynomialSet<3, N, Scalar, value_type, Simplex> basis_type;
#endif
#else // Discontinuous Basis type
#if 0
    typedef fem::Lagrange<3, N, Scalar, Discontinuous, value_type, Simplex> basis_type;
#else
    typedef OrthonormalPolynomialSet<3, N, Scalar, value_type, Simplex> basis_type;
#endif
#endif

    typedef FunctionSpace<mesh_type, basis_type, value_type > space_type;
    boost::shared_ptr<space_type> Xh( new space_type( mesh ) );

    /*matrix*/

    typedef MatrixGmm<value_type, gmm::row_major> sparse_matrix_type;
    typedef VectorUblas<value_type> vector_type;

    S_dofs << Xh.get()->nDof() << std::endl;

    timeout( t.elapsed() );

    t.restart();
    std::cout << "Integration Method construction ... ";

    IM_PK<3, 2*N, value_type, Gauss> im;

    timeout( t.elapsed() );

    typename space_type::element_type u( Xh.get() );
    typename space_type::element_type v( Xh.get() );


    std::cout << "Mass matrix ... ";
    t.restart();

    sparse_matrix_type Mass;
    BilinearForm<space_type, space_type, sparse_matrix_type> m( Xh, Xh, Mass );
    m =  integrate( elements( *mesh ), im, idt( u )*id( v ) );
    Mass.close();

    //    Mass.printMatlab( "Mass.m" );
    timeout( t.elapsed() );

    std::cout << "Stiffness matrix ... ";
    t.restart();

    sparse_matrix_type Stiff;
    BilinearForm<space_type, space_type, sparse_matrix_type> sti( Xh, Xh, Stiff );
    sti =  integrate( elements( *mesh ), im, dxt( u )*dx( v ) + dyt( u )*dy( v ) + dzt( u )*dz( v ) );

    Stiff.close();
    timeout( t.elapsed() );

    std::cout << "Global matrix ... ";
    t.restart();

    sparse_matrix_type Global;
    BilinearForm<space_type, space_type, sparse_matrix_type> g( Xh, Xh, Global );
    g =  integrate( elements( *mesh ), im, idt( u )*id( v ) + dxt( u )*dx( v ) + dyt( u )*dy( v ) + dzt( u )*dz( v ) );
    Global.close();
    timeout( t.elapsed() );

    // -- PROBLEM PARAMETERS --

#if 1 //sin(pi x)cos(pi y)cos(pi z)

    value_type pi = 4.0 * math::atan( value_type( 1.0 ) );

    __typeof__( sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() ) )
    exact_sol = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );

    __typeof__( ( 3.0*pi*pi+1.0 )*sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() ) )
    f = ( 3.0*pi*pi+1.0 )*sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );

    __typeof__( pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() ) )
    grad_x = pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );

    __typeof__( -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() ) )
    grad_y =  -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() );

    __typeof__( -pi*sin( pi*Px() )*cos( pi*Py() )*sin( pi*Pz() ) )
    grad_z =  -pi*sin( pi*Px() )*cos( pi*Py() )*sin( pi*Pz() );

    double tol = 10.0;
#else
    __typeof__( exp( -Px()-Py()-Pz() )  )
    exact_sol = exp( -Px()-Py()-Pz() );

    __typeof__( -2.0*exp( -Px()-Py()-Pz() ) )
    f = -2.0*exp( -Px()-Py()-Pz() );

    __typeof__( -exp( -Px()-Py()-Pz() ) )
    grad_x = -exp( -Px()-Py()-Pz() );

    __typeof__( -exp( -Px()-Py()-Pz() ) )
    grad_y =  -exp( -Px()-Py()-Pz() );

    __typeof__( -exp( -Px()-Py()-Pz() ) )
    grad_z =  -exp( -Px()-Py()-Pz() );

    double tol = 10.0;
#endif

    value_type s3 = 1.0/sqrt( 3.0 );

    std::cout << "Right hand side ... ";
    t.restart();

    vector_type F( Xh->nDof() );
    LinearForm<space_type, vector_type> l( Xh, F );

    l = integrate( elements( *mesh ), im, f*id( v ) )
#if 1
        + integrate( markedfaces( *mesh,21 ), im, ( Nx()*grad_x + Ny()*grad_y + Nz()*grad_z )*id( v ) )
        + integrate( markedfaces( *mesh,22 ), im, ( Nx()*grad_x + Ny()*grad_y + Nz()*grad_z )*id( v ) )
        + integrate( markedfaces( *mesh,23 ), im, ( Nx()*grad_x + Ny()*grad_y + Nz()*grad_z )*id( v ) )
        + integrate( markedfaces( *mesh,24 ), im, ( Nx()*grad_x + Ny()*grad_y + Nz()*grad_z )*id( v ) );
#else
        + integrate( markedfaces( *mesh,21 ), im, ( s3*( grad_x + grad_y + grad_z ) )*id( v ) )
        + integrate( markedfaces( *mesh,22 ), im, ( -grad_z )*id( v ) )
        + integrate( markedfaces( *mesh,23 ), im, ( -grad_y )*id( v ) )
        + integrate( markedfaces( *mesh,24 ), im, ( -grad_x )*id( v ) );
#endif
    timeout( t.elapsed() );

    vector_type ProjF( Xh->nDof() );
    LinearForm<space_type, vector_type > pro( Xh, ProjF );

    pro = integrate( elements( *mesh ), im, exact_sol*id( v ) );
    ProjF.close();

    std::cout << "Projection Solver ... ";

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
        std::cout << time << " seconds.\n";

    else if ( time < 3600 )
    {
        int_type min = ( int_type )floor( time )/60;
        std::cout << min << " minutes " << ( time-60*min ) << " seconds.\n";
    }

    else
    {
        int_type hour = ( int_type )floor( time )/3600;
        int_type min = ( int_type )floor( time-hour*3600 )/60;
        std::cout << hour << " hours " << min << " minutes " << ( time-3600*hour-60*min ) << " seconds.\n";
    }
}

