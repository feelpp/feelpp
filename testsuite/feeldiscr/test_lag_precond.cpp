/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
       Date: 2006-04-28

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
   \file test_lag_precond.cpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2006-04-28
 */
#include <iostream>
#include <sstream>
#include <fstream>

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
#include <feel/feelfilters/exporterensight.hpp>

#include <gmm_dense_lu.h>
#include <gmm_iter_solvers.h>

using boost::unit_test::test_suite;
using boost::timer;

/**
 * The definition of a new preconditioneur
 * for the iterative solver in gmm needs
 * the definition of the function mult.
 * Here is an example.
 * Warning the p1-lagrange associated problem
 * seems to furnish singular caracteristic matrices.
 **/

namespace gmm
{

template<typename Matrix> struct lag_precond
{

    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef typename linalg_traits<Matrix>::this_type matrix_type;

    matrix_type mat;

    void build_with( const Matrix &M )
    {
        mat.init_with( M );
    }

    lag_precond( const Matrix &M )
    {
        build_with( M );
    }
    lag_precond( void ) {}
};


template <typename Matrix, typename V1, typename V2> inline
void mult( const lag_precond<Matrix>& P, const V1 &v1, V2 &v2 )
{
    //   copy(v1, v2);

    gmm::iteration liter( 2.0E-15 );
    liter.set_noisy( 1 );
    liter.set_maxiter( 100 );

    gmm::identity_matrix PREC;

    size_t restart = 100;

    /** If gmres do not converge, check the matrix P.mat **/

    gmm::gmres( P.mat, v2, v1, PREC, restart, liter );
}
}

using namespace std;
using namespace Feel;

/** Streams **/

ofstream S_error( "Data/Precond/error_2D.dat" );
ofstream S_timing( "Data/Precond/timing_2D.dat" );
ofstream S_dofs( "Data/Precond/dofs_2D.dat" );
ofstream S_iter( "Data/Precond/solver_iter_2D.dat" );
ofstream S_h( "Data/Precond/h_2D.dat" );

typedef int16_type int_type;

template<int_type N, typename T> T Poisson();

void timeout( double time );

template<int_type P>
void add_tests( test_suite* test )
{
    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Poisson<P,value_type> ) ) );

    add_tests<P+2>( test );
}

template<>
void add_tests<4>( test_suite* test )
{
    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Poisson<4,value_type> ) ) );
}


/* Main test function */

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv[]*/ )
{
    test_suite* test = BOOST_TEST_SUITE( "P1-Lagrange preconditionner" );

    add_tests<2>( test );
    return test;
}


template<int_type N, typename T>
T Poisson()
{
    cout << "N = " << N << endl;

    typedef T value_type;

    timer t_tot;
    timer t;

    const int_type Order = N;

#if 0
    double h=( 4.0 )/value_type( N );
#else
    double h=1.0;
#endif

    S_h << h << std::endl;
    Gmsh __gmsh;
    string fname;

    ostringstream ostr;

    ostr << "h=" << h << ";\n"
         << "Point(1) = {-1.0,-1.0,0.0,h};\n"
         << "Point(2) = { 1.0,-1.0,0.0,h};\n"
         << "Point(3) = {-1.0, 1.0,0.0,h};\n"
         << "Line(1) = {2,3};\n"
         << "Line(2) = {3,1};\n"
         << "Line(3) = {1,2};\n"
         << "Line Loop(4) = {1,2,3};\n"
         << "Plane Surface(5) = {4};\n"
         << "Physical Surface(30) = {5};\n"
         << "Physical Line(31) = {1};\n"
         << "Physical Line(32) = {2};\n"
         << "Physical Line(33) = {3};\n";
    t.restart();
    cout <<"Mesh generation ... ";
    fname = __gmsh.generate( "triangle", ostr.str() );

    /*mesh*/

    typedef Mesh<GeoEntity<Simplex<2, 1> > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    mesh_ptrtype mesh( new mesh_type );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    mesh->partition();
    timeout( t.elapsed() );

    /** Approximation Space **/

    using namespace Feel::vf;

    cout << "Space and forms generation ...\n";
    t.restart();
    typedef BoundaryAdaptedPolynomialSet<2, Order, Scalar, value_type, Simplex> basis_type;

    typedef FunctionSpace<mesh_type, basis_type, value_type > space_type;
    boost::shared_ptr<space_type> Xh( new space_type( mesh ) );

    /* p1-Lagrange associated space */

    typedef typename space_type::P1Lagrange p1_lag_type;
    typedef typename space_type::P1Lagrange::p1_type p1_type;

    boost::shared_ptr<p1_lag_type> over_p1lag( new p1_lag_type( Xh ) );
    boost::shared_ptr<p1_type> p1lag = over_p1lag->space();

    typedef MatrixGmm<value_type, gmm::row_major> gmm_matrix_type;

    /** Dof debugging **/

#if 0
    Xh.get()->dof().get()->showMe( cout, true );
#endif

    S_dofs << Xh.get()->nDof() << endl;

    cout << "Integration method ...\n";
    IM_PK<2, 2*Order, value_type, Gauss> im;

    typename space_type::element_type u( Xh.get() );
    typename space_type::element_type v( Xh.get() );

    typename p1_type::element_type u_L( p1lag.get() );
    typename p1_type::element_type v_L( p1lag.get() );

    cout << "Mass matrix ...\n";
    gmm_matrix_type Mass;
    BilinearForm<space_type, space_type, gmm_matrix_type> m( Xh, Xh, Mass );
    m =  integrate( elements( *mesh ), im, idt( u )*id( v ) );
    Mass.close();

    cout << "Lagrange Mass Preconditionner ...\n";
    gmm_matrix_type L_Mass;
    BilinearForm<p1_type, p1_type, gmm_matrix_type> m_lag( p1lag, p1lag, L_Mass );
    m_lag =  integrate( elements( *mesh ), im, idt( u_L )*id( v_L ) );
    L_Mass.close();

    cout << "Stiffness matrix ...\n";
    gmm_matrix_type Stiff;
    BilinearForm<space_type, space_type, gmm_matrix_type> st( Xh, Xh, Stiff );
    st =  integrate( elements( *mesh ), im, dxt( u )*dx( v ) + dyt( u )*dy( v ) );
    Stiff.close();

#if 0
    Mass.printMatlab( "Mass" );
#endif

    cout << "Global matrix ...\n";
    gmm_matrix_type Global;
    BilinearForm<space_type, space_type, gmm_matrix_type> g( Xh, Xh, Global );
    g =  integrate( elements( *mesh ), im, idt( u )*id( v ) + dxt( u )*dx( v ) + dyt( u )*dy( v ) );
    Global.close();

    cout << "Lagrange Global Preconditionner ...\n";
    gmm_matrix_type L_Global;
    BilinearForm<p1_type, p1_type, gmm_matrix_type> g_lag( p1lag, p1lag, L_Global );
    g_lag =  integrate( elements( *mesh ), im, idt( u_L )*id( v_L ) + dxt( u_L )*dx( v_L ) + dyt( u_L )*dy( v_L ) );
    L_Global.close();

    // -- PROBLEM PARAMETERS --

    value_type pi = 4.0 * math::atan( value_type( 1.0 ) );

    __typeof__( sin( pi*Px() )*cos( pi*Py() ) )
    exact_sol = sin( pi*Px() )*cos( pi*Py() );

    __typeof__( ( 2.0*pi*pi+1.0 )*sin( pi*Px() )*cos( pi*Py() ) )
    f = ( 2.0*pi*pi+1.0 )*sin( pi*Px() )*cos( pi*Py() );

    __typeof__( pi*cos( pi*Px() )*cos( pi*Py() ) )
    grad_x = pi*cos( pi*Px() )*cos( pi*Py() );

    __typeof__( -pi*sin( pi*Px() )*sin( pi*Py() ) )
    grad_y = -pi*sin( pi*Px() )*sin( pi*Py() );

    VectorUblas<value_type> F( Xh->nDof() );
    LinearForm<space_type, VectorUblas<value_type> > l( Xh, F );

    l = integrate( elements( *mesh ), im, f*id( v ) )
        + integrate( markedfaces( *mesh,31 ), im, ( Nx()*grad_x + Ny()*grad_y )*id( v ) )
        + integrate( markedfaces( *mesh,32 ), im, ( Nx()*grad_x + Ny()*grad_y )*id( v ) )
        + integrate( markedfaces( *mesh,33 ), im, ( Nx()*grad_x + Ny()*grad_y )*id( v ) );

    timeout( t.elapsed() );

    VectorUblas<value_type> ProjF( Xh->nDof() );
    LinearForm<space_type, VectorUblas<value_type> > pro( Xh, ProjF );

    pro = integrate( elements( *mesh ), im, exact_sol*id( v ) );

    gmm::iteration projiter( 2.0E-15 );
    projiter.set_noisy( 0 );
    projiter.set_maxiter( 5000 );

    typename space_type::element_type sol_proj( Xh.get() );

#if 1 // The lagrange Mass Matrix seems to be singular
    cout << "Spectral Mass Matrix : " << Mass.mat() << std::endl;
    cout << "Lagrange Mass Matrix : " << L_Mass.mat() << std::endl;
    gmm::lag_precond<typename gmm_matrix_type::matrix_type> PR2( L_Mass.mat() );
#else
    gmm::identity_matrix PR2;
    //  gmm::ildltt_precond<typename gmm_matrix_type::matrix_type> PR2(L_Mass.mat(), 2, 1e-3);
#endif

    size_t restart = 200;
    gmm::gmres( Mass.mat(), sol_proj.container(), ProjF, PR2, restart, projiter );

    if ( !projiter.converged() )
    {
        cout << "Projection Solver didn't converge" << endl;
    }

    else
    {
        cout << "Projection Solver converged in " << projiter.get_iteration() << " iterations and ";
        timeout( t.elapsed() );
    }

    // cout << sol_proj.container() << endl;


    /** System Resolution **/

    cout << "System iterative resolution ...";
    t.restart();
    typename space_type::element_type sol_ap( Xh.get() );

    gmm::iteration iter( 2.0E-15 );
    iter.set_noisy( 0 );
    iter.set_maxiter( 5000 );

#if 1
    cout << "Lagrange Global Matrix : " << L_Global.mat() << std::endl;
    gmm::lag_precond<typename gmm_matrix_type::matrix_type> PR( L_Global.mat() );
#else
    gmm::identity_matrix PR;
    //  gmm::ildltt_precond<typename gmm_matrix_type::matrix_type> PR(L_Global.mat(), 2, 1e-3);
#endif

    size_t restart2= 100;
    gmm::gmres( Global.mat(), sol_ap.container(), F, PR, restart2,  iter );

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

    typename space_type::element_type error( Xh.get() );
    error = sol_proj - sol_ap;

    value_type ErrH1 = sqrt( g( error,error ) );
    value_type ErrL2 = sqrt( m( error,error ) );

    cout << "Error in H1-norm = " << ErrH1 << endl;
    cout << "Error in L2-norm = " << ErrL2 << endl;

    S_error << ErrH1 << endl;

    double time_tot( t_tot.elapsed() );

    cout << "Total resolution time = ";
    timeout( time_tot );

    S_timing << time_tot << endl;

    return ErrH1;

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

