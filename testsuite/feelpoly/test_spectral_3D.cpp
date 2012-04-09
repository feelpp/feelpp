/**
   \file spectral_3D.cpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \date 2006-03-02
 */

/**
    The aim of this file is to test the spectral resolution in 3D on the tetraedra.
    We are particularly interested in the comparison of the differant basis (orthonormal, orthogonal, boundary adapted),
    in the performance of the code, and in the observation of the spectral accuracy.
**/


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
- \Lambda u + u = f
$$
**/

/** Streams **/

ofstream B_error( "Data/3D/Boundary_error_3D.dat" );
ofstream D_error( "Data/3D/Dubiner_error_3D.dat" );
ofstream B_error_L2( "Data/3D/Boundary_error_L2_3D.dat" );
ofstream D_error_L2( "Data/3D/Dubiner_error_L2_3D.dat" );
ofstream B_cond( "Data/3D/Boundary_cond_3D.dat" );
ofstream D_cond( "Data/3D/Dubiner_cond_3D.dat" );
ofstream B_timing( "Data/3D/Boundary_timing_3D.dat" );
ofstream D_timing( "Data/3D/Dubiner_timing_3D.dat" );

template<int_type N, typename T> T Boundary();
template<int_type N, typename T> T Ortho();

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

#define ORDER 5
#else

#define ORDER 9

#endif

template<>
void add_tests<ORDER>( test_suite* test )
{
    using namespace Feel;
    using namespace std;

    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Boundary<ORDER,value_type> ) ) );
    test->add( BOOST_TEST_CASE( ( Ortho<ORDER,value_type> ) ) );
}

#define BOUNDARY 1

/* Main test function */

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv[]*/ )
{
    test_suite* test = BOOST_TEST_SUITE( "3D spectral tests" );
#if 0
    unsigned int old_cw;
    fpu_fix_start( &old_cw );
#endif
    add_tests<1>( test );
    return test;
}

// sin(pi x) cos(pi y) cos(pi z)

#define pi M_PI

template<typename T>
T f( typename node<T>::type const& t )
{
    return ( sin( pi*t[0] )*cos( pi*t[1] )*cos( pi*t[2] )*( 3*pi*pi+1 ) ) ;
}

template<typename T>
T g0( typename node<T>::type const& t )
{
    return ( pi*( cos( pi*t[0] )*cos( pi*t[1] )*cos( pi*t[2] ) - sin( pi*t[0] )*sin( pi*t[1] )*cos( pi*t[2] ) - sin( pi*t[0] )*cos( pi*t[1] )*sin( pi*t[2] ) )/( sqrt( 3.0 ) ) );
}

template<typename T>
T g1( typename node<T>::type const& t )
{
    return ( -pi*( cos( pi*t[0] )*cos( pi*t[1] )*cos( pi*t[2] ) ) );
}

template<typename T>
T g2( typename node<T>::type const& t )
{
    return ( pi*( sin( pi*t[0] )*sin( pi*t[1] )*cos( pi*t[2] ) ) );
}

template<typename T>
T g3( typename node<T>::type const& t )
{
    return ( pi*( sin( pi*t[0] )*cos( pi*t[1] )*sin( pi*t[2] ) ) );
}

template<typename T>
T exact_sol( typename node<T>::type const& t )
{
    return ( sin( pi*t[0] )*cos( pi*t[1] )*cos( pi*t[2] ) );
}

template<typename T>
T dexact_sol_x( typename node<T>::type const& t )
{
    return ( pi*cos( pi*t[0] )*cos( pi*t[1] )*cos( pi*t[2] ) );
}

template<typename T>
T dexact_sol_y( typename node<T>::type const& t )
{
    return ( -pi*sin( pi*t[0] )*sin( pi*t[1] )*cos( pi*t[2] ) );
}

template<typename T>
T dexact_sol_z( typename node<T>::type const& t )
{
    return ( -pi*sin( pi*t[0] )*cos( pi*t[1] )*sin( pi*t[2] ) );
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
void PrintMatlab( ublas::matrix<T> const & M )
{
    ofstream str( "Data/3D/Mat.m" );

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

/** Dubiner resolution **/

template<int_type N, typename T>
T Ortho()
{
    cout << "\nDubiner Basis !\n";
    cout << "N = " << N << endl;

    typedef T value_type;

    cout << "Construction of the Integration Method ... ";

    timer t_tot;
    timer t;

    const Gauss<Simplex<3,1>, 2*N, value_type> Quad;
    cout << t.elapsed() << " seconds.\n";

    cout << "Basis Construction ... ";
    t.restart();

#if 0
    const BoundaryAdaptedPolynomialSet<3, N, Scalar, value_type, Simplex> Basis;
#else
    const OrthogonalPolynomialSet<3, N, Scalar, value_type, Simplex> Basis;
#endif

    timeout( t.elapsed() );

    const int_type Dim_Basis = Basis.basis().coeff().size1();

    /** Matrix construction **/

    /** Stiffness **/
    cout << "\n/********************************/" << endl;
    cout << "Construction of Stiffness Matrix..."  <<endl;

    cout <<"Basis Derivation... ";
    t.restart();
    ublas::vector<ublas::matrix<value_type> > Diff( Basis.derivate( Quad.points() ) );
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

    ublas::matrix<value_type> Stiff( A + B + C );
    timeout( t.elapsed() );

    // PrintMatlab<value_type>(Stiff);

    /** Mass **/

    cout << "\n/********************************/" << endl;
    cout << "Construction of the Mass Matrix..." << endl;

    cout <<"Basis Evaluation... ";
    t.restart();
    ublas::matrix<value_type> Psi = Basis.evaluate( Quad.points() );
    timeout( t.elapsed() );

    cout << "Mass numeric construction ... ";
    t.restart();
    ublas::matrix<value_type> Mass ( ublas::prod( Psi, Wgts ) );
    Mass = ublas::prod( Mass, ublas::trans( Psi ) );

    //  PrintMatlab<value_type>(Mass);

    timeout( t.elapsed() );


    /** Test the mass matrix : We compute the volume of the tetrahedra **/

    ublas::unit_vector<value_type> one( Mass.size1() , 0 );

    ublas::vector<value_type> tmp( ublas::prod( Mass,one ) );

#if 1
    cout << "\nVolume tetra = " << ublas::inner_prod( one,tmp ) << endl;
#endif

    /** Global **/

    cout << "\n/********************************/" << endl;
    cout << "Assembly of the Global matrix ..."<<endl;
    ublas::matrix<value_type> Global( Stiff + Mass );

    glas::clean<ublas::matrix<value_type> >( Global, 1e-15 );


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
    ublas::vector<value_type> G3( Quad.fpoints( 3 ).size2() );

    ublas::diagonal_matrix<value_type> Wgts0( Quad.weights( 0 ).size() );
    ublas::diagonal_matrix<value_type> Wgts1( Quad.weights( 1 ).size() );
    ublas::diagonal_matrix<value_type> Wgts2( Quad.weights( 2 ).size() );
    ublas::diagonal_matrix<value_type> Wgts3( Quad.weights( 3 ).size() );

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
    G0= ublas::prod( Basis.evaluate( Quad.fpoints( 0 ) ),G0 );

    for ( int_type i=0; i < Wgts1.size1() ; ++i )
    {
        Wgts1( i,i ) = Quad.weights( 1 )( i );
        G1( i ) = g1<value_type>( ublas::column( Quad.fpoints( 1 ), i ) );
    }

    G1= ublas::prod( Wgts1,G1 );
    G1= ublas::prod( Basis.evaluate( Quad.fpoints( 1 ) ),G1 );

    for ( int_type i=0; i < Wgts2.size1() ; ++i )
    {
        Wgts2( i,i ) = Quad.weights( 2 )( i );
        G2( i ) = g2<value_type>( ublas::column( Quad.fpoints( 2 ), i ) );
    }

    G2= ublas::prod( Wgts2,G2 );
    G2= ublas::prod( Basis.evaluate( Quad.fpoints( 2 ) ),G2 );

    for ( int_type i=0; i < Wgts3.size1() ; ++i )
    {
        Wgts3( i,i ) = Quad.weights( 3 )( i );
        G3( i ) = g3<value_type>( ublas::column( Quad.fpoints( 3 ), i ) );
    }

    G3= ublas::prod( Wgts3,G3 );
    G3= ublas::prod( Basis.evaluate( Quad.fpoints( 3 ) ),G3 );

    /** Global right-hand side assembly **/

    F+=G0+G1+G2+G3;

    //  glas::clean<ublas::vector<value_type> >(F, 1e-15);
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
    D_error << err_H1 << endl;

    cout << "\n-------------------------------\n";
    cout << "Total runtime = ";
    double timing( t_tot.elapsed() );

    timeout( timing );
    cout << "-------------------------------\n";
    cout << "\n\n";
    D_timing << timing << endl;

    return err_H1;
}

/** Boundary Adapted resolution **/

template<int_type N, typename T>
T Boundary()
{
    cout << "\nBoundary Adapted Basis !\n";
    cout << "N = " << N << endl;

    typedef T value_type;

    cout << "Construction of the Integration Method ... ";

    timer t_tot;
    timer t;

    const Gauss<Simplex<3,1>, 2*N, value_type> Quad;
    cout << t.elapsed() << " seconds.\n";

    cout << "Basis Construction ... ";
    t.restart();

#if 1
    const BoundaryAdaptedPolynomialSet<3, N, Scalar, value_type, Simplex> Basis;
#else
    const OrthogonalPolynomialSet<3, N, Scalar, value_type, Simplex> Basis;
#endif

    timeout( t.elapsed() );

    const int_type Dim_Basis = Basis.basis().coeff().size1();

    /** Matrix construction **/

    /** Stiffness **/
    cout << "\n/********************************/" << endl;
    cout << "Construction of Stiffness Matrix..."  <<endl;

    cout <<"Basis Derivation... ";
    t.restart();
    ublas::vector<ublas::matrix<value_type> > Diff( Basis.derivate( Quad.points() ) );
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

    ublas::matrix<value_type> Stiff( A + B + C );
    timeout( t.elapsed() );

    // PrintMatlab<value_type>(Stiff);

    /** Mass **/

    cout << "\n/********************************/" << endl;
    cout << "Construction of the Mass Matrix..." << endl;

    cout <<"Basis Evaluation... ";
    t.restart();
    ublas::matrix<value_type> Psi = Basis.evaluate( Quad.points() );
    //  glas::clean<ublas::matrix<value_type> >(Psi, 1e-15);
    timeout( t.elapsed() );

    cout << "Mass numeric construction ... ";
    t.restart();
    ublas::matrix<value_type> Mass ( ublas::prod( Psi, Wgts ) );
    Mass = ublas::prod( Mass, ublas::trans( Psi ) );

    glas::clean<ublas::matrix<value_type> >( Mass, 1e-15 );

    if ( N== 4 )
        PrintMatlab<value_type>( Mass );

    timeout( t.elapsed() );

#if 0
    /** Test the mass matrix : We compute the volume of the tetrahedra **/

    ublas::unit_vector<value_type> one( Mass.size1() , 0 );

    ublas::vector<value_type> tmp( ublas::prod( Mass,one ) );

    cout << "\nVolume tetra = " << ublas::inner_prod( one,tmp ) << endl;
#endif

    /** Global **/

    cout << "\n/********************************/" << endl;
    cout << "Assembly of the Global matrix ..."<<endl;
    ublas::matrix<value_type> Global( Stiff + Mass );

    glas::clean<ublas::matrix<value_type> >( Global, 1e-15 );

    /** Condition number estimation **/

    cout << "Condition number estimation ... ";
    t.restart();
    value_type condition = cond2<value_type>( Global );
    timeout( t.elapsed() );

    cout << "Condition number = " << condition << endl;
    B_cond << condition << endl;

    /** Right hand side **/

    cout << "Construction of the RHS ... ";
    t.restart();
    ublas::vector<value_type> F( Quad.points().size2() );
    ublas::vector<value_type> G0( Quad.fpoints( 0 ).size2() );
    ublas::vector<value_type> G1( Quad.fpoints( 1 ).size2() );
    ublas::vector<value_type> G2( Quad.fpoints( 2 ).size2() );
    ublas::vector<value_type> G3( Quad.fpoints( 3 ).size2() );

    ublas::diagonal_matrix<value_type> Wgts0( Quad.weights( 0 ).size() );
    ublas::diagonal_matrix<value_type> Wgts1( Quad.weights( 1 ).size() );
    ublas::diagonal_matrix<value_type> Wgts2( Quad.weights( 2 ).size() );
    ublas::diagonal_matrix<value_type> Wgts3( Quad.weights( 3 ).size() );

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
    G0= ublas::prod( Basis.evaluate( Quad.fpoints( 0 ) ),G0 );

    for ( int_type i=0; i < Wgts1.size1() ; ++i )
    {
        Wgts1( i,i ) = Quad.weights( 1 )( i );
        G1( i ) = g1<value_type>( ublas::column( Quad.fpoints( 1 ), i ) );
    }

    G1= ublas::prod( Wgts1,G1 );
    G1= ublas::prod( Basis.evaluate( Quad.fpoints( 1 ) ),G1 );

    for ( int_type i=0; i < Wgts2.size1() ; ++i )
    {
        Wgts2( i,i ) = Quad.weights( 2 )( i );
        G2( i ) = g2<value_type>( ublas::column( Quad.fpoints( 2 ), i ) );
    }

    G2= ublas::prod( Wgts2,G2 );
    G2= ublas::prod( Basis.evaluate( Quad.fpoints( 2 ) ),G2 );

    for ( int_type i=0; i < Wgts3.size1() ; ++i )
    {
        Wgts3( i,i ) = Quad.weights( 3 )( i );
        G3( i ) = g3<value_type>( ublas::column( Quad.fpoints( 3 ), i ) );
    }

    G3= ublas::prod( Wgts3,G3 );
    G3= ublas::prod( Basis.evaluate( Quad.fpoints( 3 ) ),G3 );

    /** Global right-hand side assembly **/

    F+=G0+G1+G2+G3;

    //  glas::clean<ublas::vector<value_type> >(F, 1e-15);
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

    cout << "\n-------------------------------\n";
    cout << "Total runtime = ";
    double timing( t_tot.elapsed() );

    timeout( timing );
    cout << "-------------------------------\n";
    cout << "\n\n";
    B_timing << timing << endl;

    return err_H1;
}
