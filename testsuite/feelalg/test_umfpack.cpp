/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-10-28

  Copyright (C) 2004 EPFL

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
   \file test_umfpack.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-10-28
 */
#include <stdio.h>

#include <boost/mpl/and.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <feel/feelconfig.h>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>

#if defined(FEELPP_HAS_BOOST_TEST)
// Boost.Test
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>

using boost::unit_test_framework::test_suite;

#if defined(FEELPP_HAS_UMFPACK_H)

#include <feel/feelalg/solverumfpack.hpp>


void test_csr_boost()
{
    int n = 3;
    using namespace boost::numeric::ublas;

    compressed_matrix<double, column_major> m ( n, n, n*n );

    for ( unsigned i = 0; i < m.size1 (); ++ i )
        for ( unsigned j = 0; j < m.size2 (); ++ j )
            m ( i, j ) = 3 * i + j;

    std::cout << m << std::endl;

    int status;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];
    void *Symbolic, *Numeric ;

    status = umfpack_dl_symbolic ( n, n, ( const long int* )&m.index1_data()[0], ( const long int* )&m.index2_data()[0], &m.value_data()[0],
                                   &Symbolic, Control, Info ) ;
    status = umfpack_dl_numeric ( ( const long int* )&m.index1_data()[0], ( const long int* )&m.index2_data()[0], &m.value_data()[0],
                                  Symbolic, &Numeric, Control, Info ) ;
    umfpack_di_defaults ( Control ) ;

    vector<double> x( n );
    x = scalar_vector<double>( n, 1 );
    vector<double> b( n );

    b = prod( m, x );
    DVLOG(2) << "b = A * x\n";

    status = umfpack_dl_solve ( UMFPACK_A, ( const long int* )&m.index1_data()[0], ( const long int* )&m.index2_data()[0], &m.value_data()[0],
                                &x[0], &b[0], Numeric, Control, Info ) ;


    DVLOG(2) << "solver A x = b\n";

    std::cout << "B = " << b << "\n";
    std::cout << "AX = " << prod( m, x ) << "\n";
    std::cout << "x = " << x << "\n";
    DVLOG(2) << "norm_2( Ax - b ) = " << norm_2( prod( m, x ) - b ) << "\n";

}
void test_convdiff()
{
    int n = 4;
    long Ap[] = {0, 3, 6, 9, 12};
    long Ai[] = { 0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3};

    /*
      36 -10  -9  0
      -7  36  0  -9
      -9  0  36 -10
      0   -9  -7  36
      A=[[36 -10  -9  0]
      [-10  36  0  -9]
      [-9   0  36 -10]
      [0   -9  -10  36]]
     */
    /*
      double Ax[] = { 36, -10, -9,
                    -10,  36,     -9,
                    -9,      36, -10,
                        -9,  -10, 36};
    */
    double Ax[] = { 36, 0, 0,
                    1,  36,     0,
                    0,      36, 0,
                    0,  0, 36
                  };
    using namespace Feel;

    std::vector<uint> __ap( n+1 ), __ai( 12 );
    std::copy( Ap, Ap+n+1, __ap.begin() );
    std::copy( Ai, Ap+12, __ai.begin() );
    std::vector<double> __ax( 12 );
    std::copy( Ax, Ax+12,  __ax.begin() );

    DVLOG(2) << "copying done\n";
    ublas::compressed_matrix<double, ublas::column_major> __matrix( n, n );

    for ( int i = 0; i < n; ++i )
    {
        int j = Ap[i];

        while ( j != Ap[i+1] )
        {
            __matrix( Ai[j], i ) = Ax[j];
            ++j;
        }
    }

#if 0
    SolverUMFPACK __solver;
    __solver.setMatrix( n, Ap, Ai, Ax );
#else
    int status;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];
    void *Symbolic, *Numeric ;

    status = umfpack_dl_symbolic ( n, n, Ap, Ai, Ax, &Symbolic, Control, Info ) ;
    status = umfpack_dl_numeric ( Ap, Ai, Ax, Symbolic, &Numeric, Control, Info ) ;
    umfpack_di_defaults ( Control ) ;
#endif

    DVLOG(2) << "UMFPACK solver created\n";

    ublas::vector<double> x( n );
    x = ublas::scalar_vector<double>( n, 1 );
    ublas::vector<double> b( n );

    b = ublas::prod( __matrix, x );
    //DVLOG(2) << "b = A * x\n";

    // solution should be x=[1 1 1 1 1]^T

    //__solver.solve( x, b );
    status = umfpack_dl_solve ( UMFPACK_A, Ap, Ai, Ax, &x[0], &b[0], Numeric, Control, Info ) ;


    DVLOG(2) << "solver A x = b\n";

    ublas::vector<double> AX( n );
    AX = ublas::prod( __matrix, x );
    std::cout << "B = " << b << "\n";
    std::cout << "AX = " << AX << "\n";
    std::cout << "x = " << x << "\n";
    DVLOG(2) << "norm_2( Ax - b ) = " << norm_2( AX - b ) << "\n";

    ublas::vector<double> sol( n );
    sol[0]=1;
    sol[1]=1;
    sol[2]=1;
    sol[3]=1;
    DVLOG(2) << "norm_2( sol - x ) = " << norm_2( sol - x ) << "\n";
    BOOST_REQUIRE( norm_2( sol - x ) < 1e-10 );
}
void test_umfpack()
{
    using namespace Feel;

    int    n=5;

    /**
       2  3  0 0 0
       3  0  4 0 6
    A= 0 -1 -3 2 0 .
       0  0  1 0 0
       0  4  2 0 1
    */
    int    Ap [ ] = {0, 2, 5, 9, 10, 12} ;
    int   Ai [ ] = { 0, 1, 0,     2, 4, 1, 2, 3,       4, 2, 1, 4} ;
    double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;

    SolverUMFPACK __solver;

    ublas::compressed_matrix<double, ublas::column_major> __matrix( n, n );

    for ( int i = 0; i < n; ++i )
    {
        int j = Ap[i];

        while ( j != Ap[i+1] )
        {
            __matrix( Ai[j], i ) = Ax[j];
            ++j;
        }
    }

    __solver.setMatrix( __matrix );
    ublas::vector<double> x( n );
    ublas::vector<double> b( n );
    b[0] = 8.;
    b[1] = 45.;
    b[2] = -3.;
    b[3] = 3.;
    b[4] = 19. ;

    // solution should be x=[1 2 3 4 5]^T

    __solver.solve( x, b );

    ublas::vector<double> sol( n );
    sol[0]=1;
    sol[1]=2;
    sol[2]=3;
    sol[3]=4;
    sol[4]=5;
    DVLOG(2) << "norm_2( sol - x ) = " << norm_2( sol - x ) << "\n";
    BOOST_REQUIRE( norm_2( sol - x ) < 1e-10 );
}

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{
    test_suite* test= BOOST_TEST_SUITE( "UMFPACK Unit Test" );

    // this example will pass cause we know ahead of time number of expected failures
    test->add( BOOST_TEST_CASE( &test_umfpack ), 0 );
    test->add( BOOST_TEST_CASE( &test_csr_boost ), 0 );
    //test->add( BOOST_TEST_CASE( &test_convdiff ), 0 );

    return test;
}
#else
test_suite*
init_unit_test_suite( int argc, char** argv )
{
    test_suite* test= BOOST_TEST_SUITE( "UMFPACK Unit Test" );
    return test;
}
#endif  /* FEELPP_HAS_UMFPACK_H */
#else
int main()
{
    return EXIT_SUCCESS;
}
#endif

