/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-17

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
   \file test_poly.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-17
 */
#define USE_TEST
// Boost.Test
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <feel/feelpoly/polynomial.hpp>
#include <feel/feelpoly/operations.hpp>

using namespace Feel;
struct x2
{
    typedef double value_type;
    ublas::vector<double>
    operator()( ublas::matrix<double,ublas::column_major> const& m ) const
    {
        ublas::vector<double> v( m.size2() );
        v.clear();

        for ( int i = 0; i < m.size2(); ++i )
        {
            for ( int j = 0; j < m.size1(); ++j )
                v( i ) += ( j+2 )*m( j,i )*m( j,i );

            if ( m.size1() == 2 )
                v( i ) += 100*m( 0,i )*m( 1,i )+200*m( 1,i )*m( 0,i );

            //v( i ) = m(0,i)*m(0,i);
        }

        return v;
    }
};
BOOST_AUTO_TEST_CASE( test_polynomial1D )
{


    ublas::matrix<double,ublas::column_major> pts( 1, 3 );
    pts( 0,0 ) = 0;
    //pts(1,0 ) = 0;
    pts( 0,1 ) = -1;
    pts( 0,2 ) = 1;
    detail::OrthonormalPolynomialSet<1,1,2,Scalar> ps;
    Polynomial<detail::OrthonormalPolynomialSet<1,1,2,Scalar>, Scalar> p1 = project( ps, x2(), IM<1,3>() );
    PolynomialSet<detail::OrthonormalPolynomialSet<1,1,2,Scalar>, Scalar> p( ps );
    p.insert( p1.toSet( true ) );
    p.insert( p1.toSet( true ) );
    std::cout << "p=" << p.coeff() << "\n"
              << "grad(p)(0)=" << p.gradient().evaluate( pts ) << "\n";
    std::cout << "hess(p)=" << p.gradient().gradient().coeff() << "\n"
              << "hess(p)(0)=" << p.gradient().gradient().evaluate( pts ) << "\n";

    std::cout << "hess(p1) = " << p1.toSet( true ).gradient().gradient().evaluate( pts ) << "\n";
}

BOOST_AUTO_TEST_CASE( test_polynomial2D )
{


    ublas::matrix<double,ublas::column_major> pts( 2, 4 );
    pts( 0,0 ) = -1;
    pts( 1,0 ) = -1;
    pts( 0,1 ) =  1;
    pts( 1,1 ) = -1;
    pts( 0,2 ) = -1;
    pts( 1,2 ) =  1;
    pts( 0,3 ) =  0;
    pts( 1,3 ) =  0;
    detail::OrthonormalPolynomialSet<2,2,2,Scalar> ps;
    Polynomial<detail::OrthonormalPolynomialSet<2,2,2,Scalar>, Scalar> p1 = project( ps, x2(), IM<2,3>() );
    PolynomialSet<detail::OrthonormalPolynomialSet<2,2,2,Scalar>, Scalar> p( ps );
    p.insert( p1.toSet( true ) );
    p.insert( p1.toSet( true ) );
    std::cout << "p=" << p.coeff() << "\n"
              << "grad(p)(0)=" << p.gradient().evaluate( pts ) << "\n";
    std::cout << "hess(p)=" << p.gradient().gradient().coeff() << "\n"
              << "hess(p)(0)=" << p.gradient().gradient().evaluate( pts ) << "\n"
              << "hess(p1) = " << p1.toSet( true ).gradient().gradient().evaluate( pts ) << "\n";

    ublas::matrix<double,ublas::row_major> m = p.gradient().gradient().evaluate( pts );
    std::cout << "m=" << m << "\n";
    const int I = p.coeff().size1();
    const int nDim = 2;
    const int Q = m.size2();
    boost::multi_array<double,4> hessian( boost::extents[I][nDim][nDim][Q] );
    typedef boost::multi_array<double,4>::index index;

    for ( index i = 0; i < I; ++i )
        for ( index j = 0; j < nDim; ++j )
            for ( index k = 0; k < nDim; ++k )
                for ( index q = 0; q < Q; ++q )
                {
#if 1
                    std::cout << "[precompute] hessian["
                              << i << "]["
                              << j << "]["
                              << k << "]["
                              << q << "]["
                              << nDim*nDim*I*( nDim*k+j )+nDim*nDim*i+nDim*j+k << "]="
                              << "\n";
#endif // 0

                    hessian[i][j][k][q] = m( nDim*nDim*I*( nDim*k+j )+nDim*nDim*i+nDim*j+k, q );

                    std::cout << hessian[i][j][k][q] << "\n";
                }


}

BOOST_AUTO_TEST_CASE( test_polynomial3D )
{


    ublas::matrix<double,ublas::column_major> pts( 3, 4 );

    //1
    pts( 0,0 ) = -1;
    pts( 1,0 ) = -1;
    pts( 2,0 ) = -1;
    //2
    pts( 0,1 ) =  1;
    pts( 1,1 ) = -1;
    pts( 2,1 ) = -1;
    //3
    pts( 0,2 ) = -1;
    pts( 1,2 ) =  1;
    pts( 2,2 ) = -1;
    // 4
    pts( 0,3 ) =  0;
    pts( 1,3 ) =  0;
    pts( 2,3 ) =  0;
    detail::OrthonormalPolynomialSet<3,3,3,Scalar> ps;
    Polynomial<detail::OrthonormalPolynomialSet<3,3,3,Scalar>, Scalar> p = project( ps, x2(), IM<3,5>() );
    std::cout << "p=" << p.toSet().coeff() << "\n"
              << "grad(p)(0)=" << p.toSet().gradient().evaluate( pts ) << "\n";
    std::cout << "hess(p)=" << p.toSet().gradient().gradient().coeff() << "\n"
              << "hess(p)(0)=" << p.toSet().gradient().gradient().evaluate( pts ) << "\n";


}

#if 0
#if defined(USE_TEST)
test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv*/ )
{

    test_suite* test = BOOST_TEST_SUITE( "Polynomial test suite" );


    test->add( BOOST_TEST_CASE( test_polynomial1D ) );
    test->add( BOOST_TEST_CASE( test_polynomial2D ) );
    //test->add( BOOST_TEST_CASE( test_polynomial3D ) );

    return test;
}
#else
int main( int argc, char** argv )
{
    test_polynomial1D();
    //    test_polynomial2D();
    //test_polynomial3D();
}
#endif
#endif
