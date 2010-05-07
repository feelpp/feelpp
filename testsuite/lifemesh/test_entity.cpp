/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-08-10

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2009 Université de Grenoble 1 (Joseph Fourier)

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
   \file test_entity.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-08-10
 */
#define DO_TEST 1
// Boost.Test
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <life/lifecore/life.hpp>
#include <life/lifemesh/simplex.hpp>
#include <life/lifemesh/simplexproduct.hpp>
#include <life/lifemesh/refentity.hpp>

template<int Dim, int Order>
void
test_entity_isin()
{
    using namespace Life;

    Reference<Simplex<Dim, Order, Dim>, Dim, Order, Dim > ref;

    typename node<double>::type n( Dim );
    std::fill( n.begin(), n.end(), double( 2 ) );

    std::cout << "n=" << n << " isin = " << boost::get<0>( ref.isIn( n ) ) << "\n";

    if ( Dim == 3 )
        {
            n[0]=-0.9;n[1]=-0.9;n[2]=-1;
            std::cout << "n=" << n << " isin = " << boost::get<0>( ref.isIn( n ) ) << "\n";
        }
}
template<int Dim, int Order>
void
test_entity_range()
{
    using namespace Life;

    Simplex<Dim, Order> s;
    EntityRange<Simplex<Dim, Order> > r;
    r.setTopologicalDimension( 0 );
    //BOOST_CHECK( r.end() == s.numPoints );
    r.setTopologicalDimension( 1 );
    //BOOST_CHECK( r.end() == s.numEdges );
    r.setTopologicalDimension( 2 );
    //BOOST_CHECK( r.end() == s.numGeometricFaces );
    r.setTopologicalDimension( 3 );
    //BOOST_CHECK( r.end() == s.numVolumes );

    Reference<Simplex<Dim, Order,Dim>,Dim,Order,Dim> rs;
    std::cout << "S = " << rs << "\n";
    toPython( rs );
    for ( int d = 0;d <= rs.topologicalDimension();++d )
        for ( int e = rs.entityRange( d ).begin();
              e < rs.entityRange( d ).end();
              ++e )
        {
            std::cerr << "(" << d << "," << e << "): \n";
            Life::matrix_node<double>::type G(  rs.makePoints( d, e ) );
            std::cerr  << G << "\n";
        }
}
BOOST_AUTO_TEST_CASE( test_entity_isin_ )
{
    test_entity_isin<2,1>();
    test_entity_isin<3,1>();

}

BOOST_AUTO_TEST_CASE( test_simplex_ref )
{
    using namespace Life;

    // interval
    Reference<Simplex<1, 1, 1>, 1, 1, 1, double> intervalref;
    ublas::vector<double> G( 1 ); G( 0 ) = 0;
    BOOST_CHECK_SMALL( ublas::norm_2( intervalref.barycenter()-G ), 1e-14 );
    G( 0 ) = -1;
    BOOST_CHECK_SMALL( ublas::norm_2( intervalref.faceBarycenter(0)-G ), 1e-14 );
    G( 0 ) = 1;
    BOOST_CHECK_SMALL( ublas::norm_2( intervalref.faceBarycenter(1)-G ), 1e-14 );

    // triangle
    Reference<Simplex<2, 1, 2>, 2, 1, 2, double> triangleref;
    ublas::vector<double> G2( 2 );
    G2( 0 ) = -1./3.;G2( 1 ) = -1./3.;
    BOOST_CHECK_SMALL( ublas::norm_2( triangleref.barycenter()-G2 ), 1e-14 );
    G2( 0 ) = 0.;G2( 1 ) = 0;
    BOOST_CHECK_SMALL( ublas::norm_2( triangleref.faceBarycenter(0)-G2 ), 1e-14 );
    G2( 0 ) = -1.;G2( 1 ) = 0;
    BOOST_CHECK_SMALL( ublas::norm_2( triangleref.faceBarycenter(1)-G2 ), 1e-14 );
    G2( 0 ) = 0.;G2( 1 ) = -1;
    BOOST_CHECK_SMALL( ublas::norm_2( triangleref.faceBarycenter(2)-G2 ), 1e-14 );

    // tetra
    Reference<Simplex<3, 1, 3>, 3, 1, 3, double> tetraref;
    ublas::vector<double> G3( 3 );
    G3( 0 ) = -1./2.;G3( 1 ) = -1./2.;G3( 2 ) = -1./2.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetraref.barycenter()-G3 ), 1e-14 );
    G3( 0 ) = -1./3.;G3( 1 ) = -1./3.;G3( 2 ) = -1./3.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetraref.faceBarycenter(0)-G3 ), 1e-14 );
    G3( 0 ) = -1.;G3( 1 ) = -1./3.;G3( 2 ) = -1./3.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetraref.faceBarycenter(1)-G3 ), 1e-14 );
    G3( 0 ) = -1./3.;G3( 1 ) = -1.;G3( 2 ) = -1./3.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetraref.faceBarycenter(2)-G3 ), 1e-14 );
    G3( 0 ) = -1./3.;G3( 1 ) = -1./3.;G3( 2 ) = -1.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetraref.faceBarycenter(3)-G3 ), 1e-14 );
}

BOOST_AUTO_TEST_CASE( test_entity_range_ )
{
    //test_entity_range<1,1>();
}
#if 0
#if DO_TEST
test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
    Life::Assert::setLog( "assertions.log");
    test_suite* test = BOOST_TEST_SUITE( "2D Generic finite element solver test suite" );

    test->add( BOOST_TEST_CASE( ( test_entity_isin<2,1> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_isin<3,1> ) ) );
#if 0
    test->add( BOOST_TEST_CASE( ( test_entity_range<1,1> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<1,2> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<1,3> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<1,4> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<2,1> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<2,2> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<2,3> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<2,4> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<3,1> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<3,2> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<3,3> ) ) );
    test->add( BOOST_TEST_CASE( ( test_entity_range<3,4> ) ) );
#endif
    return test;
}
#else
int main()
{
    test_entity_range<1,1>();
    test_entity_range<1,2>();
    test_entity_range<2,1>();
    test_entity_range<2,2>();
    test_entity_range<2,3>();
    test_entity_range<3,1>();
    test_entity_range<3,2>();
}
#endif
#endif
