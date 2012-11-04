/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-10
 */
#define DO_TEST 1
// Boost.Test
#define BOOST_TEST_MAIN
#include <testsuite/testsuite.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/geoelement.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_CASE( test_connectivity1 )
{
    using namespace Feel;
    GeoND<2,Simplex<2> > c1;
    BOOST_CHECK_EQUAL( c1.e2p( 0,0 ), 1 );
    BOOST_CHECK_EQUAL( c1.e2p( 0,1 ), 2 );
    BOOST_CHECK_EQUAL( c1.e2p( 1,0 ), 2 );
    BOOST_CHECK_EQUAL( c1.e2p( 1,1 ), 0 );
    BOOST_CHECK_EQUAL( c1.e2p( 2,0 ), 0 );
    BOOST_CHECK_EQUAL( c1.e2p( 2,1 ), 1 );
}
BOOST_AUTO_TEST_CASE( test_connectivity2 )
{
    using namespace Feel;
    GeoND<3,Simplex<3> > c1;
    BOOST_CHECK_EQUAL( c1.e2p( 0,0 ), 1 );
    BOOST_CHECK_EQUAL( c1.e2p( 0,1 ), 2 );
    BOOST_CHECK_EQUAL( c1.e2p( 1,0 ), 2 );
    BOOST_CHECK_EQUAL( c1.e2p( 1,1 ), 0 );
    BOOST_CHECK_EQUAL( c1.e2p( 2,0 ), 0 );
    BOOST_CHECK_EQUAL( c1.e2p( 2,1 ), 1 );
    BOOST_CHECK_EQUAL( c1.e2p( 3,0 ), 0 );
    BOOST_CHECK_EQUAL( c1.e2p( 3,1 ), 3 );
    BOOST_CHECK_EQUAL( c1.e2p( 4,0 ), 1 );
    BOOST_CHECK_EQUAL( c1.e2p( 4,1 ), 3 );
    BOOST_CHECK_EQUAL( c1.e2p( 5,0 ), 2 );
    BOOST_CHECK_EQUAL( c1.e2p( 5,1 ), 3 );
}
template<int Dim, int Order>
void
test_entity_isin()
{
    using namespace Feel;

    Reference<Simplex<Dim, Order, Dim>, Dim, Order, Dim > ref;

    typename node<double>::type n( Dim );
    std::fill( n.begin(), n.end(), double( 2 ) );

    std::cout << "n=" << n << " isin = " << boost::get<0>( ref.isIn( n ) ) << "\n";

    if ( Dim == 3 )
    {
        n[0]=-0.9;
        n[1]=-0.9;
        n[2]=-1;
        std::cout << "n=" << n << " isin = " << boost::get<0>( ref.isIn( n ) ) << "\n";
    }
}
template<int Dim, int Order>
void
test_entity_range()
{
    using namespace Feel;

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

    for ( int d = 0; d <= rs.topologicalDimension(); ++d )
        for ( int e = rs.entityRange( d ).begin();
                e < rs.entityRange( d ).end();
                ++e )
        {
            std::cerr << "(" << d << "," << e << "): \n";
            Feel::matrix_node<double>::type G(  rs.makePoints( d, e ) );
            std::cerr  << G << "\n";
        }
}
BOOST_AUTO_TEST_CASE( test_entity_isin_ )
{
    test_entity_isin<2,1>();
    test_entity_isin<3,1>();

}

typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<4>,boost::mpl::int_<5> > test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_refelem, T, test_types )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "o- TestRefElem<" << T::value << ">\n" );
    Reference<Simplex<3, T::value, 3>, 3, T::value, 3, double> tetraref;
    toPython( tetraref );
}
BOOST_AUTO_TEST_CASE( test_simplex_ref )
{
    using namespace Feel;

    // interval
    Reference<Simplex<1, 1, 1>, 1, 1, 1, double> intervalref;
    ublas::vector<double> G( 1 );
    G( 0 ) = 0;
    BOOST_CHECK_SMALL( ublas::norm_2( intervalref.barycenter()-G ), 1e-14 );
    G( 0 ) = -1;
    BOOST_CHECK_SMALL( ublas::norm_2( intervalref.faceBarycenter( 0 )-G ), 1e-14 );
    G( 0 ) = 1;
    BOOST_CHECK_SMALL( ublas::norm_2( intervalref.faceBarycenter( 1 )-G ), 1e-14 );

    // triangle
    Reference<Simplex<2, 1, 2>, 2, 1, 2, double> triangleref;
    ublas::vector<double> G2( 2 );
    G2( 0 ) = -1./3.;
    G2( 1 ) = -1./3.;
    BOOST_CHECK_SMALL( ublas::norm_2( triangleref.barycenter()-G2 ), 1e-14 );
    G2( 0 ) = 0.;
    G2( 1 ) = 0;
    BOOST_CHECK_SMALL( ublas::norm_2( triangleref.faceBarycenter( 0 )-G2 ), 1e-14 );
    G2( 0 ) = -1.;
    G2( 1 ) = 0;
    BOOST_CHECK_SMALL( ublas::norm_2( triangleref.faceBarycenter( 1 )-G2 ), 1e-14 );
    G2( 0 ) = 0.;
    G2( 1 ) = -1;
    BOOST_CHECK_SMALL( ublas::norm_2( triangleref.faceBarycenter( 2 )-G2 ), 1e-14 );

    // tetra
    Reference<Simplex<3, 1, 3>, 3, 1, 3, double> tetraref;
    ublas::vector<double> G3( 3 );
    G3( 0 ) = -1./2.;
    G3( 1 ) = -1./2.;
    G3( 2 ) = -1./2.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetraref.barycenter()-G3 ), 1e-14 );
    G3( 0 ) = -1./3.;
    G3( 1 ) = -1./3.;
    G3( 2 ) = -1./3.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetraref.faceBarycenter( 0 )-G3 ), 1e-14 );
    G3( 0 ) = -1.;
    G3( 1 ) = -1./3.;
    G3( 2 ) = -1./3.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetraref.faceBarycenter( 1 )-G3 ), 1e-14 );
    G3( 0 ) = -1./3.;
    G3( 1 ) = -1.;
    G3( 2 ) = -1./3.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetraref.faceBarycenter( 2 )-G3 ), 1e-14 );
    G3( 0 ) = -1./3.;
    G3( 1 ) = -1./3.;
    G3( 2 ) = -1.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetraref.faceBarycenter( 3 )-G3 ), 1e-14 );
}
BOOST_AUTO_TEST_CASE( test_interval )
{
    using namespace Feel;

    typedef GeoND<1,Simplex<1, 1, 1> >::point_type point_type;
    // interval
    GeoND<1,Simplex<1, 1, 1> > interval;
    point_type V1( 1 );
    V1( 0 )=1;
    point_type V2( 1 );
    V2( 0 )=3;
    interval.setPoint( 0, V1 );
    interval.setPoint( 1, V2 );
    ublas::vector<double> G1( 1 );
    G1( 0 )=2;
    interval.update();
    std::cout << "[interval] barycenter = " << interval.barycenter() << "\n";
    BOOST_CHECK_SMALL( ublas::norm_2( interval.barycenter()-G1 ), 1e-14 );
    BOOST_CHECK_CLOSE( interval.measure(), 2, 1e-14 );
    G1( 0 )=-1;
    BOOST_CHECK_SMALL( ublas::norm_2( interval.normal( 0 )-G1 ), 1e-14 );
    G1( 0 )=1;
    BOOST_CHECK_SMALL( ublas::norm_2( interval.normal( 1 )-G1 ), 1e-14 );

    BOOST_CHECK_SMALL( ( double )ublas::norm_frobenius( interval.G() - interval.vertices() ), 1e-14 );

}

BOOST_AUTO_TEST_CASE( test_triangle )
{
    using namespace Feel;

    typedef GeoND<2,Simplex<2, 1, 2> >::point_type point_type;
    // interval
    GeoND<2,Simplex<2, 1, 2> > tria;
    point_type V1;
    V1( 0 )=1;
    V1( 1 )=0;
    point_type V2;
    V2( 0 )=3;
    V2( 1 )=0;
    point_type V3;
    V3( 0 )=2;
    V3( 1 )=1;
    tria.setPoint( 0, V1 );
    tria.setPoint( 1, V2 );
    tria.setPoint( 2, V3 );
    ublas::vector<double> G1( 2 );
    G1( 0 )=2;
    G1( 1 )=1./3.;
    tria.update();
    std::cout << "[tria] barycenter = " << tria.barycenter() << "\n";
    BOOST_CHECK_SMALL( ublas::norm_2( tria.barycenter()-G1 ), 1e-14 );
    BOOST_CHECK_CLOSE( tria.measure(), 1, 1e-14 );
    // check normals
    G1( 0 )=1./math::sqrt( 2. );
    G1( 1 )=1./math::sqrt( 2. );
    BOOST_CHECK_SMALL( ublas::norm_2( tria.normal( 0 )-G1 ), 1e-14 );
    G1( 0 )=-1./math::sqrt( 2. );
    G1( 1 )=1./math::sqrt( 2. );
    BOOST_CHECK_SMALL( ublas::norm_2( tria.normal( 1 )-G1 ), 1e-14 );
    G1( 0 )=0;
    G1( 1 )=-1.;
    BOOST_CHECK_SMALL( ublas::norm_2( tria.normal( 2 )-G1 ), 1e-14 );
    // check face measures
    BOOST_CHECK_CLOSE( tria.faceMeasure( 0 ), math::sqrt( 2. ), 2e-14 );
    BOOST_CHECK_CLOSE( tria.faceMeasure( 1 ), math::sqrt( 2. ), 1e-14 );
    BOOST_CHECK_CLOSE( tria.faceMeasure( 2 ), 2, 1e-14 );

    BOOST_CHECK_SMALL( ( double )ublas::norm_frobenius( tria.G() - tria.vertices() ), 1e-14 );
}
BOOST_AUTO_TEST_CASE( test_tetra )
{
    using namespace Feel;

    typedef GeoND<3,Simplex<3, 1, 3> >::point_type point_type;
    // interval
    GeoND<3,Simplex<3, 1, 3> > tetra;
    point_type V1;
    V1( 0 )=1;
    V1( 1 )=0;
    V1( 2 )=0;
    point_type V2;
    V2( 0 )=3;
    V2( 1 )=0;
    V2( 2 )=0;
    point_type V3;
    V3( 0 )=2;
    V3( 1 )=1;
    V3( 2 )=0;
    point_type V4;
    V4( 0 )=2;
    V4( 1 )=1;
    V4( 2 )=1;
    tetra.setPoint( 0, V1 );
    tetra.setPoint( 1, V2 );
    tetra.setPoint( 2, V3 );
    tetra.setPoint( 3, V4 );
    ublas::vector<double> G1( 3 );
    G1( 0 )=2;
    G1( 1 )=1./2.;
    G1( 2 )=1./4.;
    tetra.update();
    std::cout << "[tetra] barycenter = " << tetra.barycenter() << "\n";
    BOOST_CHECK_SMALL( ublas::norm_2( tetra.barycenter()-G1 ), 1e-14 );
    BOOST_CHECK_CLOSE( tetra.measure(), 1./3., 1e-13 );

    BOOST_CHECK_SMALL( ( double )ublas::norm_frobenius( tetra.G() - tetra.vertices() ), 1e-14 );
#if 0
    // check normals
    G1 = cross( V3-V2, V4-V2 );
    std::cout << "G1 = " << G1 << "\n";

    BOOST_CHECK_SMALL( ublas::norm_2( tetra.normal( 0 )-G1 ), 1e-14 );
    G1( 0 )=-1./math::sqrt( 2. );
    G1( 1 )=1./math::sqrt( 2. );
    BOOST_CHECK_SMALL( ublas::norm_2( tetra.normal( 1 )-G1 ), 1e-14 );
    G1( 0 )=0;
    G1( 1 )=-1.;
    BOOST_CHECK_SMALL( ublas::norm_2( tetra.normal( 2 )-G1 ), 1e-14 );
    // check face measures
    BOOST_CHECK_CLOSE( tetra.faceMeasure( 0 ), math::sqrt( 2. ), 2e-14 );
    BOOST_CHECK_CLOSE( tetra.faceMeasure( 1 ), math::sqrt( 2. ), 1e-14 );
    BOOST_CHECK_CLOSE( tetra.faceMeasure( 2 ), 2, 1e-14 );
#endif
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
    Feel::Assert::setLog( "assertions.log" );
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
