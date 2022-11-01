/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 21 Dec 2019

 Copyright (C) 2019 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#define BOOST_TEST_MODULE test_convex
#include <boost/mp11/list.hpp>
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelmesh/simplex.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/makemesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/geometricdata.hpp>


using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( convex_suite )

    
BOOST_AUTO_TEST_CASE( t0 )
{

    Simplex<1,1,1> s1;
    BOOST_CHECK_EQUAL( s1.numberOfVertices(), 2 );
    Simplex<2,1,2> s2;
    BOOST_CHECK_EQUAL( s2.numberOfVertices(), 3 );
    Simplex<3,1,3> s3;
    BOOST_CHECK_EQUAL( s3.numberOfVertices(), 4 );
    Simplex<1,1,2> s4;
    BOOST_CHECK_EQUAL( s4.numberOfVertices(), 2 );
    Simplex<2,1,3> s5;
    BOOST_CHECK_EQUAL( s5.numberOfVertices(), 3 );
    Simplex<1,1,3> s6;
    BOOST_CHECK_EQUAL( s6.numberOfVertices(), 2 );
}
    
BOOST_AUTO_TEST_CASE( t01 )
{

    Simplex<1,Dynamic,1> s1{1};
    BOOST_CHECK_EQUAL( s1.numberOfVertices(), 2 );
    Simplex<2,Dynamic,2> s2{1};
    BOOST_CHECK_EQUAL( s2.numberOfVertices(), 3 );
    Simplex<3,Dynamic,3> s3{1};
    BOOST_CHECK_EQUAL( s3.numberOfVertices(), 4 );
    BOOST_CHECK_EQUAL( s3.numberOfPoints(), 4 );
    BOOST_CHECK_EQUAL( s3.numberOfPointsAtCompileTime(), -1 );

    Simplex<3,Dynamic,3> s30{0};
    Simplex<3,0,3> s30_ct;
    BOOST_CHECK_EQUAL( s30.numberOfVertices(), 4 );
    BOOST_CHECK_EQUAL( s30.numberOfPoints(), 1 );
    BOOST_CHECK_EQUAL( s30.numberOfPoints(), s30_ct.numberOfPoints() );
    BOOST_CHECK_EQUAL( s30.numberOfPointsAtCompileTime(), -1 );

    Simplex<3,Dynamic,3> s32{2};
    BOOST_CHECK_EQUAL( s32.numberOfVertices(), 4 );
    BOOST_CHECK_EQUAL( s32.numberOfPoints(), 10 );
    BOOST_CHECK_EQUAL( s32.numberOfPointsAtCompileTime(), -1 );

    Simplex<1,Dynamic,2> s4{1};
    BOOST_CHECK_EQUAL( s4.numberOfVertices(), 2 );
    Simplex<2,Dynamic,3> s5{1};
    BOOST_CHECK_EQUAL( s5.numberOfVertices(), 3 );
    Simplex<1,Dynamic,3> s6{1};
    BOOST_CHECK_EQUAL( s6.numberOfVertices(), 2 );
}

using order_t = boost::mp11::mp_list_c<int, 0,1,2,3,4,5,6,10>;
BOOST_AUTO_TEST_CASE_TEMPLATE( t02, T, order_t )
{
    PointSetEquiSpaced<Simplex<3,Dynamic,3>,Dynamic,double> pts{T::value,0};
    PointSetEquiSpaced<Simplex<3,1,3>,T::value,double> pts_ct;
    BOOST_CHECK_EQUAL(pts.numberOfPoints(),pts_ct.numberOfPoints());
    BOOST_CHECK_SMALL((em(pts.points())-em(pts_ct.points())).norm(),1e-12);

}
BOOST_AUTO_TEST_CASE( t1 )
{
#if 0
    Mesh<Simplex<1>> m1;
    Mesh<Simplex<2>> m2;
    Mesh<Simplex<3>> m3;
    BOOST_CHECK_EQUAL( m3.nDim, 3 );
#endif
}

using dim_t = boost::mp11::mp_list_c<int, 2,3>;
BOOST_AUTO_TEST_CASE_TEMPLATE( t2, T, dim_t )
{
#if 0    
    auto mesh = makeMesh<Simplex<T::value>>();
    if constexpr ( T::value == 3 )
        mesh = unitCube();
    else
        mesh = unitSquare();
    auto Xh = Pch<2>( mesh );
    auto v = Xh->element(Px());
#endif    
}

BOOST_AUTO_TEST_CASE( t3 )
{
}
BOOST_AUTO_TEST_CASE( t4 )
{

}
BOOST_AUTO_TEST_SUITE_END()


