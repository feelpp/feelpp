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
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelmesh/simplex.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/makemesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/thch.hpp>
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

BOOST_AUTO_TEST_CASE( t1 )
{
    Mesh<Simplex<1>> m1;
    Mesh<Simplex<2>> m2;
    Mesh<Simplex<3>> m3;
    BOOST_CHECK_EQUAL( m3.nDim, 3 );
}

typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( t2, T, dim_types )
{
    auto mesh = makeMesh<Simplex<T::value>>();
    if constexpr ( T::value == 3 )
        mesh = unitCube();
    else
        mesh = unitSquare();
    tic();
    auto Xh = Pch<2>( mesh );
    toc("Xh built", true );
    auto v = Xh->element(Px());
    BOOST_CHECK_EQUAL( Xh->basisOrder()[0], 2 );
    BOOST_CHECK( Xh->basisName() == "lagrange" );
    BOOST_TEST_MESSAGE(  "nDof=" << Xh->nDof() << " family:" << Xh->basisName() << " Order: " << Xh->basisOrder()[0] );
    auto Vh = THch<1>( mesh );
    BOOST_CHECK_EQUAL( Vh->basisOrder()[0], 2 );
    BOOST_CHECK_EQUAL( Vh->basisOrder()[1], 1 );
    BOOST_CHECK( Vh->basisName() == "lagrange_lagrange" );
    BOOST_TEST_MESSAGE(  "nDof=" << Vh->nDof() << " family:" << Vh->basisName() << " Order: " << Vh->basisOrder() );
}

BOOST_AUTO_TEST_CASE( t3 )
{
    
}
BOOST_AUTO_TEST_CASE( t4 )
{

}
BOOST_AUTO_TEST_SUITE_END()


