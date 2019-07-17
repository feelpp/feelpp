/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 27 Sep 2014

   Copyright (C) 2014-2016 Feel++ Consortium

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
// Boost.Test
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE mesh iterator operations testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <feel/feelcore/testsuite.hpp>


#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/filters.hpp>
#include <feel/feelmesh/concatenate.hpp>
#include <feel/feelmesh/intersect.hpp>
#include <feel/feelvf/measure.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( mesh_filters_operations )

BOOST_AUTO_TEST_CASE( test_concatenate_elements )
{
    using namespace Feel;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2,1>>);
    BOOST_CHECK( mesh->hasMarkers( {"left","right"} ) );

    auto e3 = concatenate( markedelements(mesh,"left"), markedelements(mesh,"right") );
    auto submesh1 = createSubmesh( mesh, e3 );
    BOOST_CHECK_EQUAL( nelements(elements(submesh1)), nelements(elements(mesh)) );

    auto submesh2 = createSubmesh( mesh, markedelements(mesh,{"left","right"}) );
    BOOST_CHECK_EQUAL( nelements(elements(submesh2)), nelements(elements(mesh)) );

    auto e4 = concatenate( markedelements(mesh,"left"), elements(mesh) );
    BOOST_CHECK_EQUAL( nelements(elements(mesh)), nelements(e4) );
    auto e5 = concatenate( markedelements(mesh,"left"), markedelements(mesh,"right"), elements(mesh) );
    BOOST_CHECK_EQUAL( nelements(elements(mesh)), nelements(e5) );
}

BOOST_AUTO_TEST_CASE( test_intersect_elements )
{
    using namespace Feel;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2,1>>);
    auto rLeft = markedelements(mesh,"left");
    auto rRight = markedelements(mesh,"right");
    auto rAll = elements(mesh);
    auto rBe = boundaryelements(mesh);
    auto intersect1 = intersect( rLeft,rBe );
    auto intersect2 = intersect( rRight,rBe,rAll );
    BOOST_CHECK_EQUAL( nelements(intersect1)+nelements(intersect2), nelements(rBe) );
    auto submesh1 = createSubmesh( mesh, intersect1 );
    auto submesh2 = createSubmesh( mesh, intersect2 );
    double measSubmesh1 = measure(_range=elements(submesh1) );
    double measSubmesh2 = measure(_range=elements(submesh2) );
    double measMeshBe = measure(_range=rBe );
    BOOST_CHECK_EQUAL( measSubmesh1 + measSubmesh2, measMeshBe );
}


BOOST_AUTO_TEST_SUITE_END()

