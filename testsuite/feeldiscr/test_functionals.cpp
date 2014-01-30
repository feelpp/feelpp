/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-30

  Copyright (C) 2014 Feel++ Consortium

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
/**
   \file test_functionals.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-01-30
 */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE functionals
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/pchv.hpp>

#if defined( USE_BOOST_TEST )

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( functionals )

BOOST_AUTO_TEST_CASE( test_projection_hdiv_rt )
{
    BOOST_TEST_MESSAGE( "test Directional component integration" );

    // only one element
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>>, _filename="ReferenceTriangle.geo" );

    // RTk
    auto Rth = Dh<1>( mesh );
    auto urt = Xh->element();
    urt.zero();

    // Lagrange
    auto Xh = Pchv<1>( mesh );
    auto ul = Xh->element();
    ul = vf::project( _space=Xh, _range=elements(mesh), _expr=P() );

    for( auto const& element : elements(mesh) )
    {
        auto mesh_element = createSubmesh( mesh, idedelement(element.id()), Environment::worldCommSeq() );
        auto mesh_face = createSubmesh( mesh_element, boundaryfaces( mesh ), Environment::worldCommSeq() );
        auto P = Odh<0>( mesh_face );
        auto p = P->element();
        auto l = form1( _test=P );
        l = integrate( _range=elements(mesh_face), _expr=trans(idv(u))*N()*id(p) );
        urt.add( l.vector() );
    }


    BOOST_TEST_MESSAGE( "test Directional component integration done" );
}

BOOST_AUTO_TEST_SUITE_END()
#else

int
main( int argc, char* argv[] )
{
    Feel::Environment env( argc,argv,
                           makeAbout(), makeOptions() );


}

#endif
