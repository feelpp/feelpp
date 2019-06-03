/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-30

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
/**
   \file test_functionals.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-01-30
 */
// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE functionals
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/odh.hpp>
#include <feel/feeldiscr/pchv.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( functionals )

BOOST_AUTO_TEST_CASE( test_projection_hdiv_rt )
{
    BOOST_TEST_MESSAGE( "test Directional component integration" );

    using namespace Feel;

    // only one element
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>>, _filename="ReferenceTriangle.geo" );

    // RTk
    auto RTh = Dh<0>( mesh );
    auto urt = RTh->element();
    urt.zero();
    BOOST_TEST_MESSAGE( "Rth defined, dimension: " << RTh->nLocalDof() );
    // Lagrange
    auto Xh = Pchv<1>( mesh );
    auto ul = Xh->element();
    ul = vf::project( _space=Xh, _range=elements(mesh), _expr=P() );
    //ul.printMatlab( "ul.m" );
    BOOST_TEST_MESSAGE( "Xh defined, dimension: " << Xh->nLocalDof() );
    BOOST_TEST_MESSAGE( "n elements: " << nelements( elements(mesh)  ) );
    for( auto const& elementRange : elements(mesh) )
    {
        auto const& element = boost::unwrap_ref( elementRange );
        LOG(INFO) << "element : " << element.id();
        auto mesh_element = createSubmesh( _mesh=mesh, _range=idedelements(mesh,element.id()), _worldcomm=Environment::worldCommSeqPtr() );
        BOOST_TEST_MESSAGE( "n elements of extracted element " << element.id() << " : " << nelements( elements(mesh_element)  ) );
        auto mesh_face = createSubmesh( _mesh=mesh_element, _range=boundaryfaces( mesh_element ), _worldcomm=Environment::worldCommSeqPtr() );
        BOOST_TEST_MESSAGE( "n boundary faces in extracted element " << element.id() << " : " << nelements( elements(mesh_face)  ) );
        auto Ph = Odh<0>( mesh_face );
        BOOST_TEST_MESSAGE( "dimension of Ph : " << Ph->nLocalDof() );
        auto p = Ph->element();
        auto l = form1( _test=Ph );
        //l = integrate( _range=boundaryfaces(mesh_element), _expr=trans(print(idv(ul),"ul"))*print(N(),"n")*print(id(p),"p" ) );
        l = integrate( _range=boundaryfaces(mesh_element), _expr=trans(idv(ul))*N()*id(p) );
        //l.vector().printMatlab( "l.m" );
        for( auto const& dof: RTh->dof()->localDof( element.id() ) )
        {
            LOG(INFO) << "dof[ " << element.id() << ", " << dof.first.localDof() << "]=" << dof.second.index();
            urt[dof.second.index()] = l.vector()(dof.first.localDof());
            LOG(INFO) << "value[dof.second.index()] = " << urt[dof.second.index()];
        }
    }
    //urt.printMatlab( "urt.m" );

    BOOST_TEST_MESSAGE( "test Directional component integration done" );
}

BOOST_AUTO_TEST_SUITE_END()
