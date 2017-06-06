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
#define BOOST_TEST_MODULE function space on range testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/concatenate.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( space_on_range )

BOOST_AUTO_TEST_CASE( test_2d )
{
    using namespace Feel;
    typedef Mesh<Simplex<2,1>> mesh_type;
    auto mesh = loadMesh(_mesh=new mesh_type );

    typedef FunctionSpace<mesh_type,bases<Lagrange<2,Scalar> > > space_type;

    auto r1 = elements(mesh, pow(Px()-0.4,2)+pow(Py()-0.5,2) < pow(cst(0.23),2) );
    auto r2 = elements(mesh, pow(Px()-0.7,2)+pow(Py()-0.5,2) < pow(cst(0.15),2) );
    auto therange = concatenate(r1,r2);

    auto VhPS = space_type::New(_mesh=mesh,_range=therange);
    BOOST_TEST_MESSAGE( "VhPS->nDof() : " << VhPS->nDof() );
    auto VhFS = space_type::New(_mesh=mesh);
    BOOST_TEST_MESSAGE( "VhFS->nDof() : " << VhFS->nDof());

    // test element of space
    auto chiShapeFS = VhFS->element();
    chiShapeFS.setConstant(5.);
    chiShapeFS.on(_range=therange,_expr=cst(3.));
    auto chiShapePS = VhPS->element();
    chiShapePS.setConstant(3.);
    double err = normL2(_range=therange,_expr=idv(chiShapeFS)-idv(chiShapePS));
    BOOST_CHECK_SMALL( err,1e-12 );
    chiShapeFS.on(_range=therange,_expr=Px()*Py());
    chiShapePS.on(_range=therange,_expr=Px()*Py());
    err = normL2(_range=therange,_expr=idv(chiShapeFS)-idv(chiShapePS));
    BOOST_CHECK_SMALL( err,1e-12 );
    std::cout << "err="<<err<<"\n";
    double chiFSnorm2 = chiShapeFS.l2Norm();
    double chiPSnorm2 = chiShapePS.l2Norm();
    BOOST_CHECK( chiFSnorm2 > chiPSnorm2 );

    // test a laplacian solve
    auto u = VhPS->element("u");
    auto v = VhPS->element("v");
    BOOST_CHECK( VhPS->dof()->hasMeshSupport() );
    VhPS->dof()->meshSupport()->updateBoundaryFaces();
    auto myboundaryfaces = VhPS->dof()->meshSupport()->rangeBoundaryFaces();
    auto l = form1( _test=VhPS );
    l = integrate(_range=therange,
                  _expr=id(v));
    auto a = form2( _trial=VhPS, _test=VhPS);
    a = integrate(_range=therange,
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=myboundaryfaces, _rhs=l, _element=u, _expr=cst(0.) );
    a.solve(_rhs=l,_solution=u);

    // export results
    auto e = exporter( _mesh=mesh,_name="test2d" );
    e->addRegions();
    e->add( "u", u );
    e->add( "chi-partial-support", chiShapeFS );
    e->add( "chi-full-support", chiShapePS );
    e->save();

}

BOOST_AUTO_TEST_SUITE_END()

