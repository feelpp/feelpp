/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannese@feelpp.org>
       Date: 2016-06-30

  Copyright (C) 2016 Feel++ Consortium

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
   \file test_vf_operators.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-06-30
 */

#define BOOST_TEST_MODULE vf_operators testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/geotool.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_vf_operators )

BOOST_AUTO_TEST_CASE( test_vf_operators_grad_dn )
{
    using namespace Feel;
    BOOST_TEST_MESSAGE( "test_operators_grad_dn" );

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Rectangle R( doption(_name="gmsh.hsize"),"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="Boundary1",_marker3=true);
    R.setMarker(_type="line",_name="Boundary2",_marker1=true,_marker2=true,_marker4=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new Mesh<Simplex<2> >,
                             _name="domainRectangle" );
    // scalar
    auto Vh = Pch<5>( mesh );
    auto uScalar = Vh->element( Px()*Px()*Px()*Py()*Py() );
    auto gradScalarExact = trans( vec( 3*Px()*Px()*Py()*Py(), Px()*Px()*Px()*2*Py() ) );
    double resScalarExact = integrate(_range=boundaryfaces(mesh),_expr=gradScalarExact*N() ).evaluate()(0,0);
    double resScalarGradv = integrate(_range=boundaryfaces(mesh),_expr=gradv(uScalar)*N() ).evaluate()(0,0);
    double resScalarDnv = integrate(_range=boundaryfaces(mesh),_expr=dnv(uScalar) ).evaluate()(0,0);
    std::cout << "resScalarExact " << resScalarExact << "\n";
    std::cout << "resScalarGradv " << resScalarGradv << "\n";
    std::cout << "resScalarDnv " << resScalarDnv << "\n";
    BOOST_CHECK_CLOSE( resScalarGradv, resScalarExact, 1e-10 );
    BOOST_CHECK_CLOSE( resScalarDnv, resScalarExact, 1e-10 );

    // vectorial
    auto Wh = Pchv<5>( mesh );
    auto uVectorial = Wh->element(vec(Px()*Px()*Px()*Py()*Py(), Px()*Px() ) );
    auto gradVectorialExact = mat<2,2>( 3*Px()*Px()*Py()*Py(), Px()*Px()*Px()*2*Py(),
                                        2*Px(), cst(0.) );
    double resVectorialExact = integrate(_range=boundaryfaces(mesh),_expr= inner(gradVectorialExact*N(),gradVectorialExact*N()) ).evaluate()(0,0);
    double resVectorialGradv = integrate(_range=boundaryfaces(mesh),_expr= inner(gradv(uVectorial)*N(),gradv(uVectorial)*N()) ).evaluate()(0,0);
    double resVectorialDnv = integrate(_range=boundaryfaces(mesh),_expr=inner(dnv(uVectorial),dnv(uVectorial)) ).evaluate()(0,0);
    std::cout << "resVectorialExact : " << resVectorialExact << "\n";
    std::cout << "resVectorialGradv : " << resVectorialGradv << "\n";
    std::cout << "resVectorialDnv : " << resVectorialDnv << "\n";
    BOOST_CHECK_CLOSE( resVectorialGradv, resVectorialExact, 1e-10 );
    BOOST_CHECK_CLOSE( resVectorialDnv, resVectorialExact, 1e-10 );

}
BOOST_AUTO_TEST_SUITE_END()
