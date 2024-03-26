
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2015-04-25

  Copyright (C) 2015 UJF

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
   \file test_mortar.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2015-04-25
 */
#define BOOST_TEST_MODULE test_mortar2
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/moch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

namespace test_mortar2
{
using namespace Feel;

void run( int benchId, bool useVariantIntegrate = false )
{
    // case 1 : lm dirichlet conformal with submesh
    // case 2 : lm dirichlet conformal with no submesh relation
    // case 3 : lm dirichlet non conformal with no submesh relation

    typedef Mesh<Simplex<2>> mesh_type;
    double meshSize = doption(_name="gmsh.hsize");
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 4,1 );
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="Boundary1",_marker1=true);
    R.setMarker(_type="line",_name="Boundary2",_marker2=true,_marker3=true,_marker4=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new mesh_type,
                             _name="domainRectangle" );

    GeoTool::Node xl1( (benchId==2)? 0.0 : 1.0, (benchId==2)? 0.0 : 0.2 );
    GeoTool::Node xl2( (benchId==2)? 4.0 : 3.0, (benchId==2)? 0.0 : 0.7 );
    GeoTool::Line L( meshSize,"OMEGAL",xl1,xl2 );

    using trace_mesh_type = trace_mesh_t<mesh_type>;

    auto submesh = (benchId==1)?
        createSubmesh( _mesh=mesh, _range=markedfaces(mesh,"Boundary1") ) :
        L.createMesh(_mesh=new trace_mesh_type,_name="domainLine" );

    auto Vh = Pch<1>( mesh );
    auto u = Vh->elementPtr();
    std::cout << "Vh->nDof() " << Vh->nDof()<<"\n";

    auto Mh = Moch<1>( submesh );
    auto lambda = Mh->elementPtr();
    std::cout << "Mh->nDof() " << Mh->nDof()<<"\n";

    BlocksBaseGraphCSR myblockGraph(2,2);
    myblockGraph(0,0) = stencil(_test=Vh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    //myblockGraph(0,1) = stencil(_test=Vh,_trial=Mh, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(0,1) = stencil(_test=Mh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph()->transpose();
    myblockGraph(1,0) = stencil(_test=Mh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
    auto A = backend()->newBlockMatrix(_block=myblockGraph);

    BlocksBaseVector<double> myblockVec(2);
    myblockVec(0,0) = backend()->newVector( Vh );
    myblockVec(1,0) = backend()->newVector( Mh );
    auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

    BlocksBaseVector<double> myblockVecSol(2);
    myblockVecSol(0,0) = u;
    myblockVecSol(1,0) = lambda;
    auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);

    std::string fExpStr = "2*pi*pi*cos(pi*x)*sin(pi*y):x:y";
    std::string gExpStr = "cos(pi*x)*sin(pi*y):x:y";
    std::string hExpStr = gExpStr;
    if ( benchId == 3 )
    {
        fExpStr = "1";
        gExpStr = "x*(x-4):x";
        hExpStr = "0";
    }

    auto f = expr( fExpStr );
    auto g = expr( gExpStr );
    auto h = expr( hExpStr );

    form1( _test=Vh,_vector=F )
        += integrate(_range=elements(mesh),
                     _expr=f*id(u));

    form2( _trial=Vh, _test=Vh, _matrix=A)
        += integrate(_range=elements(mesh),
                     _expr=gradt(u)*trans(grad(u)) );

    if ( !useVariantIntegrate )
    {
        form1( _test=Mh, _vector=F,_rowstart=1)
            += integrate(_range=elements(submesh),
                         _expr=g*id(lambda) );
        form2( _trial=Mh, _test=Vh, _matrix=A,_colstart=1)
            += integrate( _range=elements(submesh),
                          _expr=idt(lambda)*id(u) );
        form2( _trial=Vh, _test=Mh, _matrix=A,_rowstart=1)
            += integrate(_range=elements(submesh),
                         _expr=id(lambda)*idt(u) );
    }
    else
    {
        form1( _test=Mh, _vector=F,_rowstart=1)
            += integrate(_range=markedfaces(mesh,"Boundary1"),
                         _expr=g*id(lambda) );
        form2( _trial=Mh, _test=Vh, _matrix=A,_colstart=1)
            += integrate( _range=markedfaces(mesh,"Boundary1"),
                          _expr=idt(lambda)*id(u) );
        form2( _trial=Vh, _test=Mh, _matrix=A,_rowstart=1)
            += integrate(_range=markedfaces(mesh,"Boundary1"),
                         _expr=id(lambda)*idt(u) );
    }

    std::list<std::string> listMarkOn = { "Boundary2" };
    if ( benchId == 3 )
        listMarkOn.push_back("Boundary1");
    form2( _trial=Vh, _test=Vh, _matrix=A)
        +=on(_range=markedfaces(mesh,listMarkOn), _rhs=F, _element=*u, _expr=h/*cst(0.)*/ );

    backend(_rebuild=true)->solve(_matrix=A,_rhs=F,_solution=U);

    myblockVecSol.localize(U);

    auto uExact = Vh->element(g);

    if ( benchId == 1 || benchId == 2 )
    {
        double thenormL2 = normL2(_range=elements(mesh),_expr=idv(u)-idv(uExact) );
        std::cout << "normL2 " << thenormL2 << "\n";
        BOOST_CHECK_SMALL( thenormL2, 1e-3 );
    }

    std::string tagName = (boost::format("BenchId%1%")%benchId).str();
    if (useVariantIntegrate) tagName += "Bis";
    auto e = exporter( _mesh=mesh,_name="export"+tagName );
    e->addRegions();
    e->add( "u"+tagName, *u );
    if ( benchId == 1 || benchId == 2 )
        e->add( "uExact"+tagName, uExact );
    e->save();


}

} // namespace test_mortar2

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( mortar2 )

BOOST_AUTO_TEST_CASE( mortar2_case1 )
{
    test_mortar2::run(1);
    test_mortar2::run(1,true);
}
BOOST_AUTO_TEST_CASE( mortar2_case2 )
{
    test_mortar2::run(2);
    test_mortar2::run(2,true);
}
BOOST_AUTO_TEST_CASE( mortar2_case3 )
{
    test_mortar2::run(3);
}

BOOST_AUTO_TEST_SUITE_END()
