/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#define BOOST_TEST_MODULE test_vf_chi

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( test_vf_chi )

BOOST_AUTO_TEST_CASE( test_vf_chi_1 )
{
    using namespace Feel;

    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R( doption(_name="gmsh.hsize"),"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="Boundary",_markerAll=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new Mesh<Simplex<2>>,_name="domainRectangle" );


    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();
    auto myVec = vec( cst(2.),cst(-3.) );

    double evalX = integrate(_range=elements(mesh),
                             _expr=chi( (myVec(0,0)) > cst(0.) ) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE( evalX, 1., 1e-9 );
    double evalY = integrate(_range=elements(mesh),
                             _expr=chi( (myVec(1,0)) > cst(0.) ) ).evaluate()(0,0);
    BOOST_CHECK_CLOSE( evalY, 0., 1e-9 );

    form1(_test=Vh) =
        integrate(_range=elements(mesh),
                  _expr=chi( (myVec(0,0)) > cst(0.) )*id(u) );
    form2(_test=Vh,_trial=Vh) =
        integrate(_range=elements(mesh),
                  _expr=chi( (myVec(0,0)) > cst(0.) )*id(u)*idt(u) );

}
BOOST_AUTO_TEST_SUITE_END()
