/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Chabannes Vincent <vincent.chabannes@feelpp.org>
   Date: 2014-09-04

   Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)
   Copyright (C) 2011-2016 Feel++ Consortium

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//#define USE_BOOST_TEST 1

#define BOOST_TEST_MODULE integration_relatedmesh testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>
//#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelvf/vf.hpp>
//#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>


using namespace Feel;
using namespace Feel::vf;

void run( bool useSMD )
{
    size_type ctxRelationLoc = (useSMD)? Feel::EXTRACTION_KEEP_MESH_RELATION : 0 ;

    // mesh
    double meshSize=doption(_name="gmsh.hsize");
    GeoTool::Rectangle R( meshSize,"Rectangle",
                          GeoTool::Node(-1,-1),
                          GeoTool::Node( 1, 1) );
    R.setMarker(_type="line",_name="Dirichlet",_markerAll=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);

    GeoTool::Rectangle Circ(meshSize,"OmegaCirc",
                            GeoTool::Node(-0.3,-0.3),
                            GeoTool::Node( 0.3,0.3));
    Circ.setMarker(_type="line",_name="cylinder",_markerAll=true);
    Circ.setMarker(_type="surface",_name="LocInterior",_markerAll=true);
    GeoTool::Rectangle Circ2(meshSize,"OmegaCirc2",
                             GeoTool::Node(-0.5,-0.5),
                             GeoTool::Node( 0.5,0.5) );
    Circ2.setMarker(_type="line",_name="FictitiousBoundary",_markerAll=true);
    Circ2.setMarker(_type="surface",_name="LocExterior",_markerAll=true);

    auto mesh1 = ((R-Circ2)+(Circ2-Circ)+Circ).createMesh(_mesh= new Mesh<Simplex<2,1,2>>,
                                                         _name="lomesh" );

    auto mesh2 = createSubmesh( _mesh=mesh1,_range=markedelements(mesh1,"LocInterior"),_context=ctxRelationLoc );

    std::list<std::string> listMarkers = { "LocInterior", "LocExterior" };
    auto mesh3 = createSubmesh( _mesh=mesh1,_range=markedelements(mesh1,listMarkers),_context=ctxRelationLoc );

    // function space
    auto Xh1 = Pch<2>(mesh1,true);
    auto Xh2 = Pch<3>(mesh2);
    auto Xh3 = Pch<2>(mesh3,true);
    auto u1 = Xh1->element(cst(1.));
    auto u2 = Xh2->element(cst(1.));
    auto u3 = Xh3->element(cst(1.));

    //-----------------------------------------------------------//
    //-----------------------------------------------------------//
    //-----------------------------------------------------------//
    // check bilinear form integration
    //-----------------------------------------------------------//
    //-----------------------------------------------------------//
    //-----------------------------------------------------------//

    auto rElementA = stencilRange<0,0>(markedelements(mesh1,"LocInterior"));
    auto graphElementA = ( ctxRelationLoc )?
        stencil(_test=Xh2,_trial=Xh1,_range=stencilRangeMap(rElementA) )->graph() :
        stencil(_test=Xh2,_trial=Xh1 )->graph();
    auto matElementA = backend()->newMatrix(0,0,0,0,graphElementA);
    auto bfElementA = form2(_test=Xh2,_trial=Xh1,_matrix=matElementA);
    bfElementA = integrate(_range=markedelements(mesh1,"LocInterior"),
                           _expr=idt(u1)*id(u2) );
    double energyElementA = bfElementA(u2,u1);
    BOOST_CHECK_CLOSE( energyElementA, 0.36, 1e-12 );
    Feel::cout << "energyElementA " << energyElementA << " [0.36]\n";

    //-----------------------------------------------------------//
    auto rElementB = stencilRange<0,0>(elements(mesh2));
    auto rElementBB = stencilRange<0,0>(markedelements(mesh1,"LocInterior"));
    auto graphElementB = ( ctxRelationLoc )?
        stencil(_test=Xh1,_trial=Xh2,_range=stencilRangeMap(rElementB) )->graph() :
        stencil(_test=Xh1,_trial=Xh2,_range=stencilRangeMap(rElementBB) )->graph();
    auto matElementB = backend()->newMatrix(0,0,0,0,graphElementB);
    auto bfElementB = form2(_test=Xh1,_trial=Xh2,_matrix=matElementB);
    bfElementB = integrate(_range=elements(mesh2),
                           _expr=idt(u2)*id(u1) );
    double energyElementB = bfElementB(u1,u2);
    BOOST_CHECK_CLOSE( energyElementB, 0.36, 1e-12 );
    Feel::cout << "energyElementB " << energyElementB << " [0.36]\n";

    //-----------------------------------------------------------//
    auto rElementC = stencilRange<0,0>(markedelements(mesh1,listMarkers));
    auto graphElementC = ( ctxRelationLoc )?
        stencil(_test=Xh3,_trial=Xh1,_range=stencilRangeMap(rElementC) )->graph():
        stencil(_test=Xh3,_trial=Xh1)->graph();
    auto matElementC = backend()->newMatrix(0,0,0,0,graphElementC);
    auto bfElementC = form2(_test=Xh3,_trial=Xh1,_matrix=matElementC);
    bfElementC = integrate(_range=markedelements(mesh1,listMarkers),
                           _expr=idt(u1)*id(u3) );
    double energyElementC = bfElementC(u3,u1);
    BOOST_CHECK_CLOSE( energyElementC, 1., 1e-12 );
    Feel::cout << "energyElementC " << energyElementC << " [1]\n";


    auto matElementD = backend()->newMatrix(_test=Xh3,_trial=Xh3);
    auto bfElementD = form2(_test=Xh3,_trial=Xh3,_matrix=matElementD);
    bfElementD = integrate(_range=elements(mesh3),
                           _expr=idt(u3)*id(u3) );
    double energyElementD = bfElementD(u3,u3);
    Feel::cout << "energyElementD " << energyElementD << " [1]\n";

    //-----------------------------------------------------------//
    //-----------------------------------------------------------//
    //-----------------------------------------------------------//

    auto rFacesStandard = stencilRange<0,0>(boundaryfaces(mesh2));
    auto graphFacesStandard = stencil(_test=Xh2,_trial=Xh2,
                                      _range=stencilRangeMap(rFacesStandard) )->graph();
    auto matFacesStandard = backend()->newMatrix(0,0,0,0,graphFacesStandard);
    auto bfFacesStandard = form2(_test=Xh2,_trial=Xh2,_matrix=matFacesStandard);
    bfFacesStandard = integrate(_range=boundaryfaces(mesh2),
                                _expr=idt(u2)*id(u2) );
    double energyFacesStandard = bfFacesStandard(u2,u2);
    BOOST_CHECK_CLOSE( energyFacesStandard, 2.4, 1e-12 );
    Feel::cout << "energyFacesStandard " << energyFacesStandard << " [2.4]\n";

    //-----------------------------------------------------------//
    auto rFacesNonStandardA = stencilRange<0,0>( boundaryfaces(mesh2) );
#if 0
    //auto rFacesNonStandardAA = stencilRange<0,0>( markedfaces(mesh1,"cylinder") ); // Bug with this!!
    auto rFacesNonStandardAA = stencilRange<0,0>(markedelements(mesh1,"LocInterior"));
    auto graphFacesNonStandardA = ( ctxRelationLoc )?
        stencil(_test=Xh1,_trial=Xh2,_range=stencilRangeMap(rFacesNonStandardA) )->graph() :
        stencil(_test=Xh1,_trial=Xh2,_range=stencilRangeMap(rFacesNonStandardAA) )->graph() ;
#else
    auto graphFacesNonStandardA = ( ctxRelationLoc )?
        stencil(_test=Xh1,_trial=Xh2,_range=stencilRangeMap(rFacesNonStandardA) )->graph() :
        stencil(_test=Xh2,_trial=Xh1,_range=stencilRangeMap(rFacesNonStandardA) )->graph()->transpose();
#endif
    auto matFacesNonStandardA = backend()->newMatrix(0,0,0,0,graphFacesNonStandardA);
    auto bfFacesNonStandardA = form2(_test=Xh1,_trial=Xh2,_matrix=matFacesNonStandardA);
    bfFacesNonStandardA = integrate(_range=boundaryfaces(mesh2),
                                    _expr=idt(u2)*id(u1) );
    double energyFacesNonStandardA = bfFacesNonStandardA(u1,u2);
    BOOST_CHECK_CLOSE( energyFacesNonStandardA, 2.4, 1e-12 );
    Feel::cout << "energyFacesNonStandardA " << energyFacesNonStandardA << " [2.4]\n";

    //-----------------------------------------------------------//

    auto rFacesNonStandardB = stencilRange<0,0>(markedfaces(mesh1,"cylinder"));// initial
    //auto rFacesNonStandardB = stencilRange<0,0>(boundaryfaces(mesh2)); // bad also
    //auto rFacesNonStandardB = stencilRange<0,0>(markedelements(mesh1,"LocInterior"));
    auto graphFacesNonStandardB = stencil(_test=Xh1,_trial=Xh2,_range=stencilRangeMap(rFacesNonStandardB) )->graph();
    auto matFacesNonStandardB = backend()->newMatrix(0,0,0,0,graphFacesNonStandardB);
    auto bfFacesNonStandardB = form2(_test=Xh1,_trial=Xh2,_matrix=matFacesNonStandardB);
    bfFacesNonStandardB = integrate(_range=markedfaces(mesh1,"cylinder"),//boundaryfaces(mesh2),
                                    _expr=idt(u2)*id(u1) );
    double energyFacesNonStandardB = bfFacesNonStandardB(u1,u2);
    BOOST_CHECK_CLOSE( energyFacesNonStandardB, 2.4, 1e-12 );
    Feel::cout << "energyFacesNonStandardB " << energyFacesNonStandardB << " [2.4]\n";

    //-----------------------------------------------------------//
    auto rFacesNonStandardC = stencilRange<0,0>(markedfaces(mesh1,"cylinder"));
    auto rFacesNonStandardCC = stencilRange<0,0>(boundaryfaces(mesh2));
    auto graphFacesNonStandardC = ( ctxRelationLoc )?
        stencil(_test=Xh2,_trial=Xh1,_range=stencilRangeMap(rFacesNonStandardC) )->graph() :
        stencil(_test=Xh2,_trial=Xh1,_range=stencilRangeMap(rFacesNonStandardCC) )->graph();
    auto matFacesNonStandardC = backend()->newMatrix(0,0,0,0,graphFacesNonStandardC);
    auto bfFacesNonStandardC = form2(_test=Xh2,_trial=Xh1,_matrix=matFacesNonStandardC );
    bfFacesNonStandardC = integrate(_range=markedfaces(mesh1,"cylinder"),
                                    _expr=idt(u1)*id(u2) );
    double energyFacesNonStandardC = bfFacesNonStandardC(u2,u1);
    BOOST_CHECK_CLOSE( energyFacesNonStandardC, 2.4, 1e-12 );
    BOOST_TEST_MESSAGE( fmt::format("energyFacesNonStandardC computed value [{}] should be equal to [2.4]",energyFacesNonStandardC) );

    //-----------------------------------------------------------//
    if ( ctxRelationLoc )
    {
        auto rFacesNonStandardD = stencilRange<0,0>(markedfaces(mesh3,"cylinder"));
        auto graphFacesNonStandardD = stencil(_test=Xh1,_trial=Xh3,
                                              _pattern=(size_type)Pattern::EXTENDED,
                                              _range=stencilRangeMap(rFacesNonStandardD) )->graph();
        auto matFacesNonStandardD = backend()->newMatrix(0,0,0,0,graphFacesNonStandardD);
        auto bfFacesNonStandardD = form2(_test=Xh1,_trial=Xh3,_matrix=matFacesNonStandardD);
        bfFacesNonStandardD = integrate(_range=markedfaces(mesh3,"cylinder"),//boundaryfaces(mesh2),
                                        _expr=idt(u3)*id(u1) );
        double energyFacesNonStandardD = bfFacesNonStandardD(u1,u3);
        double v = integrate(_range=markedfaces(mesh3,"cylinder"),//boundaryfaces(mesh2),
                             _expr=cst(1.) ).evaluate()(0,0);
        // faces are connected to 2 elements (left and right)
        double v1 = integrate(_range=markedfaces(mesh3,"cylinder"),//boundaryfaces(mesh2),
                              _expr=leftfacev(cst(1.))+rightfacev(cst(1.)) ).evaluate()(0,0);
        BOOST_CHECK_CLOSE( energyFacesNonStandardD, v1, 1e-12 ); // use to be 9.6
        Feel::cout << "energyFacesNonStandardD " << energyFacesNonStandardD << " [4.8] v=" << v << " v1=" << v1 << "\n";
    }

    //-----------------------------------------------------------//

    if ( ctxRelationLoc )
    {
        auto rFacesNonStandardE = stencilRange<0,0>(boundaryfaces(mesh2));
        auto rFacesNonStandardEE = stencilRange<0,0>(markedfaces(mesh1,"cylinder"));
        auto graphFacesNonStandardE = ( ctxRelationLoc )?
            stencil(_test=Xh1,_trial=Xh3,_range=stencilRangeMap(rFacesNonStandardE) )->graph() :
            stencil(_test=Xh1,_trial=Xh3,_range=stencilRangeMap(rFacesNonStandardEE) )->graph();
        auto matFacesNonStandardE = backend()->newMatrix(0,0,0,0,graphFacesNonStandardE);
        auto bfFacesNonStandardE = form2(_test=Xh1,_trial=Xh3,_matrix=matFacesNonStandardE);
        bfFacesNonStandardE = integrate(_range=boundaryfaces(mesh2),//markedfaces(mesh3,"cylinder"),//
                                        _expr=idt(u3)*id(u1) );
        double energyFacesNonStandardE = bfFacesNonStandardE(u1,u3);
        BOOST_CHECK_CLOSE( energyFacesNonStandardE, 2.4, 1e-12 );
        Feel::cout << "energyFacesNonStandardE " << energyFacesNonStandardE << " [2.4]\n";
    }

}

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( integration_relatedmesh )

BOOST_AUTO_TEST_CASE( integration_relatedmesh_useSMD )
{
    run(true);
}

BOOST_AUTO_TEST_CASE( integration_relatedmesh_notuseSMD )
{
    if ( Environment::worldComm().size() == 1 )
        run(false);
}


BOOST_AUTO_TEST_SUITE_END()
