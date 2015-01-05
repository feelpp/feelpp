/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Chabannes Vincent <vincent.chabannes@feelpp.org>
   Date: 2014-09-04

   Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)
   Copyright (C) 2011-2015 Feel++ Consortium

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

#include <testsuite/testsuite.hpp>

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

    auto mesh2 = createSubmesh( mesh1,markedelements(mesh1,"LocInterior"),ctxRelationLoc );

    std::list<std::string> listMarkers = { "LocInterior", "LocExterior" };
    auto mesh3 = createSubmesh( mesh1,markedelements(mesh1,listMarkers),ctxRelationLoc );

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
    if ( Environment::isMasterRank() )
        std::cout << "energyElementA " << energyElementA << " [0.36]\n";

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
    if ( Environment::isMasterRank() )
        std::cout << "energyElementB " << energyElementB << " [0.36]\n";

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
    if ( Environment::isMasterRank() )
        std::cout << "energyElementC " << energyElementC << " [1]\n";


    auto matElementD = backend()->newMatrix(_test=Xh3,_trial=Xh3);
    auto bfElementD = form2(_test=Xh3,_trial=Xh3,_matrix=matElementD);
    bfElementD = integrate(_range=elements(mesh3),
                           _expr=idt(u3)*id(u3) );
    double energyElementD = bfElementD(u3,u3);
    if ( Environment::isMasterRank() )
        std::cout << "energyElementD " << energyElementD << " [1]\n";

    //-----------------------------------------------------------//
    //-----------------------------------------------------------//
    //-----------------------------------------------------------//

    auto rFacesStandart = stencilRange<0,0>(boundaryfaces(mesh2));
    auto graphFacesStandart = stencil(_test=Xh2,_trial=Xh2,
                                      _range=stencilRangeMap(rFacesStandart) )->graph();
    auto matFacesStandart = backend()->newMatrix(0,0,0,0,graphFacesStandart);
    auto bfFacesStandart = form2(_test=Xh2,_trial=Xh2,_matrix=matFacesStandart);
    bfFacesStandart = integrate(_range=boundaryfaces(mesh2),
                                _expr=idt(u2)*id(u2) );
    double energyFacesStandart = bfFacesStandart(u2,u2);
    BOOST_CHECK_CLOSE( energyFacesStandart, 2.4, 1e-12 );
    if ( Environment::isMasterRank() )
        std::cout << "energyFacesStandart " << energyFacesStandart << " [2.4]\n";

    //-----------------------------------------------------------//
    auto rFacesNonStandartA = stencilRange<0,0>( boundaryfaces(mesh2) );
#if 0
    //auto rFacesNonStandartAA = stencilRange<0,0>( markedfaces(mesh1,"cylinder") ); // Bug with this!!
    auto rFacesNonStandartAA = stencilRange<0,0>(markedelements(mesh1,"LocInterior"));
    auto graphFacesNonStandartA = ( ctxRelationLoc )?
        stencil(_test=Xh1,_trial=Xh2,_range=stencilRangeMap(rFacesNonStandartA) )->graph() :
        stencil(_test=Xh1,_trial=Xh2,_range=stencilRangeMap(rFacesNonStandartAA) )->graph() ;
#else
    auto graphFacesNonStandartA = ( ctxRelationLoc )?
        stencil(_test=Xh1,_trial=Xh2,_range=stencilRangeMap(rFacesNonStandartA) )->graph() :
        stencil(_test=Xh2,_trial=Xh1,_range=stencilRangeMap(rFacesNonStandartA) )->graph()->transpose();
#endif
    auto matFacesNonStandartA = backend()->newMatrix(0,0,0,0,graphFacesNonStandartA);
    auto bfFacesNonStandartA = form2(_test=Xh1,_trial=Xh2,_matrix=matFacesNonStandartA);
    bfFacesNonStandartA = integrate(_range=boundaryfaces(mesh2),
                                    _expr=idt(u2)*id(u1) );
    double energyFacesNonStandartA = bfFacesNonStandartA(u1,u2);
    BOOST_CHECK_CLOSE( energyFacesNonStandartA, 2.4, 1e-12 );
    if ( Environment::isMasterRank() )
        std::cout << "energyFacesNonStandartA " << energyFacesNonStandartA << " [2.4]\n";

    //-----------------------------------------------------------//

    auto rFacesNonStandartB = stencilRange<0,0>(markedfaces(mesh1,"cylinder"));// initial
    //auto rFacesNonStandartB = stencilRange<0,0>(boundaryfaces(mesh2)); // bad also
    //auto rFacesNonStandartB = stencilRange<0,0>(markedelements(mesh1,"LocInterior"));
    auto graphFacesNonStandartB = stencil(_test=Xh1,_trial=Xh2,_range=stencilRangeMap(rFacesNonStandartB) )->graph();
    auto matFacesNonStandartB = backend()->newMatrix(0,0,0,0,graphFacesNonStandartB);
    auto bfFacesNonStandartB = form2(_test=Xh1,_trial=Xh2,_matrix=matFacesNonStandartB);
    bfFacesNonStandartB = integrate(_range=markedfaces(mesh1,"cylinder"),//boundaryfaces(mesh2),
                                    _expr=idt(u2)*id(u1) );
    double energyFacesNonStandartB = bfFacesNonStandartB(u1,u2);
    BOOST_CHECK_CLOSE( energyFacesNonStandartB, 2.4, 1e-12 );
    if ( Environment::isMasterRank() )
        std::cout << "energyFacesNonStandartB " << energyFacesNonStandartB << " [2.4]\n";

    //-----------------------------------------------------------//
    auto rFacesNonStandartC = stencilRange<0,0>(markedfaces(mesh1,"cylinder"));
    auto rFacesNonStandartCC = stencilRange<0,0>(boundaryfaces(mesh2));
    auto graphFacesNonStandartC = ( ctxRelationLoc )?
        stencil(_test=Xh2,_trial=Xh1,_range=stencilRangeMap(rFacesNonStandartC) )->graph() :
        stencil(_test=Xh2,_trial=Xh1,_range=stencilRangeMap(rFacesNonStandartCC) )->graph();
    auto matFacesNonStandartC = backend()->newMatrix(0,0,0,0,graphFacesNonStandartC);
    auto bfFacesNonStandartC = form2(_test=Xh2,_trial=Xh1,_matrix=matFacesNonStandartC );
    bfFacesNonStandartC = integrate(_range=markedfaces(mesh1,"cylinder"),
                                    _expr=idt(u1)*id(u2) );
    double energyFacesNonStandartC = bfFacesNonStandartC(u2,u1);
    BOOST_CHECK_CLOSE( energyFacesNonStandartC, 2.4, 1e-12 );
    if ( Environment::isMasterRank() )
        std::cout << "energyFacesNonStandartC " << energyFacesNonStandartC << " [2.4]\n";

    //-----------------------------------------------------------//
    if ( ctxRelationLoc )
    {
        auto rFacesNonStandartD = stencilRange<0,0>(markedfaces(mesh3,"cylinder"));
        auto graphFacesNonStandartD = stencil(_test=Xh1,_trial=Xh3,
                                              _pattern=(size_type)Pattern::EXTENDED,
                                              _range=stencilRangeMap(rFacesNonStandartD) )->graph();
        auto matFacesNonStandartD = backend()->newMatrix(0,0,0,0,graphFacesNonStandartD);
        auto bfFacesNonStandartD = form2(_test=Xh1,_trial=Xh3,_matrix=matFacesNonStandartD);
        bfFacesNonStandartD = integrate(_range=markedfaces(mesh3,"cylinder"),//boundaryfaces(mesh2),
                                        _expr=idt(u3)*id(u1) );
        double energyFacesNonStandartD = bfFacesNonStandartD(u1,u3);
        BOOST_CHECK_CLOSE( energyFacesNonStandartD, 9.6, 1e-12 );
        if ( Environment::isMasterRank() )
            std::cout << "energyFacesNonStandartD " << energyFacesNonStandartD << " [9.6]\n";
    }

    //-----------------------------------------------------------//

    if ( ctxRelationLoc )
    {
        auto rFacesNonStandartE = stencilRange<0,0>(boundaryfaces(mesh2));
        auto rFacesNonStandartEE = stencilRange<0,0>(markedfaces(mesh1,"cylinder"));
        auto graphFacesNonStandartE = ( ctxRelationLoc )?
            stencil(_test=Xh1,_trial=Xh3,_range=stencilRangeMap(rFacesNonStandartE) )->graph() :
            stencil(_test=Xh1,_trial=Xh3,_range=stencilRangeMap(rFacesNonStandartEE) )->graph();
        auto matFacesNonStandartE = backend()->newMatrix(0,0,0,0,graphFacesNonStandartE);
        auto bfFacesNonStandartE = form2(_test=Xh1,_trial=Xh3,_matrix=matFacesNonStandartE);
        bfFacesNonStandartE = integrate(_range=boundaryfaces(mesh2),//markedfaces(mesh3,"cylinder"),//
                                        _expr=idt(u3)*id(u1) );
        double energyFacesNonStandartE = bfFacesNonStandartE(u1,u3);
        BOOST_CHECK_CLOSE( energyFacesNonStandartE, 2.4, 1e-12 );
        if ( Environment::isMasterRank() )
            std::cout << "energyFacesNonStandartE " << energyFacesNonStandartE << " [2.4]\n";
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

