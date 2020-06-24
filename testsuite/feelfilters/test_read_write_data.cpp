/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 24 june 2020

 Copyright (C) 2020 Feel++ Consortium

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

// give a name to the testsuite
#define BOOST_TEST_MODULE read write data testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pch.hpp>
//#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelvf/vf.hpp>
// #include <feel/feelfilters/exporter.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( read_write_data )

BOOST_AUTO_TEST_CASE( test_mesh )
{
    using namespace Feel;

    using mesh_type = Mesh<Simplex<2>>;
    //auto mesh = loadMesh(_mesh=new mesh_type);
    double meshSize = doption(_name="gmsh.hsize");
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="Boundary1",_marker1=true);
    R.setMarker(_type="line",_name="Boundary2",_marker2=true);
    R.setMarker(_type="line",_name="Boundary3",_marker3=true);
    R.setMarker(_type="line",_name="Boundary4",_marker4=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new mesh_type,
                             _name="domainRectangle" );
    std::vector<std::string> namesOfMarkedFaces = { "Boundary1","Boundary2","Boundary3","Boundary4" };


    size_type nMeshElt =  mesh->numGlobalElements();
    double areaElements = measure(_range=elements(mesh));
    double areaBoundaryFaces = measure(_range=boundaryfaces(mesh));
    std::vector<double> areaMarkedFaces(namesOfMarkedFaces.size());
    for ( int k=0;k<namesOfMarkedFaces.size();++k )
        areaMarkedFaces[k] = measure(_range=markedfaces(mesh,namesOfMarkedFaces[k]));

    //std::cout << "nGhost : " << nelements( elements(mesh,EntityProcessType::GHOST_ONLY) ) << std::endl;

    std::string filenameOfMeshSaved = (boost::format("mymesh_%1%.json")%mesh_type::shape_type::name() ).str();
    std::string fullpathOfMeshSaved = (fs::current_path()/filenameOfMeshSaved).string();
    mesh->saveHDF5( fullpathOfMeshSaved );

    auto meshSeqReloaded = std::make_shared<mesh_type>( Environment::worldCommSeqPtr() );
    meshSeqReloaded->loadHDF5( fullpathOfMeshSaved );

    BOOST_CHECK( nMeshElt == meshSeqReloaded->numGlobalElements() );

    double areaElementsReloaded = integrate(_range=elements(meshSeqReloaded),_expr=cst(1.)).evaluate( true, meshSeqReloaded->worldCommPtr() )(0,0);
    double areaBoundaryFacesReloaded = integrate(_range=boundaryfaces(meshSeqReloaded),_expr=cst(1.)).evaluate( true, meshSeqReloaded->worldCommPtr() )(0,0);
    BOOST_CHECK_CLOSE( areaElements, areaElementsReloaded, 1e-10 );
    BOOST_CHECK_CLOSE( areaBoundaryFaces, areaBoundaryFacesReloaded, 1e-10 );
    for ( int k=0;k<namesOfMarkedFaces.size();++k )
    {
        double areaMarkedFacesReloaded = integrate(_range=markedfaces(meshSeqReloaded,namesOfMarkedFaces[k]),_expr=cst(1.)).evaluate( true, meshSeqReloaded->worldCommPtr() )(0,0);
        BOOST_CHECK_CLOSE( areaMarkedFaces[k], areaMarkedFacesReloaded, 1e-10 );
    }

#if 0
    if ( Environment::isMaterRank() )
    {
        auto myexporter1 = exporter( _mesh=meshSeqReloaded, _name="meshSeqReloaded" );
        myexporter1->save();
    }

     auto myexporter2 = exporter( _mesh=mesh, _name="mesh" );
     myexporter2->save();
#endif
}


BOOST_AUTO_TEST_SUITE_END()
