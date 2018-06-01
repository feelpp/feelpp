/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
   Date: 11 Dec 2017

   Copyright (C) 2017 Feel++ Consortium

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

#define BOOST_TEST_MODULE mesh mover testsuite

#include <testsuite.hpp>

#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

using namespace Feel;

BOOST_AUTO_TEST_SUITE( test_meshmover )

typedef boost::mpl::list<Mesh<Simplex<2,1>>,Mesh<Simplex<2,2>>,
                         Mesh<Simplex<3,1>>,Mesh<Simplex<3,2>>
                         //Mesh<Hypercube<2,1>>,Mesh<Hypercube<3,1>>
                         > test_mesh_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( test_allelements, MeshType, test_mesh_types )
{
    static const uint16_type nDim = MeshType::nDim;
    static const uint16_type nGeoOrder = MeshType::nOrder;
    auto mesh = loadMesh(_mesh=new MeshType);
    auto Vh = Pchv<nGeoOrder>( mesh );
    auto disp = Vh->element(one());
    std::map<size_type,std::vector<double>> M_baryByElt;
    for ( auto const& eltWrap : allelements(mesh) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        M_baryByElt[elt.id()].resize(nDim);
        auto thebary = elt.barycenter();
        for ( uint16_type d=0;d<nDim;++d )
            M_baryByElt[elt.id()][d] = thebary[d];
    }

    meshMove( mesh, disp );

    for ( auto const& eltWrap : allelements(mesh) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        auto thebary = elt.barycenter();
        for ( uint16_type d=0;d<nDim;++d )
        {
            BOOST_CHECK_CLOSE( M_baryByElt[elt.id()][d] + 1.0, thebary[d], 1e-9);
        }
    }
}

BOOST_AUTO_TEST_CASE( test_measures )
{
    GeoTool::Node x1( 0,0 );
    GeoTool::Node x2( 1,1 );
    GeoTool::Rectangle R( doption(_name="gmsh.hsize"),"OMEGA",x1,x2 );
    R.setMarker(_type="line",_name="Boundary",_markerAll=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh( _mesh=new Mesh<Simplex<2>>,_name="domainRectangle" );

    std::unordered_map<size_type,double> initialMeasures;
    for (auto const& eltWrap : elements(mesh) )
    {
        auto const& elt = unwrap_ref(eltWrap);
        initialMeasures[elt.id()] = elt.measure();
    }
    BOOST_CHECK_CLOSE(mesh->measure(),1.0,1e-11);

    // apply moving
    auto Vh = Pchv<1>( mesh );
    auto disp = Vh->element();
    disp.on(_range=elements(mesh),_expr=vec(Px(),Py()) );
    meshMove(mesh,disp);
    BOOST_CHECK_CLOSE(mesh->measure(),4.0,1e-11);

    // revert moving
    disp.scale(-1);
    meshMove(mesh,disp);
    BOOST_CHECK_CLOSE(mesh->measure(),1.0,1e-11);

    for (auto const& eltWrap : elements(mesh) )
    {
        auto const& elt = unwrap_ref(eltWrap);
        double measure1 = elt.measure();
        double measure2 = initialMeasures.find(elt.id())->second;
        BOOST_CHECK_CLOSE(measure1,measure2,1e-11);
    }

};

BOOST_AUTO_TEST_SUITE_END()
