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
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

using namespace Feel;

BOOST_AUTO_TEST_SUITE( test_meshmover )

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
