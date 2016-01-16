/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Jan 2016
 
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
#define BOOST_TEST_MODULE test_updatemarker
#include <testsuite/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feeldiscr/pdh.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Update Marker  Options" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_updatemarker" ,
                     "test_updatemarker" ,
                     "0.2",
                     "nD(n=2,3) test updatemarker",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}



 FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
 BOOST_AUTO_TEST_SUITE( inner_suite )


BOOST_AUTO_TEST_CASE( test_elements )
{
    auto mesh2d = unitSquare();
    auto mesh = createSubmesh( mesh2d, elements(mesh2d), EXTRACTION_KEEP_MESH_RELATION );
    
    auto Xh = Pdh<0>(mesh);
    auto u = Xh->element(cst(1.));
    mesh->updateMarker2( u );
    BOOST_CHECK_EQUAL( nelements(marked2elements(mesh)), nelements(elements(mesh)) );
    u.on( _range=elements(mesh), _expr=cst(0.));
    mesh->updateMarker2( u );
    BOOST_CHECK_EQUAL( nelements(marked2elements(mesh)), 0 );
}


BOOST_AUTO_TEST_SUITE_END()

