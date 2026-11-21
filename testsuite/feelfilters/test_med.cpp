/* -*- mode: c++: coding: utf-8 -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-06-16

  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)

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
   \file test_med.cpp
   \author Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>
   \date 2016-12-01
 */

#define BOOST_TEST_MODULE med
#include <feel/feelcore/testsuite.hpp>
#include <boost/test/data/test_case.hpp>
#include <array>
#include <cmath>
#include <string>

#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;
namespace bdata = boost::unit_test::data;

namespace
{
using mesh_type = Mesh<Simplex<3>>;

const std::array<std::string,2> med_default_files = {{
    "data/geo/Cylref.med",
    "data/geo/tripod.med"
}};

void checkLoadMedMesh( std::string const& filename )
{
    std::string resolved = Environment::findFile( filename );
    BOOST_REQUIRE_MESSAGE( !resolved.empty(), "Unable to resolve MED file '" << filename << "'" );
    BOOST_REQUIRE_MESSAGE( fs::exists( resolved ), "MED file does not exist: '" << resolved << "'" );

    auto mesh = loadMesh( _mesh=new mesh_type,
                          _filename=resolved,
                          _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    BOOST_REQUIRE_MESSAGE( mesh, "MED mesh load returned null for '" << resolved << "'" );

    BOOST_CHECK_GT( mesh->numVertices(), 0 );
    BOOST_CHECK_GT( mesh->numFaces(), 0 );
    BOOST_CHECK_GT( mesh->numEdges(), 0 );
    BOOST_CHECK_GT( mesh->numElements(), 0 );

    double meshMeasure = mesh->measure();
    BOOST_CHECK( std::isfinite( meshMeasure ) );
    BOOST_CHECK_GT( meshMeasure, 0.0 );
}
}

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description medoptions("Med options");
    medoptions.add_options()
        ( "med.filename", Feel::po::value<std::string>()->default_value( "data/geo/Cylref.med" ), "name of the input MED file" )
        ;
    return medoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_med" ,
                           "test_med" ,
                           "0.1",
                           "test med integration with Feelpp",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2015-2016 Laboratoire national des Champs magnetiques Intenses");

    about.addAuthor("Christophe Trophime", "developer", "christophe.trophime@lncmi.cnrs.fr", "");
    return about;

}
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( inner_suite )

BOOST_DATA_TEST_CASE( test_default_med_files, bdata::make( med_default_files ), medfile )
{
    BOOST_TEST_CONTEXT( "medfile=" << medfile )
    {
        checkLoadMedMesh( medfile );
    }
}

BOOST_AUTO_TEST_CASE( test_configured_med_file )
{
    auto configuredMedFile = soption( _name="med.filename" );
    BOOST_TEST_CONTEXT( "med.filename=" << configuredMedFile )
    {
        checkLoadMedMesh( configuredMedFile );
    }
}

BOOST_AUTO_TEST_SUITE_END()
