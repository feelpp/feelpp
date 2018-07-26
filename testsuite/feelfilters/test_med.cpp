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

#define USE_BOOST_TEST 1
//#undef USE_BOOST_TEST
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE med
#include <testsuite.hpp>
#endif

#include <string>

#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description medoptions("Med options");
    medoptions.add_options()
        ("medfile", Feel::po::value<std::string>()->default_value( "test_med.med" ), "name of the input MED file")
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

void runTest0()
{
    typedef Mesh<Simplex<3> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh;
    std::string filename = soption(_name="gmsh.filename");
    Feel::cout << "Loading " << filename << std::endl;

    if ( !filename.empty() )
    {
        mesh = loadMesh( _mesh=new mesh_type,
			 _filename=filename,
			 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );

	if ( Environment::isMasterRank() )
	  {
	    std::cout << "Number of vertices in mesh : " << mesh->numVertices() << std::endl;
	    std::cout << "Number of faces in mesh : " << mesh->numFaces() << std::endl;
	    std::cout << "Number of edges in mesh : " << mesh->numEdges() << std::endl;
	    std::cout << "Number of elts in mesh : " << mesh->numElements() << std::endl;
	  }
	
        // Feel::cout << "Markers:" << std::endl;
        // for (auto marker:  mesh->markerNames() )
        // {
        //     auto name = marker.first;
        //     auto data = marker.second;
        //     Feel::cout << "\t" << name << " dim=" << data[1] << " (" << name.size() << ")\n";
        // }
        // Feel::cout << std::endl;
    }
}


#if defined(USE_BOOST_TEST)

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( inner_suite )
BOOST_AUTO_TEST_CASE( test_0 )
{
    runTest0();
}
BOOST_AUTO_TEST_SUITE_END()
#else

int main( int argc, char* argv[] )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );
    runTest0();
}
#endif
