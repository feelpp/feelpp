/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-06-11

  Copyright (C) 2014 Feel++ Consortium

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
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
using namespace Feel;

int main( int argc, char** argv )
{
    // initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="mesh" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    // create a mesh with GMSH using Feel++ geometry tool
    auto mesh = loadMesh(_mesh=new  Mesh<Simplex<FEELPP_DIM>>);


    if ( Environment::isMasterRank() )
    {
        std::cout << " - mesh entities" << std::endl;
        std::cout << "   number of elements : " << mesh->numGlobalElements() << std::endl;
        std::cout << "      number of faces : " << mesh->numGlobalFaces() << std::endl;
        if ( FEELPP_DIM > 2 )
            std::cout << "      number of edges : " << mesh->numGlobalEdges() << std::endl;
        std::cout << "      number of points : " << mesh->numGlobalPoints() << std::endl;
        std::cout << "    number of vertices : " << mesh->numGlobalVertices() << std::endl;
        std::cout << " - mesh sizes" << std::endl;
        std::cout << "                h max : " << mesh->hMax() << std::endl;
        std::cout << "                h min : " << mesh->hMin() << std::endl;
        std::cout << "                h avg : " << mesh->hAverage() << std::endl;
        std::cout << "              measure : " << mesh->measure() << std::endl;
    }

    for( auto marker: mesh->markerNames() )
    {
       auto name = marker.first;
       auto data = marker.second;
       if ( data[1] == mesh->dimension() )
       {
          size_type nelts = nelements( markedelements(mesh, name ), true );
          if ( Environment::isMasterRank() )
          {
            std::cout << "Marker (elements) " << name << std::endl;
            std::cout << " - number of elements " << nelts << std::endl;

          }
       }
       else if ( data[1] == mesh->dimension()-1 )
       {
          size_type nelts = nelements( markedfaces(mesh, name ), true );
          std::cout << "Marker (faces) " << name << std::endl;
            std::cout << " - number of faces " << nelts << std::endl;
       }

    }
    // export results for post processing
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->save();


}   // main
//! [all]
