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
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
using namespace Feel;

int main( int argc, char** argv )
{
    po::options_description opts ( "Mesh basic information and partition");
    opts.add_options()
        ( "numPartition", po::value<int>()->default_value(1), "Number of partitions" );

    // initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _desc=opts,
                     _about=about( _name="mesh" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    // create a mesh with GMSH using Feel++ geometry tool
    auto numPartition = ioption(_name="numPartition");
    auto mesh = loadMesh(_mesh=new  Mesh<CONVEX<FEELPP_DIM>>, _partitions=numPartition);

    size_type nbdyfaces = nelements(boundaryfaces(mesh));

    if ( Environment::isMasterRank() )
    {
        std::cout << " - mesh entities" << std::endl;
        std::cout << "      number of elements : " << mesh->numGlobalElements() << std::endl;
        std::cout << "         number of faces : " << mesh->numGlobalFaces() << std::endl;
        std::cout << "number of boundary faces : " << nbdyfaces << std::endl;
        if ( FEELPP_DIM > 2 )
            std::cout << "      number of edges : " << mesh->numGlobalEdges() << std::endl;  
        std::cout << "      number of points : " << mesh->numGlobalPoints() << std::endl;
        std::cout << "    number of vertices : " << mesh->numGlobalVertices() << std::endl;
        std::cout << " - mesh sizes" << std::endl;
        std::cout << "                h max : " << mesh->hMax() << std::endl;
        std::cout << "                h min : " << mesh->hMin() << std::endl;
        std::cout << "                h avg : " << mesh->hAverage() << std::endl;
        std::cout << "              measure : " << mesh->measure() << std::endl;

        std::cout << "Number of Partitions : " << numPartition << std::endl ;
    }

    for( auto marker: mesh->markerNames() )
    {
       auto name = marker.first;
       auto data = marker.second;
       if ( data[1] == mesh->dimension() )
       {

          auto meas1 = integrate( _range=elements( mesh ), _expr=cst(1.) ).evaluate();
          size_type nelts = nelements( markedelements(mesh, name ), true );
          auto meas = integrate( _range=markedelements( mesh, name ), _expr=cst(1.) ).evaluate();
          if ( Environment::isMasterRank() )
          {
            std::cout << " - Marker (elements) " << name << std::endl;
            std::cout << "    |- number of elements " << nelts << std::endl;
            std::cout << "    |- measure (elements(mesh) - markedelements(mesh,<<"name"<<) : " << meas << " -- " << meas1 << std::endl;

          }
       }
       else if ( data[1] == mesh->dimension()-1 )
       {
          size_type nelts = nelements( markedfaces(mesh, name ), true );
          auto meas = integrate( _range=markedfaces( mesh, name ), _expr=cst(1.) ).evaluate();
          if ( Environment::isMasterRank() )
          {
            std::cout << " - Marker (faces) " << name << std::endl;
            std::cout << "    |- number of faces " << nelts << std::endl;
            std::cout << "    |- measure : " << meas << std::endl;
          }
       }

    }
    // export results for post processing
    auto e = exporter( _mesh=mesh);
    e->addRegions();
    e->save();


}   // main
//! [all]
