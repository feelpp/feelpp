/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-06-11

  Copyright (C) 2014-2016 Feel++ Consortium

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
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/partitionmesh.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelvf/vf.hpp>
using namespace Feel;

int main( int argc, char** argv )
{
    // initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="mesh" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    tic();
    //auto mesh = loadMesh(_mesh=new  Mesh<CONVEX<FEELPP_DIM>>, _partitions=1, _savehdf5=0 );
    Feel::cout << "mesh.save.enable=" << boption("mesh.save.enable") << std::endl << std::flush;
    auto mesh = loadMesh(_mesh=new  Mesh<CONVEX<FEELPP_DIM>>,_savehdf5=boption("mesh.save.enable"), _filename=soption("mesh.filename"),
                         _update=size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES));
    toc("loading mesh done",FLAGS_v>0);

    if ( boption("mesh.partition.enable") && Environment::numberOfProcessors() == 1 )
    {
        // build a MeshPartitionSet based on a mesh partition that will feed a
        // partition io data structure to generate a parallel hdf5 file from which
        // the parallel mesh can be loaded
        using io_t = PartitionIO<mesh_t<decltype(mesh)>>;
        io_t io( fs::path(soption("mesh.filename")).stem().string()+".json" );
        io.write( partitionMesh( mesh, ioption("mesh.partition.size") ) );
        return 0;
    }

    auto Xhd0 = Pdh<0>(mesh);
    auto measures = Xhd0->element();
    measures.on(_range=elements(mesh),_expr=vf::meas());
    double measMin = measures.min();
    double measMax = measures.max();
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
        std::cout << "              measure : " << mesh->measure() << "\t" << measMin<<" : " << measMax << std::endl;

        std::cout << "Number of Partitions : " << mesh->numberOfPartitions() << std::endl ;
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
            std::cout << "    |- measure (elements(mesh) - markedelements(mesh,"
                      <<name <<") : " << meas << " -- " << meas1 << std::endl;

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
