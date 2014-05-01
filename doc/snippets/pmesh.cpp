/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-04-30

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
/**
   \file pmesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-04-30
 */
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

int main(int argc, char** argv )
{
    using namespace Feel;

    // initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="pmesh" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    if ( env.isMasterRank() )
    {
        // get number of partitions
        int nparts = ioption( "gmsh.npartitions" );

        //! [load]
        std::cout << "Loading mesh and partition it into " << nparts << " partitions\n";
        // create a mesh with GMSH using Feel++ geometry tool
        auto mesh = loadMesh(_mesh=new  Mesh<Hypercube<2>>,
                             _partitions=nparts,
                             _respect_partition=true,
                             _worldcomm=Environment::worldCommSeq() );
        //! [load]
        std::cout << "looking for first element(masters) in each partition\n";
        std::vector<int> masters( nparts, -1 );
        for( auto elt : allelements( mesh ) )
        {
            int part  = elt.processId();
            int eid  = elt.id();
            //std::cout << "element " << eid << " in partition " << part << "\n";
            if ( masters[part] == -1 )
            {
                std::cout << " - in partition " << part << ", master is " << eid << "\n";
                masters[part] = eid;
            }
            bool ok = true;
            for( int p : masters )
            {
                ok = ok && (p != -1);
            }
            if ( ok )
                break;
        }

        std::cout << "Masters(first elements in each partition): { ";
        for( int m : masters ) std::cout << m << " ";
        std::cout << "}";

        //! [export]
        // export results for post processing
        auto e = exporter( _mesh=mesh );
        e->addRegions();
        e->save();

        //! [export]
    }
}   // main
