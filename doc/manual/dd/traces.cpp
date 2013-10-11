/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-10-11

  Copyright (C) 2013 Université de Strasbourg

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
   \file traces.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-10-11
 */
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="doc_traces",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    //auto Dh = Pdh<1>( mesh );
    //auto u = Dh->element();
    //auto e = exporter( _mesh=mesh );
    BOOST_FOREACH( auto neighbor_subdomain, mesh->faceNeighborSubdomains() )
    {
        LOG(INFO) << "Extracting trace mesh from neighbor : " << neighbor_subdomain;
        auto trace = createSubmesh( mesh, interprocessfaces(mesh, neighbor_subdomain ) );
        LOG(INFO) << "number of elements in trace mesh " << nelements(elements(trace)) <<      " (" << env.rank() << " vs. " << neighbor_subdomain << ")";
        LOG(INFO) << "number of faces in trace mesh    " << nelements(boundaryfaces(trace)) << " (" << env.rank() << " vs. " << neighbor_subdomain << ")";
        auto Xh = Pch<1>( trace );
        auto l = Xh->element();
        l = vf::project( Xh, elements(trace), cst( neighbor_subdomain ) );
        auto m = mean(elements(trace), idv(l))(0,0);
        CHECK( math::abs( m -  neighbor_subdomain ) < 1e-14 ) << "problem : " << m << " != " << neighbor_subdomain;
        // need to introduce interpolation operator with current subdomaine to
        // interpolate l onto the subdomain and set on the interfaces l
        // accumulate on all interfaces in u and visualize u
    }
    //e->add( "u", u );
    //e->save();
}
