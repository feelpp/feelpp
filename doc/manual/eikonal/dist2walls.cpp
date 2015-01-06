/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file dist2walls.cpp
   \author Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
   \date 2014-01-21
 */


#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelpde/reinit_fms.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

void run()
{
    typedef Mesh< FEELPP_CONVEX<FEELPP_DIM> > mesh_type;

    auto mesh = loadMesh( _mesh=new mesh_type );

    auto Xh = Pch<1>(mesh);

    auto thefms = fms( Xh );

    auto phio = Xh->element();
    phio = vf::project(Xh, elements(mesh), h() );
    if ( !soption("marker").empty() )
        phio +=vf::project(Xh, markedfaces(mesh,soption("marker")), -idv(phio) - h()/100. );
    else
        phio +=vf::project(Xh, boundaryfaces(mesh), -idv(phio) - h()/100. );
    auto phi = thefms->march(phio);

    auto exp = exporter(_mesh=mesh, _name="disttowalls");
    exp->step(0)->add("phi", phi);
    exp->save();

}

int main( int argc, char** argv )
{
    po::options_description opts( "Dist2Walls options" );
    opts.add_options()
    ( "marker", po::value<std::string>()->default_value( "" ),
      "marker of the boundary from which to compute the distance function, if empty use the entire boundary" )
    ;
    Feel::Environment env( _argc=argc, _argv=argv, _desc=opts );
    run();
    return 0;
}
