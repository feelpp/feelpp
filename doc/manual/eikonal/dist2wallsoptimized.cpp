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
   \file dist2wallsoptimized.cpp
   \author Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
   \date 2014-01-21
 */

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelpde/reinit_fms.hpp>
#include <feel/feelvf/vf.hpp>

#define DIM 2

using namespace Feel;

void run()
{
    boost::timer chrono;

    typedef Mesh< Simplex<DIM> > mesh_type;
    auto mesh = loadMesh( _mesh=new mesh_type );

    auto Xh = Pch<1>(mesh);

    chrono.restart();
    // create the fast marching
    auto thefms = fms( Xh );
    double timeInitFastMarching = chrono.elapsed();

    // first method: let the fm search for the elements crossed by the interface
    auto phio = Xh->element();
    phio = vf::project(Xh, elements(mesh), h() );
    phio +=vf::project(Xh, boundaryfaces(mesh), -idv(phio) - h()/100. );

    chrono.restart();
    auto phi = thefms->march(phio);
    double timeFmsLocElt = chrono.elapsed();


    // second method: give to the fm the elements having the good value to start with
    auto Xh0 = Pdh<0>(mesh);

    auto phio2 = Xh->element();    
    phio2 = vf::project(Xh, boundaryelements(mesh), h() );
    phio2 += vf::project(Xh, boundaryfaces(mesh), -idv(phio2) );

    auto mark = vf::project(Xh0, boundaryelements(mesh), cst(1) );
    mesh->updateMarker2( mark );

    chrono.restart();
    auto phi2 = thefms->march(phio2, true);
    double timeFmsNoLocElt = chrono.elapsed();

    LOG(INFO) << "fast marching initialized in "<<timeInitFastMarching<<"s"<<std::endl;
    LOG(INFO) << "fast marching locating the elements done in "<<timeFmsLocElt<<"s"<<std::endl;
    LOG(INFO) << "fast marching when elements DONE are given done in "<<timeFmsNoLocElt<<"s"<<std::endl;


    auto exp = exporter(_mesh=mesh, _name="disttowalls");
    exp->step(0)->add("phio", phio);
    exp->step(0)->add("phio2", phio2);
    exp->step(0)->add("phi", phi);
    exp->step(0)->add("phi2", phi2);
    exp->step(0)->add("mark", mark);
    exp->save();

}

int main( int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv );
    run();
    return 0;
}
