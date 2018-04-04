/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
   \file dist2markedwalls.cpp
   \author Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
   \date 2014-01-21
 */

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feells/reinit_fms.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelmesh/iterator.hpp>

#define DIM 2

using namespace Feel;

void run()
{
    boost::timer chrono;

    typedef Mesh< Simplex<DIM> > mesh_type;
    auto mesh = loadMesh( _mesh=new mesh_type );

    auto Xh = Pch<1>(mesh);
    auto Xh0 = Pdh<0>(mesh);

    chrono.restart();
    // create the fast marching
    auto thefms = fms( Xh );
    double timeInitFastMarching = chrono.elapsed();

    auto phio = Xh->element();
    phio = vf::project(_space=Xh, _range=boundaryelements(mesh), _expr=h() );
    phio.on(_range=boundaryfaces(mesh), _expr=cst(0.) );

    // create a new range of elements having a marked face
    typedef boost::reference_wrapper<typename mesh_type::element_type const> element_ref_type;
    // store entities in a vector
    typedef std::vector<element_ref_type> cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    for (auto const& mface : markedfaces(mesh, {2,4}))
    {
        for (auto const& face : mface )
        {
            auto const& elt = face.element0();
            myelts->push_back(boost::cref(elt));
        }
    }

    auto myrange = boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                  myelts->begin(),myelts->end(),myelts );

    // mark the elements associated to myrange (they are considered as DONE by the fmm)
    auto mark = Xh0->element();
    mark.on( _range=myrange, _expr=cst(1) );

    mesh->updateMarker2( mark );

    chrono.restart();
    auto phi = thefms->march(phio, true);
    double timeFmsNoLocElt = chrono.elapsed();

    LOG(INFO) << "fast marching initialized in "<<timeInitFastMarching<<"s"<<std::endl;
    LOG(INFO) << "fast marching when elements DONE are given done in "<<timeFmsNoLocElt<<"s"<<std::endl;

    auto exp = exporter(_mesh=mesh, _name="disttowalls");
    exp->step(0)->add("phio", phio);
    exp->step(0)->add("phi", phi);
    exp->step(0)->add("mark", mark);
    exp->step(0)->addRegions();
    exp->save();

}

int main( int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv );
    run();
    return 0;
}
