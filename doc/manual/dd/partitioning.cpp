/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Pierre Jolivet <pierre.jolivet@imag.fr>
       Date: 2014-09-21

  Copyright (C) 2014 Université de Grenoble

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
   \file partitioning.cpp
   \author Pierre Jolivet <pierre.jolivet@imag.fr>
   \date 2014-09-21
 */
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feel.hpp>
#include <feel/options.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;

/**
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description partionningoptions( "Partitionning options" );
    partionningoptions.add_options()
    ( "use_global", po::value<bool>()->default_value( false ), "use global flag on nelements" )
      ;
    return partionningoptions.add( Feel::feel_options() );
}

namespace Feel
{
template<int Dim>
class Partitioning : public Simget
{
    typedef Simget super;
public:
    /**
     * Constructor
     */
    Partitioning()
        :
        super()
    {
    }

    void run();
}; // Partitioning

template<int Dim>
void
Partitioning<Dim>::run()
{
    int p = 2;
    int* pm = new int[p];
    for(unsigned short i = 0; i < p; ++i)
        pm[i] = i * (Environment::numberOfProcessors() / p);
    bool excluded = std::binary_search(pm, pm + p, Environment::rank());
    mpi::group new_group;
    if(excluded)
        new_group = Environment::worldComm().group().include(pm,pm+p);
    else
        new_group = Environment::worldComm().group().exclude(pm,pm+p);
    delete [] pm;

    boost::mpi::communicator bComm(Environment::worldComm(), new_group);
    std::vector<int> active( bComm.size(), true );
    WorldComm wComm(bComm,bComm,bComm,bComm.rank(),active);
    //wComm.showMe();
    boost::shared_ptr<Mesh<Simplex<Dim>>> mesh;
    if(!excluded)
    {
        std::cout << "proc " << Environment::rank() 
                  << " is not excluded and is locally rank " << wComm.rank() << " and loads mesh with " 
                  << wComm.globalSize() << " partitions\n";
        // mesh = loadMesh(_mesh = new Mesh<Simplex<Dim>>(wComm), _worldcomm=wComm );
        mesh = createGMSHMesh(_mesh = new Mesh<Simplex<Dim>>(wComm),
                              _worldcomm = wComm,
                              _desc = domain(_worldcomm = wComm, _name = "hypercube", _shape = "hypercube",
                                             _xmin = 0.0, _xmax = 1.0,
                                             _ymin = 0.0, _ymax = 1.0,
                                             _zmin = 0.0, _zmax = 1.0));
        std::cout << " - nelement(mesh)=" << nelements(elements(mesh)) << "\n";
        std::cout << " - loading space\n";
        auto Vh = Pch<2>( mesh );

        auto u = Vh->element("u");
        auto f = expr( soption(_name="functions.f"), "f", wComm );

        auto g = expr( soption(_name="functions.g"), "g", wComm );
        auto v = Vh->element( g, "g" );

        auto l = form1( _test=Vh );
        l = integrate(_range=elements(mesh),_expr=f*id(v));

        auto a = form2( _trial=Vh, _test=Vh);
        a = integrate(_range=elements(mesh),
                      _expr=gradt(u)*trans(grad(v)) );
        if(boption("gmsh.domain.usenames")) {
            if(nelements(markedfaces(mesh, "Dirichlet"), boption("use_global")) > 0)
                a+=on(_range=markedfaces(mesh, "Dirichlet"), _rhs=l, _element=u, _expr=g);
        }
        else
            a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g);
        a.solve(_rhs=l,_solution=u);
        auto e = exporter( _mesh=mesh );
        //e->addRegions();
        e->add( "u", u );
        e->add( "g", v );
        e->save();
    }
    else
    {
        std::cout << "proc " << Environment::rank() 
                  << " is excluded and does not load mesh";
        int np = 0;
        mpi::all_reduce(wComm, 1, np, std::plus<int>());
        std::cout << "proc " << Environment::rank() 
                  <<   " - nb proc = " << np << "\n";
    }
} // Partitioning::run

} // Feel

int main(int argc, char** argv) {
    using namespace Feel;
    /**
     * Initialize Feel++ Environment
     */
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="partitioning",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org") );
    Application app;
    app.add(new Partitioning<2>());
    //app.add(new Partitioning<3>());
    app.run();
}
