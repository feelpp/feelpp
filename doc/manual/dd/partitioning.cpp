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
        mesh = loadMesh(_mesh = new Mesh<Simplex<Dim>>(wComm), _worldcomm=wComm );
        //auto e = exporter( _mesh=mesh );
        //e->add("pid", regionProcess(Pdh<0>(mesh)) );
        //e->save();
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
                     _about=about(_name="doc_partitioning",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org") );
    Application app;
    app.add(new Partitioning<2>());
    //app.add(new Partitioning<3>());
    app.run();
}
