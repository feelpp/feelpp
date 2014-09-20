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
#include <feel/feel.hpp>

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
    MPI_Comm comm;
    int p = 1;
    MPI_Group orig_group, new_group;
    MPI_Comm_group(Environment::worldComm(), &orig_group);
    int* pm = new int[p];
    for(unsigned short i = 0; i < p; ++i)
        pm[i] = i * (Environment::numberOfProcessors() / p);
    bool excluded = std::binary_search(pm, pm + p, Environment::rank());
    if(excluded)
        MPI_Group_incl(orig_group, p, pm, &new_group);
    else
        MPI_Group_excl(orig_group, p, pm, &new_group);
    MPI_Comm_create(Environment::worldComm(), new_group, &comm);
    MPI_Group_free(&orig_group);
    MPI_Group_free(&new_group);
    delete [] pm;
    boost::mpi::communicator bComm(comm, boost::mpi::comm_take_ownership);
    WorldComm wComm(bComm);
    boost::shared_ptr<Mesh<Simplex<Dim>>> mesh;
    if(!excluded)
        mesh = loadMesh(_mesh = new Mesh<Simplex<Dim>>, _worldcomm = wComm);
} // Partitioning::run

} // Feel

int main(int argc, char** argv) {
    using namespace Feel;
    /**
     * Initialize Feel++ Environment
     */
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="doc_partitioning",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org") );
    Application app;
    app.add(new Partitioning<2>());
    app.add(new Partitioning<3>());
    app.run();
}
