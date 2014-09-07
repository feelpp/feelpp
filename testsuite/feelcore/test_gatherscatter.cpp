/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-03-10

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
   \file test_gatherscatter.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-03-10
 */

#define FEELPP_BUG
#ifdef FEELPP_BUG
#include <feel/feel.hpp>

namespace Feel
{
template<int Dim, int Order>
class Debugpp : public Simget {
    typedef Simget super;
public:
    Debugpp() : super() { }
    void run();
};

template<int Dim,int Order>
void Debugpp<Dim,Order>::run() {
#else
#include <mpi.h>
#include <stdio.h>
#include <iostream>
    int main(int argc, char *argv[]) {
        MPI_Init(&argc, &argv);
#endif
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        double* rhs = new double[size * 4];
        double* rhsSlave = new double[4];
        if(rank == 0)
        {
            //MPI_Scatter(rhs, 4, MPI_DOUBLE, rhsSlave, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            //MPI_Scatter(rhs, 4, MPI_DOUBLE, rhsSlave, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            //MPI_Scatter(rhs, 4, MPI_DOUBLE, rhsSlave, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            //mpi::scatter( MPI_COMM_WORLD, rhs, rhsSlave, 4, 0 );
            // there is an issue here with MPI_IN_PLACE, normally we should not have to specify
            // the number of bytes and type for recv buffer but here we have to otherwise the scatter hangs
            // it seems that it occurs only when using boost::mpi
            MPI_Scatter(rhs, 4, MPI_DOUBLE, MPI_IN_PLACE, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Scatter(NULL, 0, MPI_DATATYPE_NULL, rhsSlave, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            //mpi::scatter( MPI_COMM_WORLD, rhsSlave, 4, 0 );
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0)
            std::cout << "Done !" << std::endl;
        delete [] rhs;
        delete [] rhsSlave;
#ifdef FEELPP_BUG
    }
}

int main(int argc, char** argv) {
    using namespace Feel;
    Environment env(_argc=argc, _argv=argv,
                    _desc=feel_options(),
                    _about=about(_name="feelpp_feti",
                                 _author="Feel++ Consortium",
                                 _email="feelpp-devel@feelpp.org"));
    Application app;
    app.add(new Debugpp<2,1>());
    app.run();
}
#else
MPI_Finalize();
return 0;
}
#endif
