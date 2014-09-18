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



#if defined(USE_BOOST_MPI)
#include <boost/mpi.hpp>
#else
#include <mpi.h>
#endif
#if defined(USE_PETSC)
#include <petsc.h>
#endif
#include <stdio.h>
#include <iostream>

int main(int argc, char *argv[])
{
#if defined(USE_BOOST_MPI)
    boost::mpi::environment env(argc, argv,true);
#else
    MPI_Init(&argc, &argv);
#endif

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double* rhs = new double[size * 4];
    double* rhsSlave = new double[4];
    if(rank == 0)
    {
        // there is an issue here with MPI_IN_PLACE, normally we should not have to specify
        // the number of bytes and type for recv buffer but here we have to otherwise the scatter hangs
        // it seems that it occurs only when using boost::mpi
#if defined(USE_DATATYPE)
        MPI_Scatter(rhs, 4, MPI_DOUBLE, MPI_IN_PLACE, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD); // < this line doesn't hang
#else
        MPI_Scatter(rhs, 4, MPI_DOUBLE, MPI_IN_PLACE, 4, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD);
#endif

    }
    else
    {
        MPI_Scatter(NULL, 0, MPI_DATATYPE_NULL, rhsSlave, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
        std::cout << "Done with the scatter !" << std::endl;
    if(rank == 0)
    {
        std::fill(rhs, rhs + 4, rank);
        // there is an issue here with MPI_IN_PLACE, normally we should not have to specify
        // the number of bytes and type for send buffer but here we have to otherwise the gather returns wrong results
        // it seems that it occurs only when using boost::mpi
#if defined(USE_DATATYPE)
        MPI_Gather(MPI_IN_PLACE, 0, MPI_DOUBLE, rhs, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD); // < this line gives the correct results
#else
        MPI_Gather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rhs, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif


        for(int i = 0; i < size * 4; ++i) {
            std::cout << rhs[i] << " ";
            if((i + 1) % 4 == 0)
                std::cout << "(this line should be filled with " << i / 4 << ")" << std::endl;
        }
        std::cout << "Done with the gather !" << std::endl;
    }
    else
    {
        std::fill(rhsSlave, rhsSlave + 4, rank);
        MPI_Gather(rhsSlave, 4, MPI_DOUBLE, NULL, 0, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD);
    }
    delete [] rhs;
    delete [] rhsSlave;
    MPI_Finalize();
    return 0;
}
