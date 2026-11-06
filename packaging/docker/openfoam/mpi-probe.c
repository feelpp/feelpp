/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-License-Identifier: LGPL-3.0-or-later

/** @file mpi-probe.c @brief Check the OpenFOAM image's system MPI with two ranks. */
#include <mpi.h>
#include <stdio.h>

/** @brief Verify a collective and report the MPI library and rank placement. */
int main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );
    int rank = 0;
    int size = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    int contribution = rank + 1;
    int sum = 0;
    MPI_Allreduce( &contribution, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    char processor[MPI_MAX_PROCESSOR_NAME];
    int length = 0;
    MPI_Get_processor_name( processor, &length );
    printf( "rank=%d size=%d host=%.*s sum=%d\n", rank, size, length, processor, sum );
    if ( rank == 0 )
    {
        char version[MPI_MAX_LIBRARY_VERSION_STRING];
        MPI_Get_library_version( version, &length );
        printf( "%.*s\n", length, version );
    }
    int success = size == 2 && sum == 3;
    MPI_Finalize();
    return success ? 0 : 1;
}
