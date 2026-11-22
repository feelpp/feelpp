// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <mpi.h>
#include <stdexcept>

/** Keep MPI alive during a test, including exception unwinding. */
class ScopedMpi
{
  public:
    /** Initialize the singleton test process without application arguments. */
    ScopedMpi()
    {
        if ( MPI_Init( nullptr, nullptr ) != MPI_SUCCESS )
            throw std::runtime_error( "MPI initialization failed" );
    }
    /** Finalize unless the tested lifecycle has already finalized MPI. */
    ~ScopedMpi()
    {
        int finalized = 0;
        MPI_Finalized( &finalized );
        if ( !finalized ) MPI_Finalize();
    }
    /** The test scope has a single finalization authority. */
    ScopedMpi( ScopedMpi const& ) = delete;

    /** MPI ownership cannot be reassigned. */
    ScopedMpi& operator=( ScopedMpi const& ) = delete;
};
