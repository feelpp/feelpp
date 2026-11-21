// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-3.0-or-later
#define BOOST_TEST_MODULE kokkos_runtime
#include "runtimefixture.hpp"
#include <boost/test/unit_test.hpp>
#include <feel/feelconfig.h>
#include <feel/feelcore/kokkosruntime.hpp>
#include <mpi.h>
#include <petscsys.h>
#ifdef FEELPP_HAS_KOKKOS
#include <Kokkos_Core.hpp>
#endif
#include <iostream>
#include <stdexcept>
#include <string>

using namespace Feel;
/** Run each terminal runtime lifecycle in an independent process. */
void runRuntimeScenario( std::string const& mode )
{
    int argc = 0;
    char** argv = nullptr;
    try
    {
        if ( mode == "before-mpi" )
        {
            KokkosRuntime runtime;
            BOOST_CHECK_THROW( runtime.initialize(), std::logic_error );
            return;
        }
        ScopedMpi mpi;
        if ( mode == "after-mpi" )
        {
            MPI_Finalize();
            KokkosRuntime runtime;
            BOOST_CHECK_THROW( runtime.initialize(), std::logic_error );
            return;
        }
        auto caps = KokkosRuntime::capabilities();
        std::cout << "requested=" << mode << " available=" << caps.available
                  << " host=" << caps.hostExecutionSpace << " default=" << caps.defaultExecutionSpace
                  << " petsc-kokkos=" << caps.petsc << '\n';
        if ( mode == "off" || mode == "lazy" )
        {
            KokkosRuntime runtime( mode == "off" ? KokkosInitPolicy::off : KokkosInitPolicy::lazy );
            BOOST_CHECK( !runtime.status().initialized && runtime.status().owner == KokkosRuntimeOwner::none );
            if ( mode == "off" ) BOOST_CHECK_THROW( runtime.initialize(), std::logic_error );
        }
        else if ( !caps.available )
        {
            KokkosRuntime runtime;
            BOOST_CHECK_THROW( runtime.initialize(), std::logic_error );
        }
#ifdef FEELPP_HAS_KOKKOS
        else if ( mode == "owned" )
        {
            KokkosRuntime first;
            first.initialize();
            BOOST_CHECK( first.status().owner == KokkosRuntimeOwner::feelpp );
            {
                KokkosRuntime second( KokkosInitPolicy::eager );
                first.close();
                BOOST_CHECK( second.status().initialized && second.status().owner == KokkosRuntimeOwner::feelpp );
                // The pack is destroyed before the surviving runtime lease.
                Kokkos::View<double*, Kokkos::HostSpace> pack( "pack", 1 );
                pack( 0 ) = 7.;
                BOOST_CHECK( pack( 0 ) == 7. );
            }
            BOOST_CHECK( !Kokkos::is_initialized() && Kokkos::is_finalized() );
            KokkosRuntime restart;
            BOOST_CHECK_THROW( restart.initialize(), std::logic_error );
            BOOST_CHECK_THROW( first.initialize(), std::logic_error );
        }
        else if ( mode == "external" )
        {
            Kokkos::initialize();
            {
                KokkosRuntime first( KokkosInitPolicy::eager );
                KokkosRuntime second( KokkosInitPolicy::eager );
                BOOST_CHECK( first.status().owner == KokkosRuntimeOwner::external );
                first.close();
                BOOST_CHECK( second.status().initialized );
            }
            BOOST_CHECK( Kokkos::is_initialized() );
            Kokkos::finalize();
        }
        else if ( mode == "petsc" )
        {
            if ( caps.petsc )
            {
                BOOST_CHECK( PetscInitialize( &argc, &argv, nullptr, nullptr ) == 0 );
                {
                    KokkosRuntime runtime( KokkosInitPolicy::eager, KokkosRuntimeOwner::petsc );
                    BOOST_CHECK( runtime.status().owner == KokkosRuntimeOwner::petsc );
                }
                BOOST_CHECK( Kokkos::is_initialized() );
                BOOST_CHECK( PetscFinalize() == 0 );
            }
            else
            {
                KokkosRuntime runtime( KokkosInitPolicy::lazy, KokkosRuntimeOwner::petsc );
                BOOST_CHECK_THROW( runtime.initialize(), std::logic_error );
                BOOST_CHECK( !Kokkos::is_initialized() );
                std::cout << "PETSc Kokkos ownership unavailable; request rejected\n";
            }
        }
#endif
        else
            throw std::invalid_argument( "Unknown runtime test mode" );
        MPI_Finalize();
        return;
    }
    catch ( std::exception const& error )
    {
        BOOST_FAIL( error.what() );
    }
}

/** Verify the off runtime lifecycle in a fresh process. */
BOOST_AUTO_TEST_CASE( runtime_off )
{
    runRuntimeScenario( "off" );
}

/** Verify the lazy runtime lifecycle in a fresh process. */
BOOST_AUTO_TEST_CASE( runtime_lazy )
{
    runRuntimeScenario( "lazy" );
}

/** Verify the owned runtime lifecycle in a fresh process. */
BOOST_AUTO_TEST_CASE( runtime_owned )
{
    runRuntimeScenario( "owned" );
}

/** Verify the external runtime lifecycle in a fresh process. */
BOOST_AUTO_TEST_CASE( runtime_external )
{
    runRuntimeScenario( "external" );
}

/** Verify the petsc runtime lifecycle in a fresh process. */
BOOST_AUTO_TEST_CASE( runtime_petsc )
{
    runRuntimeScenario( "petsc" );
}

/** Verify the before-mpi runtime lifecycle in a fresh process. */
BOOST_AUTO_TEST_CASE( runtime_before_mpi )
{
    runRuntimeScenario( "before-mpi" );
}

/** Verify the after-mpi runtime lifecycle in a fresh process. */
BOOST_AUTO_TEST_CASE( runtime_after_mpi )
{
    runRuntimeScenario( "after-mpi" );
}
