// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-3.0-or-later
#define BOOST_TEST_MODULE kokkos_environment
#include "runtimefixture.hpp"
#include <boost/test/unit_test.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/kokkosruntime.hpp>
#ifdef FEELPP_HAS_KOKKOS
#include <Kokkos_Core.hpp>
#endif

/** Exercise the selected policy through the real Environment lifecycle. */
void checkEnvironment( Feel::KokkosInitPolicy policy )
{
    using namespace Feel;
    KokkosRuntime::configureEnvironment( policy );
    BOOST_CHECK_THROW( (void)KokkosRuntime::environment(), std::logic_error );
    int argc = 1;
    char app[] = "kokkos_environment";
    char* args[] = { app, nullptr };
    char** argv = args;
    {
        Environment env( _argc = argc, _argv = argv, _desc = feel_nooptions() );
        auto& runtime = KokkosRuntime::environment();
        BOOST_CHECK( runtime.status().requested == policy );
        BOOST_CHECK_THROW( KokkosRuntime::configureEnvironment( policy ), std::logic_error );
        if ( policy == KokkosInitPolicy::off )
        {
            BOOST_CHECK( !runtime.status().initialized );
            BOOST_CHECK_THROW( runtime.initialize(), std::logic_error );
        }
        else if ( policy == KokkosInitPolicy::lazy )
        {
            BOOST_CHECK( !runtime.status().initialized );
            if ( KokkosRuntime::capabilities().available )
                runtime.initialize();
            else
                BOOST_CHECK_THROW( runtime.initialize(), std::logic_error );
        }
        else
        {
            BOOST_CHECK( runtime.status().initialized );
            BOOST_CHECK( runtime.status().owner == KokkosRuntimeOwner::feelpp );
        }
    }
    BOOST_CHECK_THROW( (void)KokkosRuntime::environment(), std::logic_error );
#ifdef FEELPP_HAS_KOKKOS
    BOOST_TEST( Kokkos::is_finalized() == ( policy != KokkosInitPolicy::off ) );
#endif
}

/** Environment off policy must preserve CPU-only execution. */
BOOST_AUTO_TEST_CASE( environment_off )
{
    checkEnvironment( Feel::KokkosInitPolicy::off );
}
/** Environment lazy policy initializes only on an explicit execution request. */
BOOST_AUTO_TEST_CASE( environment_lazy )
{
    checkEnvironment( Feel::KokkosInitPolicy::lazy );
}
#ifdef FEELPP_HAS_KOKKOS
/** Environment owns eager initialization and finalizes before MPI. */
BOOST_AUTO_TEST_CASE( environment_owned )
{
    checkEnvironment( Feel::KokkosInitPolicy::eager );
}
/** Environment must leave externally initialized Kokkos alive. */
BOOST_AUTO_TEST_CASE( environment_external )
{
    using namespace Feel;
    ScopedMpi mpi;
    Kokkos::initialize();
    KokkosRuntime::configureEnvironment( KokkosInitPolicy::eager );
    int argc = 1;
    char app[] = "kokkos_external_environment";
    char* args[] = { app, nullptr };
    char** argv = args;
    {
        Environment env( _argc = argc, _argv = argv, _desc = feel_nooptions() );
        BOOST_CHECK( KokkosRuntime::environment().status().owner == KokkosRuntimeOwner::external );
    }
    BOOST_TEST( Kokkos::is_initialized() );
    Kokkos::finalize();
}
#endif
