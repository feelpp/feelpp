// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <Kokkos_Core.hpp>
#include <cmath>
#define BOOST_TEST_MODULE kokkos_smoke
#include <boost/test/unit_test.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/kokkosruntime.hpp>
#include <iostream>
#include <mpi.h>

/** Fill and reduce host data, reporting the execution space actually used. */
BOOST_AUTO_TEST_CASE( host_smoke )
{
    Feel::KokkosRuntime::configureEnvironment( Feel::KokkosInitPolicy::eager );
    int argc = 1;
    char app[] = "kokkos_host_smoke";
    char* args[] = { app, nullptr };
    char** argv = args;
    {
        Feel::Environment env( Feel::_argc = argc, Feel::_argv = argv, Feel::_desc = Feel::feel_nooptions() );
        auto& runtime = Feel::KokkosRuntime::environment();
        // Host execution is explicit even when the configured default is a GPU.
        // Views declared after the lease are destroyed before runtime finalization.
        using ExecutionSpace = Kokkos::DefaultHostExecutionSpace;
        Kokkos::View<double*, ExecutionSpace::memory_space> values( "values", 1024 );
        Kokkos::parallel_for( "fill", Kokkos::RangePolicy<ExecutionSpace>( 0, 1024 ), KOKKOS_LAMBDA( int i ) { values( i ) = 1.; } );
        double sum = 0;
        Kokkos::parallel_reduce( "sum", Kokkos::RangePolicy<ExecutionSpace>( 0, 1024 ), KOKKOS_LAMBDA( int i, double& local ) { local += values( i ); }, sum );
        ExecutionSpace().fence();
        std::cout << "requested=host actual=" << ExecutionSpace::name()
                  << " owner=" << ( runtime.status().owner == Feel::KokkosRuntimeOwner::feelpp ? "feelpp" : "borrowed" )
                  << " sum=" << sum << '\n';
        BOOST_TEST( std::abs( sum - 1024. ) < 1e-12 );
    }
}
