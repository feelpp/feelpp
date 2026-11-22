// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <feel/feelconfig.h>
#include <feel/feelcore/environmentruntime.hpp>
#include <feel/feelcore/kokkosruntime.hpp>
#include <mpi.h>
#include <mutex>
#include <petscsys.h>
#include <stdexcept>
#ifdef FEELPP_HAS_KOKKOS
#include <Kokkos_Core.hpp>
#endif

namespace Feel
{
namespace
{
/// Serializes shared lease acquisition and final release within this library.
std::mutex runtimeMutex;
/// Observes a lease owned exclusively by Environment.
std::weak_ptr<KokkosRuntime> environmentLease;

/** Enforce initialization/finalization order and MPI main-thread use. */
void requireMpi()
{
    int initialized = 0, finalized = 0, mainThread = 0;
    MPI_Initialized( &initialized );
    MPI_Finalized( &finalized );
    if ( !initialized || finalized )
        throw std::logic_error( "Kokkos runtime requires active MPI" );
    MPI_Is_thread_main( &mainThread );
    if ( !mainThread )
        throw std::logic_error( "Kokkos runtime requires the MPI main thread" );
}
} // namespace

struct KokkosRuntime::State
{
    /// Authority established by the first successful acquisition.
    KokkosRuntimeOwner owner = KokkosRuntimeOwner::none;
    /** Fence and finalize only a runtime initialized by this service. */
    ~State()
    {
#ifdef FEELPP_HAS_KOKKOS
        if ( owner == KokkosRuntimeOwner::feelpp && Kokkos::is_initialized() )
        {
            requireMpi();
            Kokkos::fence();
            Kokkos::finalize();
        }
#endif
    }
};

KokkosRuntime::KokkosRuntime( KokkosInitPolicy policy, KokkosRuntimeOwner borrowedOwner )
    : M_policy( policy ), M_borrowedOwner( borrowedOwner )
{
    if ( borrowedOwner != KokkosRuntimeOwner::external && borrowedOwner != KokkosRuntimeOwner::petsc )
        throw std::invalid_argument( "Borrowed Kokkos owner must be external or PETSc" );
    if ( policy == KokkosInitPolicy::eager ) initialize();
}
KokkosRuntime::~KokkosRuntime()
{
    close();
}

void KokkosRuntime::configureEnvironment( KokkosInitPolicy policy, KokkosRuntimeOwner borrowedOwner )
{
    auto factory = [policy, borrowedOwner]()
    {
        auto lease = std::make_shared<KokkosRuntime>( policy, borrowedOwner );
        environmentLease = lease;
        return lease;
    };
    registerEnvironmentRuntime( "kokkos", std::move( factory ) );
}

KokkosRuntime& KokkosRuntime::environment()
{
    auto lease = environmentLease.lock();
    if ( !lease ) throw std::logic_error( "No live Environment-owned Kokkos runtime" );
    return *lease;
}

void KokkosRuntime::initialize()
{
    std::lock_guard lock( runtimeMutex );
    if ( M_closed ) throw std::logic_error( "Kokkos lease is closed" );
    if ( M_policy == KokkosInitPolicy::off ) throw std::logic_error( "Kokkos initialization is disabled" );
#ifdef FEELPP_HAS_KOKKOS
    requireMpi();
    if ( Kokkos::is_finalized() ) throw std::logic_error( "Kokkos cannot be restarted after finalization" );
    if ( M_state )
    {
        if ( !Kokkos::is_initialized() ) throw std::logic_error( "Kokkos owner finalized a live lease" );
        return;
    }
    // A weak registry lets the last lease release the runtime immediately.
    // Keeping the registry in this compiled library avoids per-translation-unit owners.
    static std::weak_ptr<State> active;
    if ( auto shared = active.lock() )
    {
        if ( !Kokkos::is_initialized() ) throw std::logic_error( "Kokkos owner finalized a live lease" );
        if ( M_borrowedOwner == KokkosRuntimeOwner::petsc && shared->owner != KokkosRuntimeOwner::petsc )
            throw std::logic_error( "Active Kokkos owner is not PETSc" );
        M_state = std::move( shared );
        return;
    }
    auto state = std::make_shared<State>();
    // Borrowing is explicit: observing an initialized runtime cannot identify
    // whether PETSc or an application initialized it first.
    if ( M_borrowedOwner == KokkosRuntimeOwner::petsc )
    {
        PetscBool initialized = PETSC_FALSE;
        if ( !capabilities().petsc || PetscInitialized( &initialized ) != 0 || !initialized || !Kokkos::is_initialized() )
            throw std::logic_error( "PETSc must own an initialized Kokkos runtime" );
        state->owner = KokkosRuntimeOwner::petsc;
    }
    else if ( Kokkos::is_initialized() )
        state->owner = KokkosRuntimeOwner::external;
    else
    {
        Kokkos::initialize();
        state->owner = KokkosRuntimeOwner::feelpp;
    }
    M_state = state;
    active = state;
#else
    throw std::logic_error( "Kokkos is unavailable in this build" );
#endif
}

void KokkosRuntime::close()
{
    std::lock_guard lock( runtimeMutex );
    M_state.reset();
    M_closed = true;
}

KokkosRuntimeStatus KokkosRuntime::status() const
{
    std::lock_guard lock( runtimeMutex );
    bool initialized = false, finalized = false;
#ifdef FEELPP_HAS_KOKKOS
    initialized = Kokkos::is_initialized();
    finalized = Kokkos::is_finalized();
#endif
    return { M_policy, M_state ? M_state->owner : KokkosRuntimeOwner::none, initialized, finalized };
}

KokkosCapabilities KokkosRuntime::capabilities()
{
#ifdef FEELPP_HAS_KOKKOS
#if defined( PETSC_HAVE_KOKKOS )
    constexpr bool petsc = true;
#else
    constexpr bool petsc = false;
#endif
    return { true, petsc, Kokkos::DefaultHostExecutionSpace::name(), Kokkos::DefaultExecutionSpace::name() };
#else
    return { false, false, "unavailable", "unavailable" };
#endif
}
} // namespace Feel
