// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <memory>
#include <string>

namespace Feel
{
/** When a lease may initialize Kokkos. */
enum class KokkosInitPolicy
{
    off,  ///< Reject initialization, including explicit requests.
    lazy, ///< Initialize on the first explicit request.
    eager ///< Initialize during lease construction.
};

/** Actual finalization authority of the shared runtime. */
enum class KokkosRuntimeOwner
{
    none,     ///< This lease has not acquired the runtime.
    feelpp,   ///< The last Feel++ lease fences and finalizes Kokkos.
    external, ///< The external caller retains finalization authority.
    petsc     ///< PETSc retains finalization authority.
};

/** Compiled capabilities, independent of initialization or kernel execution. */
struct KokkosCapabilities
{
    bool available;                    ///< Whether this library supports Kokkos.
    bool petsc;                        ///< Whether PETSc was built with Kokkos.
    std::string hostExecutionSpace;    ///< Explicit host space compiled into Kokkos.
    std::string defaultExecutionSpace; ///< Default space; does not imply its use.
};

/** Requested initialization policy and observed process runtime state. */
struct KokkosRuntimeStatus
{
    KokkosInitPolicy requested; ///< Policy selected for this lease.
    KokkosRuntimeOwner owner;   ///< Authority associated with this acquired lease.
    bool initialized;           ///< Whether the process runtime is currently active.
    bool finalized;             ///< Whether the process runtime has been finalized.
};

/** Scoped runtime lease. Create, initialize and release on the MPI main thread.
 * MPI must outlive every lease; destroy execution packs before their lease.
 * Off and unused lazy leases never initialize Kokkos. External owners must
 * keep their runtime alive until all borrowed leases have been released.
 * Status describes runtime state, not execution of a kernel.
 */
class KokkosRuntime
{
  public:
    /** Create a lease, initializing only for the eager policy.
     * @param policy Initialization timing; lazy is the default.
     * @param borrowedOwner Authority to require when borrowing an existing runtime.
     * PETSc must already own initialized Kokkos when selected explicitly.
     * @throws std::invalid_argument For an invalid borrowing authority.
     * @throws std::logic_error If eager initialization cannot satisfy the contract.
     */
    explicit KokkosRuntime( KokkosInitPolicy policy = KokkosInitPolicy::lazy,
                            KokkosRuntimeOwner borrowedOwner = KokkosRuntimeOwner::external );

    /** Release this lease; the last owned lease finalizes before MPI teardown. */
    ~KokkosRuntime();

    /** Leases have explicit, noncopyable lifetimes. */
    KokkosRuntime( KokkosRuntime const& ) = delete;

    /** Assignment cannot implicitly replace a live lease. */
    KokkosRuntime& operator=( KokkosRuntime const& ) = delete;

    /** Acquire the shared runtime after MPI initialization on its main thread.
     * Repeated calls on an active lease are harmless. A closed lease or finalized
     * runtime cannot be restarted. Acquisition never initializes PETSc itself.
     * @throws std::logic_error For unavailable/disabled runtime or invalid lifecycle.
     */
    void initialize();

    /** Irreversibly release this lease; repeated calls are harmless.
     * Destroy all packs using this lease first. When Feel++ owns the runtime,
     * its last lease fences outstanding work and finalizes while MPI is active.
     */
    void close();

    /** Observe policy, lease ownership and process lifecycle without initialization. */
    [[nodiscard]] KokkosRuntimeStatus status() const;

    /** Report build capabilities without initializing Kokkos or querying a device. */
    [[nodiscard]] static KokkosCapabilities capabilities();

    /** Select the Kokkos policy before constructing Environment.
     * Environment owns the lease from MPI/solver setup until before teardown.
     * @param policy Initialization timing; lazy preserves CPU-only use.
     * @param borrowedOwner Explicit authority when borrowing, including PETSc-owned Kokkos.
     * @throws std::logic_error For duplicate or late configuration.
     */
    static void configureEnvironment( KokkosInitPolicy policy = KokkosInitPolicy::lazy,
                                      KokkosRuntimeOwner borrowedOwner = KokkosRuntimeOwner::external );

    /** Access the lease owned by the live Environment.
     * The reference expires at Environment teardown. Destroy execution packs first.
     * @throws std::logic_error Before configured setup or after teardown.
     */
    [[nodiscard]] static KokkosRuntime& environment();

  private:
    /** Shared finalization authority across simultaneously acquired leases. */
    struct State;

    std::shared_ptr<State> M_state;     ///< Empty until explicit acquisition.
    KokkosInitPolicy M_policy;          ///< Immutable per-lease initialization policy.
    KokkosRuntimeOwner M_borrowedOwner; ///< Explicit authority required when borrowing.

    bool M_closed = false; ///< Prevents acquisition after release.
};
} // namespace Feel
