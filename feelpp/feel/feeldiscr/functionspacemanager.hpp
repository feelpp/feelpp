/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

   SPDX-FileCopyrightText: 2026 University of Strasbourg

   SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
 * @file functionspacemanager.hpp
 * @brief Mesh-aware construction, reuse, invalidation, and observation of function spaces.
 */

#ifndef FEELPP_FUNCTIONSPACEMANAGER_HPP
#define FEELPP_FUNCTIONSPACEMANAGER_HPP 1

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <typeindex>

#include <boost/signals2/connection.hpp>

#include <feel/feeldiscr/enums.hpp>
#include <feel/feeldiscr/functionspacebase.hpp>
#include <feel/feelmesh/enums.hpp>
#include <feel/feelmesh/meshbase.hpp>

namespace Feel
{

/**
 * @brief Select how a whole-mesh function-space request interacts with the manager.
 */
enum class FunctionSpaceReusePolicy
{
    automatic, ///< Reuse when globally enabled; otherwise preserve always-new behavior.
    reuse,     ///< Return the retained space, or build and retain it on a cache miss.
    rebuild,   ///< Build a new space and replace the matching retained entry.
    bypass     ///< Build a new space without reading or modifying manager state.
};

/**
 * @brief Runtime configuration of the process-local function-space manager.
 *
 * Configuration that affects lookup or retention must be identical on every
 * rank of a distributed mesh communicator.
 */
struct FunctionSpaceManagerConfig
{
    bool enabled = false;                    ///< Enable reuse for @c automatic requests.
    std::size_t maxEntries = 64;             ///< Maximum entries retained in this process.
    std::size_t maxEntriesPerMesh = 16;      ///< Maximum entries retained for one mesh.
    bool mpiConsistencyDiagnostics = false;  ///< Log rank-local hit/miss disagreement.
};

/**
 * @brief Process-local counters and retained-state measurements.
 *
 * Counter values are cumulative until process shutdown. The retained values
 * are snapshots computed when @ref FunctionSpaceManager::stats is called.
 */
struct FunctionSpaceManagerStats
{
    uint64_type lookups = 0;             ///< Cacheable requests examined by the manager.
    uint64_type hits = 0;                ///< Requests returning an existing retained space.
    uint64_type misses = 0;              ///< Requests requiring construction.
    uint64_type builds = 0;              ///< Successful constructions made by the manager.
    uint64_type rebuilds = 0;            ///< Explicit @c rebuild requests.
    uint64_type evictions = 0;           ///< Entries removed by a capacity limit.
    uint64_type invalidations = 0;       ///< Entries removed after a structural mesh change.
    std::size_t retainedEntries = 0;     ///< Entries currently retained by the manager.
    uint64_type retainedDofs = 0;        ///< Sum of DOF counts for retained entries.
};

/**
 * @brief Normalized construction options that participate in a manager key.
 */
struct FunctionSpaceManagerOptions
{
    /// Effective extended-DOF-table mode.
    DofTableExtendedType extendedDofTable = DofTableExtendedType::VERTICES;
    /// Effective mesh update components.
    size_type components = MESH_RENUMBER | MESH_CHECK;
};

/**
 * @brief Own and reuse canonical whole-mesh function spaces within a process.
 *
 * Entries are keyed by mesh ownership identity, structural mesh revision,
 * concrete function-space type, and normalized construction options. The
 * manager holds strong references in bounded global and per-mesh LRU caches.
 * Coordinate-only mesh changes refresh live function spaces without replacing
 * their DOF topology.
 *
 * Calls for a distributed mesh are collective on its communicator: all ranks
 * must issue compatible requests in the same order. A cache hit is returned
 * only when every rank reports a local hit.
 */
class FEELPP_EXPORT FunctionSpaceManager
{
  public:
    /// Type-erased function-space pointer retained by the manager.
    using erased_space_ptrtype = std::shared_ptr<FunctionSpaceBase>;
    /// Callable that constructs a type-erased function space on a miss.
    using erased_factory_type = std::function<erased_space_ptrtype()>;
    /// Slot invoked when the associated mesh reports a change.
    using mesh_change_slot_type = std::function<void( MESH_CHANGES )>;
    /// Callable that connects a mesh-change slot and returns its connection.
    using mesh_connector_type = std::function<boost::signals2::connection( mesh_change_slot_type )>;

    /**
     * @brief Access the process-wide function-space manager.
     * @return the process-wide function-space manager instance
     */
    static FunctionSpaceManager& instance();

    /**
     * @brief Replace the manager configuration and enforce its capacity limits.
     * @param config new process-local configuration
     */
    void configure( FunctionSpaceManagerConfig config );

    /**
     * @brief Return the effective manager configuration.
     * @return the effective process-local manager configuration
     */
    FunctionSpaceManagerConfig config() const;

    /**
     * @brief Release entries retained for one mesh.
     *
     * Spaces that are still owned by callers remain valid.
     *
     * @tparam MeshType concrete mesh type
     * @param mesh mesh whose retained entries are released
     */
    template <typename MeshType>
    void clear( std::shared_ptr<MeshType> const& mesh )
    {
        this->clearOwner( meshOwner( mesh ) );
    }

    /**
     * @brief Release all entries retained by the manager.
     *
     * Spaces that are still owned by callers remain valid.
     */
    void clear();

    /**
     * @brief Return counters and retained-state measurements for one mesh.
     * @tparam MeshType concrete mesh type
     * @param mesh mesh for which statistics are requested
     * @return process-local statistics for @p mesh
     */
    template <typename MeshType>
    FunctionSpaceManagerStats stats( std::shared_ptr<MeshType> const& mesh ) const
    {
        return this->statsOwner( meshOwner( mesh ) );
    }

    /**
     * @brief Return aggregate statistics for all managed meshes.
     * @return aggregate process-local manager statistics
     */
    FunctionSpaceManagerStats stats() const;

    /**
     * @brief Return a compatible retained space or construct one according to a policy.
     *
     * This type-erased entry point is intended for function-space factory
     * helpers. Prefer @ref getOrCreateFunctionSpace from typed code.
     *
     * @param meshOwner weak ownership identity of the mesh
     * @param structuralRevision current function-space structural mesh revision
     * @param spaceType concrete requested function-space type
     * @param options normalized construction options included in the key
     * @param policy requested reuse policy
     * @param worldComm communicator on which request consistency is checked
     * @param factory construction callable used when a new space is required
     * @param connector optional mesh-change signal connector
     * @return the retained or newly constructed type-erased function space
     * @throws std::runtime_error if request signatures differ across MPI ranks
     */
    erased_space_ptrtype getOrCreate(
        std::weak_ptr<void> const& meshOwner,
        uint64_type structuralRevision,
        std::type_index spaceType,
        FunctionSpaceManagerOptions const& options,
        FunctionSpaceReusePolicy policy,
        WorldComm const& worldComm,
        erased_factory_type factory,
        mesh_connector_type connector );

  private:
    /// Private implementation containing keys, entries, counters, and locking.
    struct Impl;

    /// Construct the singleton and register its shutdown observer.
    FunctionSpaceManager();
    /// Release retained entries and mesh signal connections.
    ~FunctionSpaceManager();
    FunctionSpaceManager( FunctionSpaceManager const& ) = delete;
    FunctionSpaceManager& operator=( FunctionSpaceManager const& ) = delete;

    /**
     * @brief Convert a typed mesh pointer to ownership identity without retaining it.
     * @tparam MeshType concrete mesh type
     * @param mesh mesh whose shared ownership identity is requested
     * @return weak type-erased ownership identity
     */
    template <typename MeshType>
    static std::weak_ptr<void> meshOwner( std::shared_ptr<MeshType> const& mesh )
    {
        std::shared_ptr<void> owner = mesh;
        return std::weak_ptr<void>( owner );
    }

    /**
     * @brief Release entries associated with a type-erased mesh identity.
     * @param meshOwner weak mesh ownership identity
     */
    void clearOwner( std::weak_ptr<void> const& meshOwner );

    /**
     * @brief Return statistics for a type-erased mesh ownership identity.
     * @param meshOwner weak mesh ownership identity
     * @return process-local statistics associated with @p meshOwner
     */
    FunctionSpaceManagerStats statsOwner( std::weak_ptr<void> const& meshOwner ) const;

    /**
     * @brief Refresh live function spaces after a coordinate-only mesh change.
     * @param meshOwner weak mesh ownership identity
     * @param changes kind of mesh change reported by the mesh
     */
    void refreshGeometry( std::weak_ptr<void> const& meshOwner, MESH_CHANGES changes );

    /// Optionally log statistics and release all manager state exactly once.
    void shutdown();

    /// Opaque manager implementation.
    std::unique_ptr<Impl> M_impl;
};

/**
 * @brief Normalize an extended-DOF-table option for use in a manager key.
 * @param dte requested extended-DOF-table mode
 * @return @c VERTICES for @c DEFAULT, otherwise @p dte
 */
inline DofTableExtendedType
normalizeFunctionSpaceDofTable( DofTableExtendedType dte )
{
    return dte == DofTableExtendedType::DEFAULT ? DofTableExtendedType::VERTICES : dte;
}

/**
 * @brief Typed factory adapter for the process-wide function-space manager.
 *
 * Mesh-change observation is connected when the mesh type exposes the
 * expected @c meshChanged signal. The returned pointer is checked against the
 * requested concrete type.
 *
 * @tparam SpaceType concrete function-space type
 * @tparam MeshType concrete mesh type
 * @tparam Factory construction callable type
 * @param mesh mesh on which the space is defined
 * @param options normalized options that participate in the manager key
 * @param policy requested reuse policy
 * @param factory callable constructing @c std::shared_ptr<SpaceType>
 * @return retained or newly constructed function space
 */
template <typename SpaceType, typename MeshType, typename Factory>
std::shared_ptr<SpaceType>
getOrCreateFunctionSpace( std::shared_ptr<MeshType> const& mesh,
                          FunctionSpaceManagerOptions const& options,
                          FunctionSpaceReusePolicy policy,
                          Factory&& factory )
{
    std::shared_ptr<void> owner = mesh;
    auto erasedFactory = [factory = std::forward<Factory>( factory )]() mutable
    {
        return std::static_pointer_cast<FunctionSpaceBase>( factory() );
    };

    FunctionSpaceManager::mesh_connector_type connector;
#if !defined( __INTEL_COMPILER )
    if constexpr ( requires( MeshType& candidate,
                             FunctionSpaceManager::mesh_change_slot_type slot )
                   { candidate.meshChanged.connect( slot ); } )
    {
        std::weak_ptr<MeshType> weakMesh = mesh;
        connector = [weakMesh]( FunctionSpaceManager::mesh_change_slot_type slot )
        {
            if ( auto lockedMesh = weakMesh.lock() )
                return lockedMesh->meshChanged.connect( std::move( slot ) );
            return boost::signals2::connection{};
        };
    }
#endif

    auto erased = FunctionSpaceManager::instance().getOrCreate(
        std::weak_ptr<void>( owner ), mesh->functionSpaceStructuralRevision(),
        std::type_index( typeid( SpaceType ) ), options, policy, mesh->worldComm(),
        std::move( erasedFactory ), std::move( connector ) );
    auto result = std::dynamic_pointer_cast<SpaceType>( erased );
    CHECK( result ) << "function-space manager returned an incompatible type";
    return result;
}

} // namespace Feel

#endif // FEELPP_FUNCTIONSPACEMANAGER_HPP
