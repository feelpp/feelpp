/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

   SPDX-FileCopyrightText: 2026 University of Strasbourg

   SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <feel/feeldiscr/functionspacemanager.hpp>

#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/operations.hpp>

#include <feel/feelcore/environment.hpp>

namespace Feel
{
namespace
{
using owner_type = std::weak_ptr<void>;
using owner_less_type = std::owner_less<owner_type>;

/**
 * @brief Mix one value into a 64-bit request signature.
 * @param seed current signature
 * @param value value to incorporate
 * @return combined signature
 */
uint64_type hashCombine( uint64_type seed, uint64_type value )
{
    seed ^= value + UINT64_C( 0x9e3779b97f4a7c15 ) + ( seed << 6 ) + ( seed >> 2 );
    return seed;
}

/**
 * @brief Compute a deterministic 64-bit hash of a null-terminated string.
 * @param value string to hash
 * @return FNV-1a hash of @p value
 */
uint64_type hashString( char const* value )
{
    uint64_type hash = UINT64_C( 1469598103934665603 );
    for ( ; *value; ++value )
    {
        hash ^= static_cast<unsigned char>( *value );
        hash *= UINT64_C( 1099511628211 );
    }
    return hash;
}

} // namespace

/**
 * @brief Private synchronized storage and policy implementation.
 */
struct FunctionSpaceManager::Impl
{
    /** @brief Key for one function space associated with a mesh bucket. */
    struct Key
    {
        uint64_type structuralRevision;             ///< Structural mesh revision.
        std::type_index spaceType;                  ///< Concrete function-space type.
        DofTableExtendedType extendedDofTable;      ///< Normalized DOF-table mode.
        size_type components;                       ///< Effective mesh components.

        /** @return @c true when this key precedes @p other. */
        bool operator<( Key const& other ) const
        {
            return std::tie( structuralRevision, spaceType, extendedDofTable, components ) <
                   std::tie( other.structuralRevision, other.spaceType,
                             other.extendedDofTable, other.components );
        }
    };

    /** @brief Strongly retained function space and its LRU metadata. */
    struct Entry
    {
        erased_space_ptrtype space; ///< Retained function space.
        uint64_type dofs = 0;       ///< DOF count used for retained-state reporting.
        uint64_type lastAccess = 0; ///< Monotonic LRU access stamp.
    };

    /** @brief Entries, observers, and counters associated with one mesh owner. */
    struct MeshBucket
    {
        std::map<Key, Entry> entries; ///< Retained entries keyed by construction identity.
        std::vector<std::weak_ptr<FunctionSpaceBase>> liveSpaces; ///< Spaces requiring geometry refresh.
        boost::signals2::connection geometryConnection; ///< Mesh-change signal connection.
        FunctionSpaceManagerStats counters; ///< Cumulative counters for this mesh.
    };

    using buckets_type = std::map<owner_type, MeshBucket, owner_less_type>;

    mutable std::recursive_mutex mutex;
    FunctionSpaceManagerConfig config;
    bool runtimeConfigLoaded = false;
    bool logStats = false;
    bool shutdownComplete = false;
    uint64_type accessClock = 0;
    buckets_type buckets;
    FunctionSpaceManagerStats counters;

    /** @brief Load manager options from the Environment exactly once. */
    void loadRuntimeConfigOnce()
    {
        if ( runtimeConfigLoaded )
            return;
        runtimeConfigLoaded = true;
        auto const& vm = Environment::vm();
        if ( vm.count( "functionspace.manager.enable" ) )
            config.enabled = vm["functionspace.manager.enable"].as<bool>();
        if ( vm.count( "functionspace.manager.max-entries" ) )
            config.maxEntries = vm["functionspace.manager.max-entries"].as<std::size_t>();
        if ( vm.count( "functionspace.manager.max-entries-per-mesh" ) )
            config.maxEntriesPerMesh = vm["functionspace.manager.max-entries-per-mesh"].as<std::size_t>();
        if ( vm.count( "functionspace.manager.mpi-consistency-diagnostics" ) )
            config.mpiConsistencyDiagnostics =
                vm["functionspace.manager.mpi-consistency-diagnostics"].as<bool>();
        if ( vm.count( "functionspace.manager.log-stats" ) )
            logStats = vm["functionspace.manager.log-stats"].as<bool>();
    }

    /**
     * @brief Increment the same counter globally and for one mesh.
     * @param bucket mesh bucket whose counter is incremented
     * @param member counter member to update
     * @param amount increment amount
     */
    void increment( MeshBucket& bucket, uint64_type FunctionSpaceManagerStats::*member,
                    uint64_type amount = 1 )
    {
        counters.*member += amount;
        bucket.counters.*member += amount;
    }

    /** @return number of entries retained across all mesh buckets. */
    std::size_t retainedEntries() const
    {
        std::size_t result = 0;
        for ( auto const& [owner, bucket] : buckets )
            result += bucket.entries.size();
        return result;
    }

    /**
     * @brief Erase one entry and optionally account for an eviction.
     * @param bucket bucket containing @p entry
     * @param entry iterator of the entry to erase
     * @param eviction whether to increment the eviction counter
     */
    void eraseEntry( MeshBucket& bucket, typename std::map<Key, Entry>::iterator entry,
                     bool eviction )
    {
        bucket.entries.erase( entry );
        if ( eviction )
            increment( bucket, &FunctionSpaceManagerStats::evictions );
    }

    /**
     * @brief Enforce per-mesh and global LRU capacity limits.
     * @param insertedBucket bucket in which an entry was inserted
     */
    void evictIfNeeded( MeshBucket& insertedBucket )
    {
        while ( insertedBucket.entries.size() > config.maxEntriesPerMesh )
        {
            auto oldest = std::min_element(
                insertedBucket.entries.begin(), insertedBucket.entries.end(),
                []( auto const& left, auto const& right )
                { return left.second.lastAccess < right.second.lastAccess; } );
            eraseEntry( insertedBucket, oldest, true );
        }

        while ( retainedEntries() > config.maxEntries )
        {
            auto oldestBucket = buckets.end();
            typename std::map<Key, Entry>::iterator oldestEntry;
            uint64_type oldestAccess = std::numeric_limits<uint64_type>::max();
            for ( auto bucketIt = buckets.begin(); bucketIt != buckets.end(); ++bucketIt )
            {
                for ( auto entryIt = bucketIt->second.entries.begin();
                      entryIt != bucketIt->second.entries.end(); ++entryIt )
                {
                    if ( entryIt->second.lastAccess < oldestAccess )
                    {
                        oldestAccess = entryIt->second.lastAccess;
                        oldestBucket = bucketIt;
                        oldestEntry = entryIt;
                    }
                }
            }
            if ( oldestBucket == buckets.end() )
                break;
            eraseEntry( oldestBucket->second, oldestEntry, true );
        }
    }

    /** @brief Remove empty buckets whose mesh owner has expired. */
    void pruneExpiredBuckets()
    {
        for ( auto it = buckets.begin(); it != buckets.end(); )
        {
            if ( it->first.expired() && it->second.entries.empty() )
                it = buckets.erase( it );
            else
                ++it;
        }
    }
};

FunctionSpaceManager& FunctionSpaceManager::instance()
{
    static FunctionSpaceManager manager;
    return manager;
}

FunctionSpaceManager::FunctionSpaceManager()
    : M_impl( std::make_unique<Impl>() )
{
    Environment::addDeleteObserver( [this]() { this->shutdown(); } );
}

FunctionSpaceManager::~FunctionSpaceManager()
{
    this->shutdown();
}

void FunctionSpaceManager::configure( FunctionSpaceManagerConfig config )
{
    std::lock_guard lock( M_impl->mutex );
    M_impl->config = config;
    M_impl->runtimeConfigLoaded = true;
    for ( auto& [owner, bucket] : M_impl->buckets )
        M_impl->evictIfNeeded( bucket );
}

FunctionSpaceManagerConfig FunctionSpaceManager::config() const
{
    std::lock_guard lock( M_impl->mutex );
    M_impl->loadRuntimeConfigOnce();
    return M_impl->config;
}

FunctionSpaceManager::erased_space_ptrtype FunctionSpaceManager::getOrCreate(
    std::weak_ptr<void> const& meshOwner, uint64_type structuralRevision,
    std::type_index spaceType, FunctionSpaceManagerOptions const& options,
    FunctionSpaceReusePolicy policy, WorldComm const& worldComm,
    erased_factory_type factory, mesh_connector_type connector )
{
    std::lock_guard lock( M_impl->mutex );
    M_impl->loadRuntimeConfigOnce();

    auto const normalizedDofTable = normalizeFunctionSpaceDofTable( options.extendedDofTable );
    auto resolvedPolicy = policy;
    if ( resolvedPolicy == FunctionSpaceReusePolicy::automatic )
        resolvedPolicy = M_impl->config.enabled ? FunctionSpaceReusePolicy::reuse
                                                : FunctionSpaceReusePolicy::bypass;

    uint64_type requestSignature = hashString( spaceType.name() );
    requestSignature = hashCombine( requestSignature, static_cast<uint64_type>( resolvedPolicy ) );
    requestSignature = hashCombine( requestSignature, structuralRevision );
    requestSignature = hashCombine( requestSignature, static_cast<uint64_type>( normalizedDofTable ) );
    requestSignature = hashCombine( requestSignature, options.components );
    requestSignature = hashCombine( requestSignature, M_impl->config.enabled );
    requestSignature = hashCombine( requestSignature, M_impl->config.maxEntries );
    requestSignature = hashCombine( requestSignature, M_impl->config.maxEntriesPerMesh );
    requestSignature = hashCombine( requestSignature, M_impl->config.mpiConsistencyDiagnostics );

    auto const& comm = worldComm.globalComm();
    if ( comm.size() > 1 )
    {
        auto const minimumSignature = boost::mpi::all_reduce(
            comm, requestSignature, boost::mpi::minimum<uint64_type>() );
        auto const maximumSignature = boost::mpi::all_reduce(
            comm, requestSignature, boost::mpi::maximum<uint64_type>() );
        if ( minimumSignature != maximumSignature )
            throw std::runtime_error(
                "function-space manager request differs across the mesh communicator" );
    }

    if ( resolvedPolicy == FunctionSpaceReusePolicy::bypass )
        return factory();

    M_impl->pruneExpiredBuckets();
    auto [bucketIt, insertedBucket] = M_impl->buckets.try_emplace( meshOwner );
    auto& bucket = bucketIt->second;
    if ( insertedBucket && connector )
    {
        bucket.geometryConnection = connector(
            [this, meshOwner]( MESH_CHANGES changes )
            { this->refreshGeometry( meshOwner, changes ); } );
    }

    for ( auto entryIt = bucket.entries.begin(); entryIt != bucket.entries.end(); )
    {
        if ( entryIt->first.structuralRevision != structuralRevision )
        {
            entryIt = bucket.entries.erase( entryIt );
            M_impl->increment( bucket, &FunctionSpaceManagerStats::invalidations );
        }
        else
            ++entryIt;
    }

    Impl::Key key{ structuralRevision, spaceType, normalizedDofTable, options.components };
    M_impl->increment( bucket, &FunctionSpaceManagerStats::lookups );
    auto found = bucket.entries.find( key );
    bool const localHit = found != bucket.entries.end();

    bool allHit = localHit;
    bool anyHit = localHit;
    if ( comm.size() > 1 && resolvedPolicy == FunctionSpaceReusePolicy::reuse )
    {
        allHit = boost::mpi::all_reduce( comm, localHit, std::logical_and<bool>() );
        if ( M_impl->config.mpiConsistencyDiagnostics )
            anyHit = boost::mpi::all_reduce( comm, localHit, std::logical_or<bool>() );
    }

    if ( resolvedPolicy == FunctionSpaceReusePolicy::reuse && allHit )
    {
        M_impl->increment( bucket, &FunctionSpaceManagerStats::hits );
        found->second.lastAccess = ++M_impl->accessClock;
        return found->second.space;
    }

    if ( M_impl->config.mpiConsistencyDiagnostics && anyHit != allHit && worldComm.isMasterRank() )
        LOG( WARNING ) << "function-space manager local hit mismatch; rebuilding on all ranks";

    M_impl->increment( bucket, &FunctionSpaceManagerStats::misses );
    if ( resolvedPolicy == FunctionSpaceReusePolicy::rebuild )
        M_impl->increment( bucket, &FunctionSpaceManagerStats::rebuilds );

    auto built = factory();
    M_impl->increment( bucket, &FunctionSpaceManagerStats::builds );
    auto const dofs = built->mapPtr() ? static_cast<uint64_type>( built->mapPtr()->nDof() ) : 0;
    bucket.entries.insert_or_assign( key, Impl::Entry{ built, dofs, ++M_impl->accessClock } );
    auto const& vm = Environment::vm();
    bool const legacyMeshConnection = vm.count( "connect" ) && vm["connect"].as<bool>();
    if ( !legacyMeshConnection )
        bucket.liveSpaces.emplace_back( built );
    M_impl->evictIfNeeded( bucket );
    return built;
}

void FunctionSpaceManager::clearOwner( std::weak_ptr<void> const& meshOwner )
{
    std::lock_guard lock( M_impl->mutex );
    auto found = M_impl->buckets.find( meshOwner );
    if ( found != M_impl->buckets.end() )
        found->second.entries.clear();
}

void FunctionSpaceManager::clear()
{
    std::lock_guard lock( M_impl->mutex );
    for ( auto& [owner, bucket] : M_impl->buckets )
        bucket.entries.clear();
    M_impl->pruneExpiredBuckets();
}

FunctionSpaceManagerStats FunctionSpaceManager::statsOwner(
    std::weak_ptr<void> const& meshOwner ) const
{
    std::lock_guard lock( M_impl->mutex );
    auto found = M_impl->buckets.find( meshOwner );
    if ( found == M_impl->buckets.end() )
        return {};
    auto result = found->second.counters;
    result.retainedEntries = found->second.entries.size();
    for ( auto const& [key, entry] : found->second.entries )
        result.retainedDofs += entry.dofs;
    return result;
}

FunctionSpaceManagerStats FunctionSpaceManager::stats() const
{
    std::lock_guard lock( M_impl->mutex );
    auto result = M_impl->counters;
    for ( auto const& [owner, bucket] : M_impl->buckets )
    {
        result.retainedEntries += bucket.entries.size();
        for ( auto const& [key, entry] : bucket.entries )
            result.retainedDofs += entry.dofs;
    }
    return result;
}

void FunctionSpaceManager::refreshGeometry( std::weak_ptr<void> const& meshOwner,
                                               MESH_CHANGES changes )
{
    if ( changes != MESH_CHANGES_POINTS_COORDINATES )
        return;

    std::vector<std::shared_ptr<FunctionSpaceBase>> spaces;
    {
        std::lock_guard lock( M_impl->mutex );
        auto found = M_impl->buckets.find( meshOwner );
        if ( found == M_impl->buckets.end() )
            return;
        auto& liveSpaces = found->second.liveSpaces;
        for ( auto it = liveSpaces.begin(); it != liveSpaces.end(); )
        {
            if ( auto space = it->lock() )
            {
                spaces.push_back( std::move( space ) );
                ++it;
            }
            else
                it = liveSpaces.erase( it );
        }
    }
    for ( auto const& space : spaces )
        space->updateAfterMeshChange( changes );
}

void FunctionSpaceManager::shutdown()
{
    std::lock_guard lock( M_impl->mutex );
    if ( M_impl->shutdownComplete )
        return;
    M_impl->shutdownComplete = true;
    if ( M_impl->logStats && Environment::isMasterRank() )
    {
        auto const currentStats = this->stats();
        std::cout << "function-space manager: lookups=" << currentStats.lookups
                  << " hits=" << currentStats.hits
                  << " misses=" << currentStats.misses
                  << " builds=" << currentStats.builds
                  << " evictions=" << currentStats.evictions << '\n';
    }
    M_impl->buckets.clear();
}

} // namespace Feel
