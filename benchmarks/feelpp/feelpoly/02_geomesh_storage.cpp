/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file 02_geomesh_storage.cpp
 * @brief Memory, allocation, traversal, partition, and persistence baselines
 *        for geometric mesh entities.
 */

#include <benchmark/benchmark.h>

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelmesh/partitionmesh.hpp>

#include <boost/archive/binary_oarchive.hpp>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <memory>
#include <new>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>

#if defined( __GLIBC__ )
#include <malloc.h>
#endif

#if defined( __linux__ )
#include <unistd.h>
#endif

namespace geomesh_allocation_probe
{
/**
 * @struct AllocationSnapshot
 * @brief Cumulative C++ allocation activity captured by the benchmark.
 */
struct AllocationSnapshot
{
    /** Number of calls to the global C++ allocation functions. */
    std::uint64_t calls = 0;
    /** Sum of requested bytes passed to the global C++ allocation functions. */
    std::uint64_t bytes = 0;
};

/** Whether allocation recording is currently enabled. */
std::atomic<bool> enabled{ false };
/** Cumulative allocation-call counter. */
std::atomic<std::uint64_t> calls{ 0 };
/** Cumulative requested-byte counter. */
std::atomic<std::uint64_t> bytes{ 0 };

/**
 * @brief Record one C++ allocation request when probing is enabled.
 * @param size requested allocation size in bytes
 */
inline void record( std::size_t size ) noexcept
{
    if ( !enabled.load( std::memory_order_relaxed ) )
        return;
    calls.fetch_add( 1, std::memory_order_relaxed );
    bytes.fetch_add( size, std::memory_order_relaxed );
}

/** @brief Reset and enable allocation recording. */
inline void begin() noexcept
{
    calls.store( 0, std::memory_order_relaxed );
    bytes.store( 0, std::memory_order_relaxed );
    enabled.store( true, std::memory_order_release );
}

/**
 * @brief Disable allocation recording and return the cumulative snapshot.
 * @return recorded C++ allocation calls and requested bytes
 */
[[nodiscard]] inline AllocationSnapshot end() noexcept
{
    enabled.store( false, std::memory_order_release );
    return { calls.load( std::memory_order_relaxed ),
             bytes.load( std::memory_order_relaxed ) };
}

/**
 * @brief Allocate storage with the requested alignment.
 * @param size requested byte count
 * @param alignment requested alignment
 * @return allocated storage
 * @throws std::bad_alloc when allocation fails
 */
[[nodiscard]] void* allocate( std::size_t size, std::size_t alignment )
{
    record( size );
    void* storage = nullptr;
    if ( alignment <= alignof( std::max_align_t ) )
        storage = std::malloc( size );
    else if ( ::posix_memalign( &storage, alignment, size ) != 0 )
        storage = nullptr;
    if ( !storage )
        throw std::bad_alloc{};
    return storage;
}
} // namespace geomesh_allocation_probe

/**
 * @brief Global scalar allocation hook used only by this benchmark executable.
 * @param size requested byte count
 * @return allocated storage
 */
void* operator new( std::size_t size )
{
    return geomesh_allocation_probe::allocate(
        size, alignof( std::max_align_t ) );
}

/**
 * @brief Global array allocation hook used only by this benchmark executable.
 * @param size requested byte count
 * @return allocated storage
 */
void* operator new[]( std::size_t size )
{
    return geomesh_allocation_probe::allocate(
        size, alignof( std::max_align_t ) );
}

/** @brief Release scalar storage allocated by the benchmark hook. */
void operator delete( void* storage ) noexcept
{
    std::free( storage );
}

/** @brief Release array storage allocated by the benchmark hook. */
void operator delete[]( void* storage ) noexcept
{
    std::free( storage );
}

/** @brief Release sized scalar storage allocated by the benchmark hook. */
void operator delete( void* storage, std::size_t ) noexcept
{
    std::free( storage );
}

/** @brief Release sized array storage allocated by the benchmark hook. */
void operator delete[]( void* storage, std::size_t ) noexcept
{
    std::free( storage );
}

/**
 * @brief Global aligned scalar allocation hook.
 * @param size requested byte count
 * @param alignment requested alignment
 * @return allocated aligned storage
 */
void* operator new( std::size_t size, std::align_val_t alignment )
{
    return geomesh_allocation_probe::allocate(
        size, static_cast<std::size_t>( alignment ) );
}

/**
 * @brief Global aligned array allocation hook.
 * @param size requested byte count
 * @param alignment requested alignment
 * @return allocated aligned storage
 */
void* operator new[]( std::size_t size, std::align_val_t alignment )
{
    return geomesh_allocation_probe::allocate(
        size, static_cast<std::size_t>( alignment ) );
}

/** @brief Release aligned scalar storage. */
void operator delete( void* storage, std::align_val_t ) noexcept
{
    std::free( storage );
}

/** @brief Release aligned array storage. */
void operator delete[]( void* storage, std::align_val_t ) noexcept
{
    std::free( storage );
}

/** @brief Release sized aligned scalar storage. */
void operator delete(
    void* storage, std::size_t, std::align_val_t ) noexcept
{
    std::free( storage );
}

/** @brief Release sized aligned array storage. */
void operator delete[](
    void* storage, std::size_t, std::align_val_t ) noexcept
{
    std::free( storage );
}

/**
 * @brief Nothrow scalar allocation hook.
 * @param size requested byte count
 * @return allocated storage or `nullptr`
 */
void* operator new( std::size_t size, std::nothrow_t const& ) noexcept
{
    try
    {
        return ::operator new( size );
    }
    catch ( ... )
    {
        return nullptr;
    }
}

/**
 * @brief Nothrow array allocation hook.
 * @param size requested byte count
 * @return allocated storage or `nullptr`
 */
void* operator new[]( std::size_t size, std::nothrow_t const& ) noexcept
{
    try
    {
        return ::operator new[]( size );
    }
    catch ( ... )
    {
        return nullptr;
    }
}

/** @brief Release nothrow scalar storage. */
void operator delete(
    void* storage, std::nothrow_t const& ) noexcept
{
    std::free( storage );
}

/** @brief Release nothrow array storage. */
void operator delete[](
    void* storage, std::nothrow_t const& ) noexcept
{
    std::free( storage );
}

using namespace Feel;

namespace
{
/**
 * @struct MeshStorageStatistics
 * @brief Entity counts and conservative storage lower bounds for one mesh.
 */
struct MeshStorageStatistics
{
    /** Local point count. */
    std::uint64_t points = 0;
    /** Local edge count. */
    std::uint64_t edges = 0;
    /** Local face count. */
    std::uint64_t faces = 0;
    /** Local cell count. */
    std::uint64_t cells = 0;
    /** Sum of local geometric point references across cells, faces, and edges. */
    std::uint64_t geometricPointReferences = 0;
    /** Payload bytes reserved by entity point-pointer vectors. */
    std::uint64_t pointPointerCapacityBytes = 0;
    /** Inline bytes obtained from entity counts times entity `sizeof`. */
    std::uint64_t inlineEntityBytes = 0;
    /** Raw physical-coordinate scalar bytes, excluding object overhead. */
    std::uint64_t coordinateScalarBytes = 0;
    /** Inline bytes in one stored point record. */
    std::uint64_t pointInlineBytes = 0;
    /** Inline bytes in one stored edge record. */
    std::uint64_t edgeInlineBytes = 0;
    /** Inline bytes in one stored face record. */
    std::uint64_t faceInlineBytes = 0;
    /** Inline bytes in one stored cell record. */
    std::uint64_t cellInlineBytes = 0;
    /** Populated remote-partition-ID sidecar payloads. */
    std::uint64_t remoteIdPayloads = 0;
    /** Populated incidence sidecar payloads. */
    std::uint64_t incidencePayloads = 0;
    /** Populated marker sidecar payloads. */
    std::uint64_t markerPayloads = 0;
    /** Populated cached-measure sidecar payloads. */
    std::uint64_t measurePayloads = 0;
    /** Bytes in the compact trivially copyable entity core. */
    std::uint64_t entityCoreBytes = 0;
    /** Allocated compact local-index connectivity entries. */
    std::uint64_t compactConnectivityEntries = 0;
    /** Connectivity entries belonging to live entities. */
    std::uint64_t activeCompactConnectivityEntries = 0;
    /** Released connectivity entries available for reuse. */
    std::uint64_t reusableCompactConnectivityEntries = 0;
    /** Scalars in aligned mesh-owned coordinate storage. */
    std::uint64_t compactCoordinateScalars = 0;
    /** Bytes in aligned mesh-owned coordinate storage. */
    std::uint64_t compactCoordinateBytes = 0;
    /** Entries in the sparse external-point-ID lookup. */
    std::uint64_t externalPointLookupEntries = 0;

    /**
     * @brief Return the total number of stored topological entities.
     * @return sum of point, edge, face, and cell counts
     */
    [[nodiscard]] std::uint64_t entityCount() const noexcept
    {
        return points + edges + faces + cells;
    }
};

/**
 * @brief Accumulate one entity-dimension's sparse-sidecar statistics.
 * @param source source statistics
 * @param destination aggregate benchmark statistics
 */
void accumulateSparseMetadataStatistics(
    GeoEntitySparseMetadataStatistics const& source,
    MeshStorageStatistics& destination )
{
    destination.remoteIdPayloads += source.remoteIdEntries;
    destination.incidencePayloads += source.incidenceEntries;
    destination.markerPayloads += source.markerEntries;
    destination.measurePayloads += source.measureEntries;
}

/**
 * @brief Accumulate one entity-dimension's compact connectivity counters.
 * @param source source connectivity statistics
 * @param destination aggregate benchmark statistics
 */
void accumulateConnectivityStatistics(
    GeoEntityConnectivityStatistics const& source,
    MeshStorageStatistics& destination )
{
    destination.compactConnectivityEntries +=
        source.allocatedEntries;
    destination.activeCompactConnectivityEntries +=
        source.activeEntries;
    destination.reusableCompactConnectivityEntries +=
        source.reusableEntries;
}

/**
 * @brief Return the current process resident set size on Linux.
 * @return resident bytes, or zero when unavailable
 */
[[nodiscard]] std::uint64_t residentBytes()
{
#if defined( __linux__ )
    std::ifstream status( "/proc/self/statm" );
    std::uint64_t totalPages = 0;
    std::uint64_t residentPages = 0;
    status >> totalPages >> residentPages;
    if ( !status )
        return 0;
    auto const pageSize = ::sysconf( _SC_PAGESIZE );
    return pageSize > 0
               ? residentPages * static_cast<std::uint64_t>( pageSize )
               : 0;
#else
    return 0;
#endif
}

/**
 * @brief Return bytes currently allocated in the glibc heap.
 * @return live heap bytes, or zero when unavailable
 */
[[nodiscard]] std::uint64_t liveHeapBytes()
{
#if defined( __GLIBC__ )
    return static_cast<std::uint64_t>( ::mallinfo2().uordblks );
#else
    return 0;
#endif
}

/**
 * @brief Form a nonnegative difference for process-memory counters.
 * @param after counter after an operation
 * @param before counter before an operation
 * @return `max(after-before,0)`
 */
[[nodiscard]] std::uint64_t
nonnegativeDifference( std::uint64_t after, std::uint64_t before ) noexcept
{
    return after >= before ? after - before : 0;
}

/**
 * @brief Return the geometric cell-name component used in benchmark labels.
 * @tparam MeshType Feel++ mesh type
 * @return `simplex` or `hypercube`
 */
template <typename MeshType>
[[nodiscard]] constexpr std::string_view cellName() noexcept
{
    if constexpr ( MeshType::shape_type::is_simplex )
        return "simplex";
    else
        return "hypercube";
}

/**
 * @brief Construct and fully update a benchmark mesh.
 *
 * The generated mesh file is cached by Feel++; after a warmup construction,
 * repeated calls benchmark import and topology update rather than Gmsh
 * generation.
 *
 * @tparam MeshType Feel++ mesh type
 * @param resolution inverse characteristic mesh size
 * @return fully updated distributed mesh
 */
template <typename MeshType>
[[nodiscard]] std::shared_ptr<MeshType>
makeBenchmarkMesh( int resolution )
{
    auto const meshSize = 1.0 / static_cast<double>( resolution );
    auto const name =
        "geomesh-storage-" + std::string( cellName<MeshType>() ) + "-" +
        std::to_string( MeshType::nDim ) + "d-g" +
        std::to_string( MeshType::shape_type::nOrder ) + "-r" +
        std::to_string( resolution );

    return createGMSHMesh(
        _mesh = new MeshType,
        _desc = domain(
            _name = name, _shape = std::string( cellName<MeshType>() ),
            _convex = std::string( cellName<MeshType>() ),
            _dim = MeshType::nDim,
            _order = MeshType::shape_type::nOrder, _h = meshSize ),
        _h = meshSize, _verbose = 0, _straighten = false,
        _force_rebuild = false );
}

/**
 * @brief Accumulate canonical connectivity and lazy compatibility payload.
 * @tparam EntityType geometric entity type
 * @param entity geometric entity
 * @param statistics output statistics
 */
template <typename EntityType>
void accumulateEntityPointStorage(
    EntityType const& entity, MeshStorageStatistics& statistics )
{
    statistics.geometricPointReferences +=
        entity.pointIds().size();
    statistics.pointPointerCapacityBytes +=
        entity.pointPointerCacheCapacityBytes();
}

/**
 * @brief Inspect entity counts and conservative storage lower bounds.
 * @tparam MeshType Feel++ mesh type
 * @param mesh mesh to inspect
 * @return collected storage statistics
 */
template <typename MeshType>
[[nodiscard]] MeshStorageStatistics
inspectMeshStorage( std::shared_ptr<MeshType> const& mesh )
{
    MeshStorageStatistics statistics;
    statistics.points = mesh->numPoints();
    statistics.edges = mesh->numEdges();
    statistics.faces = mesh->numFaces();
    statistics.cells = mesh->numElements();

    statistics.inlineEntityBytes =
        statistics.points * sizeof( typename MeshType::point_type ) +
        statistics.edges * sizeof( typename MeshType::edge_type ) +
        statistics.faces * sizeof( typename MeshType::face_type ) +
        statistics.cells * sizeof( typename MeshType::element_type );
    statistics.coordinateScalarBytes =
        statistics.points * MeshType::nRealDim *
        sizeof( typename MeshType::value_type );
    statistics.pointInlineBytes =
        sizeof( typename MeshType::point_type );
    statistics.edgeInlineBytes =
        sizeof( typename MeshType::edge_type );
    statistics.faceInlineBytes =
        sizeof( typename MeshType::face_type );
    statistics.cellInlineBytes =
        sizeof( typename MeshType::element_type );
    statistics.entityCoreBytes =
        sizeof( typename MeshType::element_type::core_type );
    auto const coordinateStatistics =
        mesh->pointCoordinateStorageStatistics();
    statistics.compactCoordinateScalars =
        coordinateStatistics.scalarEntries;
    statistics.compactCoordinateBytes =
        coordinateStatistics.scalarEntries *
        sizeof( typename MeshType::value_type );
    statistics.externalPointLookupEntries =
        coordinateStatistics.externalLookupEntries;

    accumulateSparseMetadataStatistics(
        mesh->pointSparseMetadataStatistics(), statistics );
    accumulateSparseMetadataStatistics(
        mesh->faceSparseMetadataStatistics(), statistics );
    accumulateSparseMetadataStatistics(
        mesh->elementSparseMetadataStatistics(), statistics );
    accumulateConnectivityStatistics(
        mesh->pointConnectivityStatistics(), statistics );
    accumulateConnectivityStatistics(
        mesh->faceConnectivityStatistics(), statistics );
    accumulateConnectivityStatistics(
        mesh->elementConnectivityStatistics(), statistics );
    if constexpr ( requires {
                       mesh->edgeSparseMetadataStatistics();
                   } )
    {
        accumulateSparseMetadataStatistics(
            mesh->edgeSparseMetadataStatistics(), statistics );
        accumulateConnectivityStatistics(
            mesh->edgeConnectivityStatistics(), statistics );
    }

    for ( auto const& [id, element] : mesh->elements() )
    {
        benchmark::DoNotOptimize( id );
        accumulateEntityPointStorage( element, statistics );
    }
    for ( auto const& [id, face] : mesh->faces() )
    {
        benchmark::DoNotOptimize( id );
        accumulateEntityPointStorage( face, statistics );
    }
    for ( auto const& [id, edge] : mesh->edges() )
    {
        benchmark::DoNotOptimize( id );
        accumulateEntityPointStorage( edge, statistics );
    }
    return statistics;
}

/**
 * @brief Publish mesh counts and static storage counters.
 * @param state Google Benchmark state
 * @param statistics collected mesh statistics
 */
void publishStorageCounters(
    benchmark::State& state,
    MeshStorageStatistics const& statistics )
{
    state.counters["points"] =
        static_cast<double>( statistics.points );
    state.counters["edges"] =
        static_cast<double>( statistics.edges );
    state.counters["faces"] =
        static_cast<double>( statistics.faces );
    state.counters["cells"] =
        static_cast<double>( statistics.cells );
    state.counters["entities"] =
        static_cast<double>( statistics.entityCount() );
    state.counters["geometric-point-references"] =
        static_cast<double>( statistics.geometricPointReferences );
    state.counters["inline-entity-bytes"] =
        static_cast<double>( statistics.inlineEntityBytes );
    state.counters["point-pointer-capacity-bytes"] =
        static_cast<double>( statistics.pointPointerCapacityBytes );
    state.counters["coordinate-scalar-bytes"] =
        static_cast<double>( statistics.coordinateScalarBytes );
    state.counters["inline-bytes-per-point"] =
        static_cast<double>( statistics.pointInlineBytes );
    state.counters["inline-bytes-per-edge"] =
        static_cast<double>( statistics.edgeInlineBytes );
    state.counters["inline-bytes-per-face"] =
        static_cast<double>( statistics.faceInlineBytes );
    state.counters["inline-bytes-per-cell"] =
        static_cast<double>( statistics.cellInlineBytes );
    state.counters["compact-entity-core-bytes"] =
        static_cast<double>( statistics.entityCoreBytes );
    state.counters["remote-id-sidecar-payloads"] =
        static_cast<double>( statistics.remoteIdPayloads );
    state.counters["incidence-sidecar-payloads"] =
        static_cast<double>( statistics.incidencePayloads );
    state.counters["marker-sidecar-payloads"] =
        static_cast<double>( statistics.markerPayloads );
    state.counters["measure-sidecar-payloads"] =
        static_cast<double>( statistics.measurePayloads );
    state.counters["compact-connectivity-entries"] =
        static_cast<double>(
            statistics.compactConnectivityEntries );
    state.counters["active-compact-connectivity-entries"] =
        static_cast<double>(
            statistics.activeCompactConnectivityEntries );
    state.counters["reusable-compact-connectivity-entries"] =
        static_cast<double>(
            statistics.reusableCompactConnectivityEntries );
    state.counters["compact-connectivity-bytes"] =
        static_cast<double>(
            statistics.compactConnectivityEntries *
            sizeof( std::uint32_t ) );
    state.counters["compact-coordinate-scalars"] =
        static_cast<double>(
            statistics.compactCoordinateScalars );
    state.counters["compact-coordinate-bytes"] =
        static_cast<double>(
            statistics.compactCoordinateBytes );
    state.counters["external-point-lookup-entries"] =
        static_cast<double>(
            statistics.externalPointLookupEntries );
}

/**
 * @brief Measure retained footprint and C++ allocation activity.
 * @tparam MeshType Feel++ mesh type
 * @param state Google Benchmark state; first argument is mesh resolution
 */
template <typename MeshType>
void benchmarkMeshFootprint( benchmark::State& state )
{
    auto const heapBefore = liveHeapBytes();
    auto const rssBefore = residentBytes();
    geomesh_allocation_probe::begin();
    auto mesh = makeBenchmarkMesh<MeshType>(
        static_cast<int>( state.range( 0 ) ) );
    auto const allocations = geomesh_allocation_probe::end();
    auto const heapAfter = liveHeapBytes();
    auto const rssAfter = residentBytes();
    auto const statistics = inspectMeshStorage( mesh );

    for ( [[maybe_unused]] auto iteration : state )
        benchmark::DoNotOptimize( mesh.get() );

    publishStorageCounters( state, statistics );
    state.counters["cpp-allocation-calls"] =
        static_cast<double>( allocations.calls );
    state.counters["cpp-requested-bytes"] =
        static_cast<double>( allocations.bytes );
    state.counters["live-heap-delta-bytes"] =
        static_cast<double>(
            nonnegativeDifference( heapAfter, heapBefore ) );
    state.counters["rss-delta-bytes"] =
        static_cast<double>(
            nonnegativeDifference( rssAfter, rssBefore ) );

    auto const entityCount =
        std::max<std::uint64_t>( statistics.entityCount(), 1 );
    auto const pointCount =
        std::max<std::uint64_t>( statistics.points, 1 );
    state.counters["live-heap-bytes-per-entity"] =
        static_cast<double>(
            nonnegativeDifference( heapAfter, heapBefore ) ) /
        static_cast<double>( entityCount );
    state.counters["live-heap-bytes-per-geometric-node"] =
        static_cast<double>(
            nonnegativeDifference( heapAfter, heapBefore ) ) /
        static_cast<double>( pointCount );
    state.SetItemsProcessed(
        state.iterations() * statistics.entityCount() );
}

/**
 * @brief Benchmark cached mesh import and topology update.
 * @tparam MeshType Feel++ mesh type
 * @param state Google Benchmark state; first argument is mesh resolution
 */
template <typename MeshType>
void benchmarkMeshImportAndUpdate( benchmark::State& state )
{
    auto warmup = makeBenchmarkMesh<MeshType>(
        static_cast<int>( state.range( 0 ) ) );
    auto const warmupStatistics = inspectMeshStorage( warmup );
    warmup.reset();

    geomesh_allocation_probe::begin();
    for ( [[maybe_unused]] auto iteration : state )
    {
        auto mesh = makeBenchmarkMesh<MeshType>(
            static_cast<int>( state.range( 0 ) ) );
        benchmark::DoNotOptimize( mesh.get() );
    }
    auto const allocations = geomesh_allocation_probe::end();

    publishStorageCounters( state, warmupStatistics );
    state.counters["cpp-allocation-calls"] =
        benchmark::Counter(
            static_cast<double>( allocations.calls ),
            benchmark::Counter::kAvgIterations );
    state.counters["cpp-requested-bytes"] =
        benchmark::Counter(
            static_cast<double>( allocations.bytes ),
            benchmark::Counter::kAvgIterations );
    state.SetItemsProcessed(
        state.iterations() * warmupStatistics.entityCount() );
}

/**
 * @brief Benchmark traversal of cell topology without coordinate gathering.
 * @tparam MeshType Feel++ mesh type
 * @param state Google Benchmark state; first argument is mesh resolution
 */
template <typename MeshType>
void benchmarkElementTraversal( benchmark::State& state )
{
    auto mesh = makeBenchmarkMesh<MeshType>(
        static_cast<int>( state.range( 0 ) ) );
    std::uint64_t checksum = 0;

    for ( [[maybe_unused]] auto iteration : state )
    {
        for ( auto const& [id, element] : mesh->elements() )
        {
            checksum += id;
            checksum += element.processId();
            checksum += element.nPoints();
            checksum += element.nNeighbors();
        }
        benchmark::DoNotOptimize( checksum );
    }
    state.SetItemsProcessed(
        state.iterations() * mesh->numElements() );
}

/**
 * @brief Benchmark traversal of stored facets and their cell connectivity.
 * @tparam MeshType Feel++ mesh type
 * @param state Google Benchmark state; first argument is mesh resolution
 */
template <typename MeshType>
void benchmarkFacetTraversal( benchmark::State& state )
{
    auto mesh = makeBenchmarkMesh<MeshType>(
        static_cast<int>( state.range( 0 ) ) );
    std::uint64_t checksum = 0;

    for ( [[maybe_unused]] auto iteration : state )
    {
        for ( auto const& [id, face] : mesh->faces() )
        {
            checksum += id;
            checksum += face.processId();
            checksum += face.nPoints();
            checksum += face.isConnectedTo0();
            checksum += face.isConnectedTo1();
        }
        benchmark::DoNotOptimize( checksum );
    }
    state.SetItemsProcessed(
        state.iterations() * mesh->numFaces() );
}

/**
 * @brief Benchmark legacy allocation-based `element.G()` coordinate gathering.
 * @tparam MeshType Feel++ mesh type
 * @param state Google Benchmark state; first argument is mesh resolution
 */
template <typename MeshType>
void benchmarkCoordinateGather( benchmark::State& state )
{
    auto mesh = makeBenchmarkMesh<MeshType>(
        static_cast<int>( state.range( 0 ) ) );

    geomesh_allocation_probe::begin();
    for ( [[maybe_unused]] auto iteration : state )
    {
        for ( auto const& [id, element] : mesh->elements() )
        {
            benchmark::DoNotOptimize( id );
            auto coordinates = element.G();
            benchmark::DoNotOptimize( coordinates.data().begin() );
        }
    }
    auto const allocations = geomesh_allocation_probe::end();

    state.counters["cpp-allocation-calls"] =
        benchmark::Counter(
            static_cast<double>( allocations.calls ),
            benchmark::Counter::kAvgIterations );
    state.counters["cpp-requested-bytes"] =
        benchmark::Counter(
            static_cast<double>( allocations.bytes ),
            benchmark::Counter::kAvgIterations );
    state.SetItemsProcessed(
        state.iterations() * mesh->numElements() );
}

/**
 * @brief Benchmark compact-ID gathering from the aligned Eigen coordinate map.
 * @tparam MeshType Feel++ mesh type
 * @param state Google Benchmark state; first argument is mesh resolution
 */
template <typename MeshType>
void benchmarkCompactCoordinateGather( benchmark::State& state )
{
    auto mesh = makeBenchmarkMesh<MeshType>(
        static_cast<int>( state.range( 0 ) ) );
    auto const& constMesh = *mesh;
    auto coordinates = constMesh.pointCoordinates();
    Eigen::Matrix<
        typename MeshType::value_type, MeshType::nRealDim,
        Eigen::Dynamic, Eigen::ColMajor>
        localCoordinates(
            MeshType::nRealDim,
            MeshType::element_type::numPoints );

    geomesh_allocation_probe::begin();
    for ( [[maybe_unused]] auto iteration : state )
    {
        for ( auto const& [id, element] : mesh->elements() )
        {
            benchmark::DoNotOptimize( id );
            auto const pointIds = element.pointIds();
            for ( uint16_type localPoint = 0;
                  localPoint < pointIds.size(); ++localPoint )
                localCoordinates.col( localPoint ) =
                    coordinates.col(
                        mesh->pointLocalIndex(
                            pointIds[localPoint] ) );
            benchmark::DoNotOptimize( localCoordinates.data() );
        }
    }
    auto const allocations = geomesh_allocation_probe::end();

    state.counters["cpp-allocation-calls"] =
        benchmark::Counter(
            static_cast<double>( allocations.calls ),
            benchmark::Counter::kAvgIterations );
    state.counters["cpp-requested-bytes"] =
        benchmark::Counter(
            static_cast<double>( allocations.bytes ),
            benchmark::Counter::kAvgIterations );
    state.SetItemsProcessed(
        state.iterations() * mesh->numElements() );
}

/**
 * @brief Benchmark compact gathering through mathematical entity views.
 * @tparam MeshType Feel++ mesh type
 * @param state Google Benchmark state; first argument is mesh resolution
 */
template <typename MeshType>
void benchmarkEntityViewCoordinateGather( benchmark::State& state )
{
    auto mesh = makeBenchmarkMesh<MeshType>(
        static_cast<int>( state.range( 0 ) ) );
    Eigen::Matrix<
        typename MeshType::value_type, MeshType::nRealDim,
        Eigen::Dynamic, Eigen::ColMajor>
        localCoordinates(
            MeshType::nRealDim,
            MeshType::element_type::numPoints );

    geomesh_allocation_probe::begin();
    for ( [[maybe_unused]] auto iteration : state )
    {
        for ( auto const cell : mesh->elementViews() )
        {
            cell.updateCoordinates( localCoordinates );
            benchmark::DoNotOptimize(
                localCoordinates.data() );
        }
    }
    auto const allocations = geomesh_allocation_probe::end();

    state.counters["cpp-allocation-calls"] =
        benchmark::Counter(
            static_cast<double>( allocations.calls ),
            benchmark::Counter::kAvgIterations );
    state.counters["cpp-requested-bytes"] =
        benchmark::Counter(
            static_cast<double>( allocations.bytes ),
            benchmark::Counter::kAvgIterations );
    state.SetItemsProcessed(
        state.iterations() * mesh->numElements() );
}

/**
 * @brief Benchmark binary serialization of the current owning mesh graph.
 * @tparam MeshType Feel++ mesh type
 * @param state Google Benchmark state; first argument is mesh resolution
 */
template <typename MeshType>
void benchmarkMeshSerialization( benchmark::State& state )
{
    auto mesh = makeBenchmarkMesh<MeshType>(
        static_cast<int>( state.range( 0 ) ) );
    std::size_t serializedBytes = 0;

    for ( [[maybe_unused]] auto iteration : state )
    {
        std::ostringstream stream(
            std::ios::out | std::ios::binary );
        {
            boost::archive::binary_oarchive archive( stream );
            archive << *mesh;
        }
        serializedBytes = stream.str().size();
        benchmark::DoNotOptimize( serializedBytes );
    }
    state.counters["serialized-bytes"] =
        static_cast<double>( serializedBytes );
    state.SetBytesProcessed(
        state.iterations() * serializedBytes );
}

#if defined( FEELPP_HAS_METIS )
/**
 * @brief Benchmark serial Metis partition metadata construction.
 * @tparam MeshType Feel++ mesh type
 * @param state Google Benchmark state; arguments are resolution and partitions
 */
template <typename MeshType>
void benchmarkMeshPartitioning( benchmark::State& state )
{
    if ( Environment::worldComm().globalSize() != 1 )
    {
        state.SkipWithError(
            "partition baseline must run in a single process" );
        return;
    }

    std::uint64_t elementCount = 0;
    for ( [[maybe_unused]] auto iteration : state )
    {
        state.PauseTiming();
        auto mesh = makeBenchmarkMesh<MeshType>(
            static_cast<int>( state.range( 0 ) ) );
        elementCount = mesh->numElements();
        state.ResumeTiming();

        auto partitions = partitionMesh(
            mesh, static_cast<rank_type>( state.range( 1 ) ) );
        benchmark::DoNotOptimize( partitions.get() );
    }
    state.SetItemsProcessed(
        state.iterations() * elementCount );
}
#endif

/**
 * @brief Register every storage operation for one geometric mesh type.
 * @tparam MeshType Feel++ mesh type
 * @param name stable benchmark name prefix
 * @param resolutions mesh-resolution arguments
 */
template <typename MeshType>
void registerMeshStorageBenchmarks(
    std::string const& name,
    std::initializer_list<std::int64_t> resolutions )
{
    for ( auto const resolution : resolutions )
    {
        benchmark::RegisterBenchmark(
            ( name + "/footprint" ).c_str(),
            &benchmarkMeshFootprint<MeshType> )
            ->Arg( resolution )
            ->Iterations( 1 );
        benchmark::RegisterBenchmark(
            ( name + "/import-update" ).c_str(),
            &benchmarkMeshImportAndUpdate<MeshType> )
            ->Arg( resolution )
            ->Unit( benchmark::kMillisecond )
            ->Iterations( 3 );
        benchmark::RegisterBenchmark(
            ( name + "/element-traversal" ).c_str(),
            &benchmarkElementTraversal<MeshType> )
            ->Arg( resolution )
            ->Unit( benchmark::kMicrosecond );
        benchmark::RegisterBenchmark(
            ( name + "/facet-traversal" ).c_str(),
            &benchmarkFacetTraversal<MeshType> )
            ->Arg( resolution )
            ->Unit( benchmark::kMicrosecond );
        benchmark::RegisterBenchmark(
            ( name + "/coordinate-gather" ).c_str(),
            &benchmarkCoordinateGather<MeshType> )
            ->Arg( resolution )
            ->Unit( benchmark::kMicrosecond );
        benchmark::RegisterBenchmark(
            ( name + "/compact-coordinate-gather" ).c_str(),
            &benchmarkCompactCoordinateGather<MeshType> )
            ->Arg( resolution )
            ->Unit( benchmark::kMicrosecond );
        benchmark::RegisterBenchmark(
            ( name + "/entity-view-coordinate-gather" ).c_str(),
            &benchmarkEntityViewCoordinateGather<MeshType> )
            ->Arg( resolution )
            ->Unit( benchmark::kMicrosecond );
        benchmark::RegisterBenchmark(
            ( name + "/serialization" ).c_str(),
            &benchmarkMeshSerialization<MeshType> )
            ->Arg( resolution )
            ->Unit( benchmark::kMillisecond );
#if defined( FEELPP_HAS_METIS )
        benchmark::RegisterBenchmark(
            ( name + "/partition" ).c_str(),
            &benchmarkMeshPartitioning<MeshType> )
            ->Args( { resolution, 2 } )
            ->Unit( benchmark::kMillisecond )
            ->Iterations( 3 );
#endif
    }
}

/**
 * @brief Register the complete static mesh-storage benchmark matrix.
 * @return true after registration
 */
[[nodiscard]] bool registerStorageMatrix()
{
    registerMeshStorageBenchmarks<Mesh<Simplex<2, 1>>>(
        "geomesh/simplex/2d/g1", { 4, 8 } );
    registerMeshStorageBenchmarks<Mesh<Simplex<2, 2>>>(
        "geomesh/simplex/2d/g2", { 4, 8 } );
    registerMeshStorageBenchmarks<Mesh<Simplex<3, 1>>>(
        "geomesh/simplex/3d/g1", { 2, 4 } );
    registerMeshStorageBenchmarks<Mesh<Simplex<3, 2>>>(
        "geomesh/simplex/3d/g2", { 2, 4 } );

    registerMeshStorageBenchmarks<Mesh<Hypercube<2, 1>>>(
        "geomesh/hypercube/2d/g1", { 4, 8 } );
    registerMeshStorageBenchmarks<Mesh<Hypercube<2, 2>>>(
        "geomesh/hypercube/2d/g2", { 4, 8 } );
    registerMeshStorageBenchmarks<Mesh<Hypercube<3, 1>>>(
        "geomesh/hypercube/3d/g1", { 2, 4 } );
    registerMeshStorageBenchmarks<Mesh<Hypercube<3, 2>>>(
        "geomesh/hypercube/3d/g2", { 2, 4 } );

#if FEELPP_MESH_MAX_ORDER >= 3
    registerMeshStorageBenchmarks<Mesh<Simplex<2, 3>>>(
        "geomesh/simplex/2d/g3", { 4, 8 } );
    registerMeshStorageBenchmarks<Mesh<Simplex<3, 3>>>(
        "geomesh/simplex/3d/g3", { 2, 4 } );
#endif
    return true;
}

/** Trigger benchmark registration before `main()`. */
[[maybe_unused]] bool const storageMatrixRegistered =
    registerStorageMatrix();
} // namespace

/**
 * @brief Initialize Feel++ and execute the geometric mesh storage benchmarks.
 * @param argc command-line argument count
 * @param argv command-line argument vector
 * @return process exit code
 */
int main( int argc, char** argv )
{
    Environment environment(
        _argc = argc, _argv = argv,
        _about = about(
            _name = "feelpp_bench_geomesh_storage",
            _author = "Christophe Prud'homme",
            _email = "christophe.prudhomme@feelpp.org" ) );

    benchmark::Initialize( &argc, argv );
    benchmark::RunSpecifiedBenchmarks();
}
