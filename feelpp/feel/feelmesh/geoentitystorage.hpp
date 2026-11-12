/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file geoentitystorage.hpp
 * @brief Compact geometric-entity records and mesh-owned sparse metadata.
 */
#if !defined( FEELPP_FEELMESH_GEOENTITYSTORAGE_HPP )
#define FEELPP_FEELMESH_GEOENTITYSTORAGE_HPP 1

#include <feel/feelcore/feeltypes.hpp>

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Feel
{

/**
 * @struct GeoEntityCore
 * @brief Compact, trivially copyable hot data shared by geometric entities.
 *
 * The record contains identity, packed topology/context flags, ownership, and
 * the active connectivity span in mesh-owned contiguous index storage.
 * Legacy pointer access can coexist temporarily as a compatibility cache
 * while consumers migrate to index-based views.
 *
 * @tparam IndexT external entity identifier type
 */
template <std::integral IndexT>
struct GeoEntityCore
{
    /** Entity identifier. */
    IndexT id = 0;
    /** Packed geometric context, reference geometry, shape, and flags. */
    uint32_type entityContext = 0;
    /** Offset in mesh-owned connectivity storage. */
    uint32_type connectivityOffset = std::numeric_limits<uint32_type>::max();
    /** Owning process identifier. */
    rank_type processId = invalid_rank_type_value;
    /** Process identifier in the current partition representation. */
    rank_type processIdInPartition = invalid_rank_type_value;
    /** Highest-dimensional boundary entity touching this entity. */
    uint16_type boundaryEntityDimension = invalid_uint16_type_value;
    /** Number of entries in the mesh-owned connectivity span. */
    uint16_type connectivityCount = 0;
};

static_assert( std::is_standard_layout_v<GeoEntityCore<uint32_type>> );
static_assert( std::is_trivially_copyable_v<GeoEntityCore<uint32_type>> );

/**
 * @struct GeoEntityMeasureMetadata
 * @brief Lazily allocated cached measures for one geometric entity.
 *
 * @tparam T coordinate and measure scalar type
 */
template <typename T>
struct GeoEntityMeasureMetadata
{
    /** Cached entity measure. */
    std::optional<T> measure;
    /** Cached measures of topological faces. */
    std::vector<std::optional<T>> faceMeasures;
    /** Cached measure of the point-element neighbor patch. */
    std::optional<T> pointElementNeighborsMeasure;

    /** @brief Return true when no cached measure is present. */
    [[nodiscard]] bool empty() const noexcept
    {
        return !measure && faceMeasures.empty() &&
               !pointElementNeighborsMeasure;
    }
};

/**
 * @struct GeoEntityDetachedMetadata
 * @brief Sparse metadata temporarily owned by a detached importer entity.
 *
 * Stored mesh entities transfer this payload to a mesh-owned
 * `GeoEntityStorage`. The detached form exists so importers and
 * local construction can keep setting markers before insertion.
 *
 * @tparam IndexT entity identifier type
 * @tparam MarkerT Feel++ marker type
 * @tparam T coordinate and measure scalar type
 */
template <std::integral IndexT, typename MarkerT, typename T>
struct GeoEntityDetachedMetadata
{
    /** Local connectivity represented by external entity identifiers. */
    std::vector<IndexT> connectivity;
    /** Remote partition identifier mapping. */
    std::map<rank_type, IndexT> remoteIds;
    /** Incident element identifiers and local positions. */
    std::set<std::pair<IndexT, uint16_type>> elementIncidence;
    /** Marker values indexed by marker kind. */
    std::map<uint16_type, MarkerT> markers;
    /** Optional cached measures. */
    GeoEntityMeasureMetadata<T> measures;

    /** @brief Return true when every sparse payload is empty. */
    [[nodiscard]] bool empty() const noexcept
    {
        return connectivity.empty() && remoteIds.empty() &&
               elementIncidence.empty() && markers.empty() &&
               measures.empty();
    }
};

/**
 * @struct GeoEntitySparseMetadataStatistics
 * @brief Number of populated mesh-owned sparse metadata payloads.
 */
struct GeoEntitySparseMetadataStatistics
{
    /** Entities carrying remote-partition identifiers. */
    size_type remoteIdEntries = 0;
    /** Entities carrying point/edge-to-element incidence. */
    size_type incidenceEntries = 0;
    /** Entities carrying markers. */
    size_type markerEntries = 0;
    /** Entities carrying cached measures. */
    size_type measureEntries = 0;

    /** @brief Return the number of independent populated sidecar payloads. */
    [[nodiscard]] size_type payloadCount() const noexcept
    {
        return remoteIdEntries + incidenceEntries + markerEntries +
               measureEntries;
    }
};

/**
 * @struct GeoEntityConnectivityStatistics
 * @brief Occupancy counters for contiguous local connectivity storage.
 */
struct GeoEntityConnectivityStatistics
{
    /** Total allocated connectivity slots, including released spans. */
    size_type allocatedEntries = 0;
    /** Connectivity slots belonging to live entities. */
    size_type activeEntries = 0;
    /** Released slots available for exact-cardinality reuse. */
    size_type reusableEntries = 0;
};

/**
 * @class GeoEntityStorage
 * @brief Mesh-owned contiguous connectivity and independent sparse sidecars.
 *
 * Connectivity is stored in one contiguous local-index array. Each sparse
 * metadata category has its own associative container. Consequently an
 * entity with only a marker does not also pay for empty remote-ID and
 * incidence trees. Released connectivity spans are reused by exact
 * cardinality, preserving stable offsets for every other live entity.
 *
 * @tparam IndexT entity identifier type
 * @tparam MarkerT Feel++ marker type
 * @tparam T coordinate and measure scalar type
 */
template <std::integral IndexT, typename MarkerT, typename T>
class GeoEntityStorage
{
  public:
    /** Entity identifier type. */
    using index_type = IndexT;
    /** Marker type. */
    using marker_type = MarkerT;
    /** Scalar value type. */
    using value_type = T;
    /** Remote partition identifier mapping. */
    using remote_id_map_type = std::map<rank_type, index_type>;
    /** Point/edge-to-element incidence set. */
    using incidence_set_type =
        std::set<std::pair<index_type, uint16_type>>;
    /** Marker mapping. */
    using marker_map_type = std::map<uint16_type, marker_type>;
    /** Cached measure payload. */
    using measure_metadata_type = GeoEntityMeasureMetadata<value_type>;
    /** Detached aggregate used during construction and copying. */
    using detached_metadata_type =
        GeoEntityDetachedMetadata<index_type, marker_type, value_type>;

    GeoEntityStorage() = default;
    GeoEntityStorage( GeoEntityStorage const& ) = delete;
    GeoEntityStorage& operator=( GeoEntityStorage const& ) = delete;
    GeoEntityStorage( GeoEntityStorage&& ) noexcept = default;
    GeoEntityStorage& operator=( GeoEntityStorage&& ) noexcept = default;

    /** @brief Remove all sparse metadata payloads. */
    void clear() noexcept
    {
        M_connectivity.clear();
        M_reusableConnectivity.clear();
        M_activeConnectivityEntries = 0;
        M_remoteIds.clear();
        M_incidence.clear();
        M_markers.clear();
        M_measures.clear();
    }

    /**
     * @brief Reserve contiguous connectivity for fixed-cardinality entities.
     * @param entityCount expected entity count
     * @param cardinality connectivity entries per entity
     */
    void reserveConnectivity(
        size_type entityCount, uint16_type cardinality )
    {
        auto const entries =
            static_cast<std::uint64_t>( entityCount ) * cardinality;
        if ( entries > std::numeric_limits<uint32_type>::max() )
            throw std::length_error(
                "geometric connectivity exceeds 32-bit offsets" );
        M_connectivity.reserve( static_cast<size_type>( entries ) );
    }

    /**
     * @brief Store one entity's contiguous connectivity.
     * @param connectivity external identifiers in local mathematical order
     * @return offset of the stored span
     */
    [[nodiscard]] uint32_type allocateConnectivity(
        std::span<index_type const> connectivity )
    {
        if ( connectivity.empty() )
            return std::numeric_limits<uint32_type>::max();
        if ( connectivity.size() >
             std::numeric_limits<uint16_type>::max() )
            throw std::length_error(
                "entity connectivity exceeds 16-bit cardinality" );

        auto const cardinality =
            static_cast<uint16_type>( connectivity.size() );
        auto reusable = std::find_if(
            M_reusableConnectivity.begin(),
            M_reusableConnectivity.end(),
            [cardinality]( auto const& span )
            {
                return span.second == cardinality;
            } );
        uint32_type offset = 0;
        if ( reusable != M_reusableConnectivity.end() )
        {
            offset = reusable->first;
            std::copy(
                connectivity.begin(), connectivity.end(),
                M_connectivity.begin() + offset );
            M_reusableConnectivity.erase( reusable );
        }
        else
        {
            if ( M_connectivity.size() + connectivity.size() >
                 std::numeric_limits<uint32_type>::max() )
                throw std::length_error(
                    "geometric connectivity exceeds 32-bit offsets" );
            offset = static_cast<uint32_type>( M_connectivity.size() );
            M_connectivity.insert(
                M_connectivity.end(), connectivity.begin(),
                connectivity.end() );
        }
        M_activeConnectivityEntries += cardinality;
        return offset;
    }

    /**
     * @brief Release an entity connectivity span for later exact-size reuse.
     * @param offset first connectivity entry
     * @param cardinality number of entries
     */
    void releaseConnectivity(
        uint32_type offset, uint16_type cardinality )
    {
        if ( offset == std::numeric_limits<uint32_type>::max() ||
             cardinality == 0 )
            return;
        if ( static_cast<size_type>( offset ) + cardinality >
             M_connectivity.size() )
            throw std::out_of_range(
                "invalid geometric connectivity release span" );
        if ( cardinality > M_activeConnectivityEntries )
            throw std::logic_error(
                "geometric connectivity occupancy underflow" );
        if ( std::any_of(
                 M_reusableConnectivity.begin(),
                 M_reusableConnectivity.end(),
                 [offset, cardinality]( auto const& span )
                 {
                     return span.first == offset &&
                            span.second == cardinality;
                 } ) )
            throw std::logic_error(
                "geometric connectivity span released twice" );
        M_reusableConnectivity.emplace_back( offset, cardinality );
        M_activeConnectivityEntries -= cardinality;
    }

    /**
     * @brief Return a read-only contiguous connectivity view.
     * @param offset first connectivity entry
     * @param cardinality number of entries
     */
    [[nodiscard]] std::span<index_type const> connectivity(
        uint32_type offset, uint16_type cardinality ) const
    {
        if ( cardinality == 0 )
            return {};
        if ( offset == std::numeric_limits<uint32_type>::max() ||
             static_cast<size_type>( offset ) + cardinality >
                 M_connectivity.size() )
            throw std::out_of_range(
                "invalid geometric connectivity span" );
        return { M_connectivity.data() + offset, cardinality };
    }

    /**
     * @brief Update one entry of a live connectivity span.
     * @param offset first connectivity entry
     * @param cardinality number of entries
     * @param localIndex local connectivity index
     * @param id new external identifier
     */
    void setConnectivityEntry(
        uint32_type offset, uint16_type cardinality,
        uint16_type localIndex, index_type id )
    {
        if ( localIndex >= cardinality ||
             offset == std::numeric_limits<uint32_type>::max() ||
             static_cast<size_type>( offset ) + localIndex >=
                 M_connectivity.size() )
            throw std::out_of_range(
                "invalid local geometric connectivity index" );
        M_connectivity[offset + localIndex] = id;
    }

    /** @brief Return connectivity allocation and occupancy counters. */
    [[nodiscard]] GeoEntityConnectivityStatistics
    connectivityStatistics() const noexcept
    {
        size_type reusableEntries = 0;
        for ( auto const& [offset, cardinality] :
              M_reusableConnectivity )
        {
            (void)offset;
            reusableEntries += cardinality;
        }
        return {
            M_connectivity.size(), M_activeConnectivityEntries,
            reusableEntries
        };
    }

    /**
     * @brief Remove every sparse payload associated with an entity.
     * @param id entity identifier
     */
    void erase( index_type id ) noexcept
    {
        M_remoteIds.erase( id );
        M_incidence.erase( id );
        M_markers.erase( id );
        M_measures.erase( id );
    }

    /**
     * @brief Move metadata to a new entity identifier.
     * @param oldId current identifier
     * @param newId new identifier
     */
    void rekey( index_type oldId, index_type newId )
    {
        if ( oldId == newId )
            return;
        auto data = snapshot( oldId );
        erase( oldId );
        assign( newId, std::move( data ) );
    }

    /**
     * @brief Return a copy of every payload associated with an entity.
     * @param id entity identifier
     */
    [[nodiscard]] detached_metadata_type snapshot( index_type id ) const
    {
        detached_metadata_type result;
        if ( auto const it = M_remoteIds.find( id );
             it != M_remoteIds.end() )
            result.remoteIds = it->second;
        if ( auto const it = M_incidence.find( id );
             it != M_incidence.end() )
            result.elementIncidence = it->second;
        if ( auto const it = M_markers.find( id );
             it != M_markers.end() )
            result.markers = it->second;
        if ( auto const it = M_measures.find( id );
             it != M_measures.end() )
            result.measures = it->second;
        return result;
    }

    /**
     * @brief Replace an entity's metadata from a detached payload.
     * @param id entity identifier
     * @param data detached metadata payload
     */
    void assign( index_type id, detached_metadata_type data )
    {
        erase( id );
        if ( !data.remoteIds.empty() )
            M_remoteIds.emplace( id, std::move( data.remoteIds ) );
        if ( !data.elementIncidence.empty() )
            M_incidence.emplace(
                id, std::move( data.elementIncidence ) );
        if ( !data.markers.empty() )
            M_markers.emplace( id, std::move( data.markers ) );
        if ( !data.measures.empty() )
            M_measures.emplace( id, std::move( data.measures ) );
    }

    /** @brief Return remote IDs or a shared immutable empty mapping. */
    [[nodiscard]] remote_id_map_type const&
    remoteIds( index_type id ) const noexcept
    {
        auto const it = M_remoteIds.find( id );
        return it != M_remoteIds.end() ? it->second : emptyRemoteIds();
    }

    /** @brief Return mutable remote IDs, creating only this sidecar payload. */
    [[nodiscard]] remote_id_map_type& mutableRemoteIds( index_type id )
    {
        return M_remoteIds[id];
    }

    /** @brief Remove an entity's remote-ID payload. */
    void clearRemoteIds( index_type id ) noexcept
    {
        M_remoteIds.erase( id );
    }

    /** @brief Return incidence or a shared immutable empty set. */
    [[nodiscard]] incidence_set_type const&
    incidence( index_type id ) const noexcept
    {
        auto const it = M_incidence.find( id );
        return it != M_incidence.end() ? it->second : emptyIncidence();
    }

    /** @brief Return mutable incidence, creating only this sidecar payload. */
    [[nodiscard]] incidence_set_type& mutableIncidence( index_type id )
    {
        return M_incidence[id];
    }

    /** @brief Remove an entity's incidence payload. */
    void clearIncidence( index_type id ) noexcept
    {
        M_incidence.erase( id );
    }

    /** @brief Return markers or a shared immutable empty mapping. */
    [[nodiscard]] marker_map_type const&
    markers( index_type id ) const noexcept
    {
        auto const it = M_markers.find( id );
        return it != M_markers.end() ? it->second : emptyMarkers();
    }

    /** @brief Return mutable markers, creating only this sidecar payload. */
    [[nodiscard]] marker_map_type& mutableMarkers( index_type id )
    {
        return M_markers[id];
    }

    /** @brief Remove an entity's marker payload. */
    void clearMarkers( index_type id ) noexcept
    {
        M_markers.erase( id );
    }

    /** @brief Return measures or a shared immutable empty payload. */
    [[nodiscard]] measure_metadata_type const&
    measures( index_type id ) const noexcept
    {
        auto const it = M_measures.find( id );
        return it != M_measures.end() ? it->second : emptyMeasures();
    }

    /** @brief Return mutable measures, creating only this sidecar payload. */
    [[nodiscard]] measure_metadata_type& mutableMeasures( index_type id )
    {
        return M_measures[id];
    }

    /** @brief Remove an entity's cached-measure payload. */
    void clearMeasures( index_type id ) noexcept
    {
        M_measures.erase( id );
    }

    /** @brief Return populated sidecar payload counts. */
    [[nodiscard]] GeoEntitySparseMetadataStatistics statistics() const noexcept
    {
        return {
            M_remoteIds.size(), M_incidence.size(), M_markers.size(),
            M_measures.size()
        };
    }

  private:
    [[nodiscard]] static remote_id_map_type const& emptyRemoteIds() noexcept
    {
        static remote_id_map_type const value;
        return value;
    }

    [[nodiscard]] static incidence_set_type const& emptyIncidence() noexcept
    {
        static incidence_set_type const value;
        return value;
    }

    [[nodiscard]] static marker_map_type const& emptyMarkers() noexcept
    {
        static marker_map_type const value;
        return value;
    }

    [[nodiscard]] static measure_metadata_type const& emptyMeasures() noexcept
    {
        static measure_metadata_type const value;
        return value;
    }

    std::unordered_map<index_type, remote_id_map_type> M_remoteIds;
    std::unordered_map<index_type, incidence_set_type> M_incidence;
    std::unordered_map<index_type, marker_map_type> M_markers;
    std::unordered_map<index_type, measure_metadata_type> M_measures;
    std::vector<index_type> M_connectivity;
    std::vector<std::pair<uint32_type, uint16_type>>
        M_reusableConnectivity;
    size_type M_activeConnectivityEntries = 0;
};

/**
 * @brief Compatibility name for code referring specifically to sidecars.
 *
 * Connectivity and sparse metadata now share one mesh-owned entity storage.
 */
template <std::integral IndexT, typename MarkerT, typename T>
using GeoEntitySparseMetadataStorage =
    GeoEntityStorage<IndexT, MarkerT, T>;

} // namespace Feel

#endif // FEELPP_FEELMESH_GEOENTITYSTORAGE_HPP
