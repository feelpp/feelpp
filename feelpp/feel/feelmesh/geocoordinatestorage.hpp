/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file geocoordinatestorage.hpp
 * @brief Aligned Eigen coordinate ownership with external-ID lookup.
 */
#if !defined( FEELPP_FEELMESH_GEOCOORDINATESTORAGE_HPP )
#define FEELPP_FEELMESH_GEOCOORDINATESTORAGE_HPP 1

#include <feel/feelcore/feeltypes.hpp>

#include <Eigen/Core>

#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace Feel
{

/**
 * @struct GeoCoordinateStorageStatistics
 * @brief Occupancy and lookup counters for mesh-owned coordinates.
 */
struct GeoCoordinateStorageStatistics
{
    /** Number of live coordinate columns. */
    size_type points = 0;
    /** Number of stored coordinate scalars. */
    size_type scalarEntries = 0;
    /** External-ID lookup entries; zero for dense identity IDs. */
    size_type externalLookupEntries = 0;
};

/**
 * @class GeoCoordinateStorage
 * @brief Mesh-owned aligned structure of point-coordinate columns.
 *
 * Each column is one geometric point, so a cell's point identifiers can be
 * gathered directly into Eigen kernels. External identifiers equal to dense
 * local indices require no hash lookup. Imported sparse identifiers activate
 * one external-ID-to-local-index table.
 *
 * Maps and column expressions are invalidated by insertion, erase, clear, or
 * rebuild. Coordinate value updates do not invalidate them.
 *
 * @tparam RealDim ambient coordinate dimension
 * @tparam T coordinate scalar type
 * @tparam IndexT external point identifier type
 */
template <int RealDim, typename T, std::integral IndexT>
class GeoCoordinateStorage
{
    static_assert( RealDim > 0 );

  public:
    /** Coordinate scalar type. */
    using value_type = T;
    /** External point identifier type. */
    using index_type = IndexT;
    /** Dense local point index type. */
    using local_index_type = uint32_type;
    /**
     * Aligned point-column coordinate matrix type.
     *
     * Eigen requires a fixed one-row matrix with dynamic columns to be
     * row-major. Its scalar layout is identical to column-major storage in
     * one dimension because each point column contains one value.
     */
    using matrix_type =
        Eigen::Matrix<value_type, RealDim, Eigen::Dynamic,
                      RealDim == 1 ? Eigen::RowMajor
                                   : Eigen::ColMajor>;
    /** Read-only zero-copy coordinate map type. */
    using const_map_type = Eigen::Map<matrix_type const>;
    /** Mutable zero-copy coordinate map type. */
    using map_type = Eigen::Map<matrix_type>;

    /** @brief Return the number of live point columns. */
    [[nodiscard]] size_type size() const noexcept
    {
        return static_cast<size_type>(
            M_coordinates.size() / RealDim );
    }

    /** @brief Return true when no coordinate is stored. */
    [[nodiscard]] bool empty() const noexcept
    {
        return M_coordinates.empty();
    }

    /** @brief Remove every coordinate and external-ID lookup. */
    void clear() noexcept
    {
        M_coordinates.clear();
        M_externalIds.clear();
        M_externalToLocal.clear();
        M_denseExternalIds = true;
        M_denseExternalIdBase = 0;
    }

    /**
     * @brief Reserve columns for a future rebuild.
     * @param pointCount expected point count
     */
    void reserve( size_type pointCount )
    {
        if ( pointCount >
            std::numeric_limits<size_type>::max() / RealDim )
            throw std::length_error(
                "coordinate storage exceeds addressable size" );
        M_coordinates.reserve( pointCount * RealDim );
        M_externalIds.reserve( pointCount );
    }

    /**
     * @brief Append one point coordinate.
     * @param id external point identifier
     * @param coordinate coordinate components
     * @return stable dense local point index
     */
    [[nodiscard]] local_index_type append(
        index_type id, std::span<value_type const> coordinate )
    {
        if ( coordinate.size() != RealDim )
            throw std::invalid_argument(
                "point coordinate has an invalid ambient dimension" );
        if ( size() >= std::numeric_limits<local_index_type>::max() )
            throw std::length_error(
                "coordinate storage exceeds 32-bit local indices" );
        if ( contains( id ) )
            throw std::invalid_argument(
                "duplicate external point identifier" );

        auto const localIndex =
            static_cast<local_index_type>( size() );
        if ( localIndex == 0 )
            M_denseExternalIdBase = id;
        M_coordinates.insert(
            M_coordinates.end(), coordinate.begin(),
            coordinate.end() );

        auto const expectedDenseId =
            static_cast<index_type>(
                M_denseExternalIdBase + localIndex );
        if ( M_denseExternalIds && id != expectedDenseId )
        {
            M_denseExternalIds = false;
            M_externalIds.reserve( localIndex + 1 );
            M_externalToLocal.reserve( localIndex + 1 );
            for ( local_index_type i = 0; i < localIndex; ++i )
            {
                M_externalIds.push_back(
                    static_cast<index_type>(
                        M_denseExternalIdBase + i ) );
                M_externalToLocal.emplace(
                    static_cast<index_type>(
                        M_denseExternalIdBase + i ),
                    i );
            }
        }
        if ( !M_denseExternalIds )
        {
            M_externalIds.push_back( id );
            M_externalToLocal.emplace( id, localIndex );
        }
        return localIndex;
    }

    /**
     * @brief Return whether an external point identifier is stored.
     * @param id external point identifier
     */
    [[nodiscard]] bool contains( index_type id ) const noexcept
    {
        if ( M_denseExternalIds )
            return id >= M_denseExternalIdBase &&
                   static_cast<std::uint64_t>(
                       id - M_denseExternalIdBase ) <
                       static_cast<std::uint64_t>( size() );
        return M_externalToLocal.contains( id );
    }

    /**
     * @brief Return the dense local index of an external point identifier.
     * @param id external point identifier
     */
    [[nodiscard]] local_index_type localIndex( index_type id ) const
    {
        if ( M_denseExternalIds )
        {
            if ( !contains( id ) )
                throw std::out_of_range(
                    "unknown external point identifier" );
            return static_cast<local_index_type>(
                id - M_denseExternalIdBase );
        }
        auto const iterator = M_externalToLocal.find( id );
        if ( iterator == M_externalToLocal.end() )
            throw std::out_of_range(
                "unknown external point identifier" );
        return iterator->second;
    }

    /**
     * @brief Return an external identifier from its local index.
     * @param localIndex dense local point index
     */
    [[nodiscard]] index_type
    externalId( local_index_type localIndex ) const
    {
        if ( localIndex >= size() )
            throw std::out_of_range(
                "invalid local point index" );
        return M_denseExternalIds
                   ? static_cast<index_type>(
                         M_denseExternalIdBase + localIndex )
                   : M_externalIds[localIndex];
    }

    /** @brief Return an immutable zero-copy Eigen map of all coordinates. */
    [[nodiscard]] const_map_type coordinates() const noexcept
    {
        return { M_coordinates.data(), RealDim,
                 static_cast<Eigen::Index>( size() ) };
    }

    /** @brief Return a mutable zero-copy Eigen map of all coordinates. */
    [[nodiscard]] map_type coordinates() noexcept
    {
        return { M_coordinates.data(), RealDim,
                 static_cast<Eigen::Index>( size() ) };
    }

    /** @brief Return storage and external-lookup counters. */
    [[nodiscard]] GeoCoordinateStorageStatistics statistics() const noexcept
    {
        return {
            size(), size() * RealDim, M_externalToLocal.size()
        };
    }

  private:
    std::vector<
        value_type, Eigen::aligned_allocator<value_type>>
        M_coordinates;
    std::vector<index_type> M_externalIds;
    std::unordered_map<index_type, local_index_type>
        M_externalToLocal;
    bool M_denseExternalIds = true;
    index_type M_denseExternalIdBase = 0;
};

} // namespace Feel

#endif // FEELPP_FEELMESH_GEOCOORDINATESTORAGE_HPP
