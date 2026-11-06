/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#pragma once

#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <stdexcept>
#include <vector>

namespace Feel::detail
{
/**
 * @brief Merges compact, lexicographically sorted marker sequences.
 *
 * Each sequence stores records as `[length, marker ID...]`. Marker IDs, rather
 * than record lengths, define their ordering so fragment numbering stays
 * consistent with `std::set<Marker>`.
 */
template <typename Value>
struct MergeMarkerBuffers
{
    /**
     * @brief Produces the sorted union of two compact marker sequences.
     * @param left First sorted sequence.
     * @param right Second sorted sequence.
     * @return The deduplicated sorted union in the same compact representation.
     */
    std::vector<Value> operator()( std::vector<Value> const& left, std::vector<Value> const& right ) const
    {
        std::vector<Value> result;
        result.reserve( left.size() + right.size() );
        std::size_t i = 0, j = 0;
        while ( i < left.size() && j < right.size() )
        {
            auto a = std::span<Value const>( left ).subspan( i + 1, left[i] );
            auto b = std::span<Value const>( right ).subspan( j + 1, right[j] );
            if ( std::equal( a.begin(), a.end(), b.begin(), b.end() ) )
            {
                result.insert( result.end(), left.begin() + i, left.begin() + i + a.size() + 1 );
                i += a.size() + 1;
                j += b.size() + 1;
            }
            else if ( std::lexicographical_compare( a.begin(), a.end(), b.begin(), b.end() ) )
            {
                result.insert( result.end(), left.begin() + i, left.begin() + i + a.size() + 1 );
                i += a.size() + 1;
            }
            else
            {
                result.insert( result.end(), right.begin() + j, right.begin() + j + b.size() + 1 );
                j += b.size() + 1;
            }
        }
        result.insert( result.end(), left.begin() + i, left.end() );
        result.insert( result.end(), right.begin() + j, right.end() );
        return result;
    }
};

/**
 * @brief Builds the globally consistent marker-fragment mapping.
 *
 * Consumes local marker sets, reduces compact buffers to rank zero, and
 * broadcasts primitive arrays in bounded chunks. Only the returned public
 * mapping is replicated; no tree of sets is globally retained.
 *
 * @param localMarkers Local marker combinations; empty on return.
 * @param comm Communicator containing the participating ranks.
 * @return Fragment identifiers mapped to lexicographically ordered markers.
 * @throw std::overflow_error if the number of fragments exceeds int range.
 */
template <typename Marker>
std::map<int, Marker> globalMarkerFragments( std::set<Marker>& localMarkers, boost::mpi::communicator const& comm )
{
    std::map<int, Marker> result;
    if ( comm.size() == 1 )
    {
        while ( !localMarkers.empty() )
        {
            auto node = localMarkers.extract( localMarkers.begin() );
            result.emplace_hint( result.end(), result.size(), std::move( node.value() ) );
        }
        return result;
    }

    using value_type = typename Marker::value_type;
    std::vector<value_type> local;
    std::size_t count = localMarkers.size();
    for ( auto const& marker : localMarkers )
        count += marker.size();
    local.reserve( count );
    for ( auto const& marker : localMarkers )
    {
        local.push_back( marker.size() );
        local.insert( local.end(), marker.begin(), marker.end() );
    }
    localMarkers.clear();

    std::vector<value_type> global;
    // Reduce one buffer as an object, not Boost.MPI's element-wise vector overload.
    boost::mpi::reduce( comm, &local, 1, &global, MergeMarkerBuffers<value_type>{}, 0 );
    std::vector<value_type>().swap( local );
    std::uint64_t length = global.size();
    boost::mpi::broadcast( comm, length, 0 );
    global.resize( length );
    constexpr std::size_t chunkSize = 1024 * 1024;
    for ( std::size_t offset = 0; offset < global.size(); offset += chunkSize )
    {
        int chunk = std::min( chunkSize, global.size() - offset );
        boost::mpi::broadcast( comm, global.data() + offset, chunk, 0 );
    }

    for ( std::size_t offset = 0; offset < global.size(); )
    {
        if ( result.size() > std::size_t( std::numeric_limits<int>::max() ) )
            throw std::overflow_error( "mesh fragment ID exceeds int range" );
        auto values = std::span<value_type const>( global ).subspan( offset + 1, global[offset] );
        Marker marker;
        marker.assign( values );
        result.emplace_hint( result.end(), result.size(), std::move( marker ) );
        offset += values.size() + 1;
    }
    return result;
}
}
